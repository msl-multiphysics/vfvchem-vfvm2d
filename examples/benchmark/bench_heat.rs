#[path = "common.rs"]
mod common;
use common::*;
use vfvchem_vfvm2d::*;

/// Benchmark: nonlinear heat conduction with k(T) = k0 * (1 + beta * T).
/// Run with: `cargo run --release --example bench_heat`
///
/// Geometry: square 0 <= x, y <= L (L = 0.1 m) from Variables::new_bounds
///
/// Equations (steady, no source)
/// 0 = div( k(T) * grad(T) ); k(T) = k0 * (1 + beta * T)
///
/// Kirchhoff transform
/// phi(T) = int_0^T k / k0 dT = T + beta * T^2 / 2
/// k0 * lap(phi) = 0, so phi is any harmonic function
/// T = (sqrt(1 + 2 * beta * phi) - 1) / beta
///
/// Case "slab" (1D): T = T_hot at x = 0, T = 0 at x = L, zero flux on top and bottom
/// phi = phi_hot * (1 - x / L); heat flux q = k0 * phi_hot / L
///
/// Case "harm" (2D): T from phi on all four boundaries (Dirichlet)
/// phi = phi_hot * (0.25 * (s^2 - t^2) + 0.25 * s + 0.5); s = x / L, t = y / L
/// inward heat flow per boundary = k0 * int grad(phi) . n dl
/// left = -0.25 k0 phi_hot, right = 0.75 k0 phi_hot, bottom = 0, top = -0.5 k0 phi_hot
///
/// Property ratios k(T_hot) / k(0) = 1 + 100 beta = 2, 10, 100 (T_hot = 100 K above reference)

const LEN: f64 = 0.1;
const TEMP_HOT: f64 = 100.0;
const COND_0: f64 = 1.0;

struct Result {
    h: f64,
    num_vert: usize,
    err: (f64, f64, f64),  // L1, L2, Linf
    err_far: f64,          // Linf excluding vertices within 6 h of a corner
    hflw: [f64; 4],        // inward heat flow at bnd_l, bnd_r, bnd_b, bnd_t
    bal: f64,              // global imbalance in units of round-off
    eq_res: f64,           // max interior residual in units of round-off
    temp_min: f64,
    temp_max: f64,
    conv: bool,
    temp_vec: Vec<f64>,    // solution for the iterative-error comparison
}

fn phi_func(case: &str, phi_hot: f64, x: f64, y: f64) -> f64 {
    let (s, t) = (x / LEN, y / LEN);
    if case == "slab" { phi_hot * (1.0 - s) } else { phi_hot * (0.25 * (s * s - t * t) + 0.25 * s + 0.5) }
}

fn temp_func(beta: f64, phi: f64) -> f64 {
    ((1.0 + 2.0 * beta * phi).sqrt() - 1.0) / beta
}

fn run(case: &'static str, beta: f64, h: f64, tol: f64) -> Result {
    let phi_hot = TEMP_HOT + 0.5 * beta * TEMP_HOT * TEMP_HOT;

    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, LEN, LEN, h, h, "target/bench_output/heat").unwrap();

    // geometry
    let dom = vars.add_dom("dom").unwrap();
    let bnd_l = vars.add_bnd(dom, "bnd_l").unwrap();
    let bnd_r = vars.add_bnd(dom, "bnd_r").unwrap();
    let bnd_b = vars.add_bnd(dom, "bnd_b").unwrap();
    let bnd_t = vars.add_bnd(dom, "bnd_t").unwrap();

    // unknowns
    let temp = vars.add_scldom_unkn(dom, 0.5 * TEMP_HOT, "temp").unwrap();

    // domain values; conductivity depends on the temperature unknown
    let cond = vars.add_scldom_func(dom, Box::new(move |_t, _pos, u| COND_0 * (1.0 + beta * u[0])), vec![temp], "cond").unwrap();
    let hlin = vars.add_scldom_cnst(dom, 0.0, "hlin").unwrap();
    let hsrc = vars.add_scldom_cnst(dom, 0.0, "hsrc").unwrap();

    // physics
    let mut phys = SteadyHeat::new().unwrap();
    phys.add_heat_dom(dom, temp, cond, hlin, hsrc).unwrap();
    let mut vid_bnd = Vec::new();
    for (i, &bnd) in [bnd_l, bnd_r, bnd_b, bnd_t].iter().enumerate() {
        vid_bnd.push(bnd_vid(&vars, bnd));
        if case == "slab" && i >= 2 {
            let hflx = vars.add_sclbnd_cnst(bnd, 0.0).unwrap();
            phys.add_hflx_bnd(bnd, hflx).unwrap();
        } else {
            let tset = vars.add_sclbnd_func(bnd, Box::new(move |_t, (x, y), _u| temp_func(beta, phi_func(case, phi_hot, x, y))), vec![]).unwrap();
            phys.add_temp_bnd(bnd, tset).unwrap();
        }
    }

    // solver; picard iteration on k(T)
    let mut lin = RslabLu::new(4).unwrap();
    let conv = phys.solve(&mut vars, &mut lin, 300, tol, 1.0).is_ok();

    // errors
    let pos = vert_pos(&vars, dom);
    let mut err = Vec::new();
    let mut err_far: f64 = 0.0;
    for (vid, &(x, y)) in pos.iter().enumerate() {
        let e = vars.scldom[temp].value[vid] - temp_func(beta, phi_func(case, phi_hot, x, y));
        let dist_crn = [(0.0, 0.0), (LEN, 0.0), (LEN, LEN), (0.0, LEN)].iter().map(|&(xc, yc): &(f64, f64)| ((x - xc).powi(2) + (y - yc).powi(2)).sqrt()).fold(f64::INFINITY, f64::min);
        if dist_crn > 6.0 * h { err_far = err_far.max(e.abs()); }
        err.push(e);
    }
    let temp_vec = vars.scldom[temp].value.clone();
    let temp_min = temp_vec.iter().cloned().fold(f64::INFINITY, f64::min);
    let temp_max = temp_vec.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // boundary heat flows from the residual of the domain operator
    let diff = HeatDomDiffusion::new(&vars, dom, temp, cond).unwrap();
    let (res, scl) = get_res(&vars, &[&diff]);
    let start = scl_start(&vars, temp);
    let mut hflw = [0.0; 4];
    let mut bal_sum = 0.0;
    let mut bal_scl = 0.0;
    for i in 0..4 {
        hflw[i] = vid_bnd[i].iter().map(|vid| -res[start + vid]).sum();
        bal_sum += hflw[i];
        bal_scl += vid_bnd[i].iter().map(|vid| scl[start + vid]).sum::<f64>();
    }
    let bal = bal_sum.abs() / (EPS * bal_scl);
    let mut eq_res: f64 = 0.0;
    for vid in 0..vars.dom[dom].num_pts() { eq_res = eq_res.max(res[start + vid].abs() / (EPS * scl[start + vid])); }

    Result { h, num_vert: pos.len(), err: err_norm(&err, &vars.tri[dom].vert_area), err_far, hflw, bal, eq_res, temp_min, temp_max, conv, temp_vec }
}

fn main() {
    let h_list = [0.01, 0.005, 0.0025, 0.00125, 0.000625, 0.0003125];
    let tol_study = 1e-10;  // picard tolerance for the convergence study
    let tol_ref = 1e-12;    // tighter tolerance to estimate the iterative error of the study solution
    let tol_dflt = 1e-3;    // tolerance used in the examples (information only)
    let mut fail_list = Vec::new();

    for (case, beta_list) in [("slab", vec![0.01, 0.03, 0.09, 0.99]), ("harm", vec![0.01, 0.09, 0.99])] {
    for beta in beta_list {
        let phi_hot = TEMP_HOT + 0.5 * beta * TEMP_HOT * TEMP_HOT;
        let hflw_ex = if case == "slab" {
            [COND_0 * phi_hot, -COND_0 * phi_hot, 0.0, 0.0]
        } else {
            [-0.25 * COND_0 * phi_hot, 0.75 * COND_0 * phi_hot, 0.0, -0.5 * COND_0 * phi_hot]
        };
        let hflw_ref = 0.75 * COND_0 * phi_hot;  // scale for relative heat flow errors
        let tmin_ex = if case == "slab" { 0.0 } else { temp_func(beta, 0.25 * phi_hot) };

        // resolution length of the exact solution
        // (k / k0)^2 = 1 + 2 beta phi is linear in phi; T has a curvature length dl = L * kmin^2 / (kmax^2 - kmin^2)
        // the slab with k ratio r has dl = L / (r^2 - 1) at the cold wall; orders are only expected for h < dl
        let k2_min = 1.0 + 2.0 * beta * if case == "slab" { 0.0 } else { 0.25 * phi_hot };
        let k2_max = 1.0 + 2.0 * beta * phi_hot;
        let dl = LEN * k2_min / (k2_max - k2_min);

        let mut res_list = Vec::new();
        let mut itr_list = Vec::new();  // (L2 of T_study - T_ref, L2 of T_dflt - T_ref)
        for &h in h_list.iter() {
            let r = run(case, beta, h, tol_study);
            let r_ref = run(case, beta, h, tol_ref);
            let r_dflt = run(case, beta, h, tol_dflt);
            let area = vec![1.0; r.temp_vec.len()];
            let dif_s: Vec<f64> = r.temp_vec.iter().zip(r_ref.temp_vec.iter()).map(|(a, b)| a - b).collect();
            let dif_d: Vec<f64> = r_dflt.temp_vec.iter().zip(r_ref.temp_vec.iter()).map(|(a, b)| a - b).collect();
            itr_list.push((err_norm(&dif_s, &area).1, err_norm(&dif_d, &area).1, r_ref.conv));
            res_list.push(r);
        }

        // table
        // Qerr_* = (Q - Q_exact) / (0.75 k0 phi_hot) per boundary; bal = |sum Q| / (0.75 k0 phi_hot)
        // it_L2 = L2 of T(tol 1e-10) - T(tol 1e-12); it_dflt = L2 of T(tol 1e-3) - T(tol 1e-12)
        println!("\n=== bench_heat case {case} beta = {beta} (k ratio {:.0}; resolution length {:.3e} m) ===", 1.0 + beta * TEMP_HOT, dl);
        println!("{:>9} {:>6} {:>10} {:>10} {:>10} {:>6} {:>6} {:>6} {:>10} {:>10} {:>10} {:>10} {:>10} {:>9} {:>9} {:>9} {:>9} {:>5}", "h", "N", "L1", "L2", "Linf", "o_L1", "o_L2", "o_Linf", "Linf_far", "Qerr_l", "Qerr_r", "Qerr_b", "Qerr_t", "bal", "overshoot", "it_L2", "it_dflt", "conv");
        for (i, r) in res_list.iter().enumerate() {
            let (o1, o2, oi) = if i == 0 { (f64::NAN, f64::NAN, f64::NAN) } else {
                let p = &res_list[i - 1];
                (get_order(p.err.0, r.err.0, p.h, r.h), get_order(p.err.1, r.err.1, p.h, r.h), get_order(p.err.2, r.err.2, p.h, r.h))
            };
            let over = (r.temp_max - TEMP_HOT).max(tmin_ex - r.temp_min).max(0.0);
            let bal = r.hflw.iter().sum::<f64>().abs() / hflw_ref;
            println!("{:>9.3e} {:>6} {:>10.3e} {:>10.3e} {:>10.3e} {:>6.2} {:>6.2} {:>6.2} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>5}",
                r.h, r.num_vert, r.err.0, r.err.1, r.err.2, o1, o2, oi, r.err_far,
                (r.hflw[0] - hflw_ex[0]) / hflw_ref, (r.hflw[1] - hflw_ex[1]) / hflw_ref, (r.hflw[2] - hflw_ex[2]) / hflw_ref, (r.hflw[3] - hflw_ex[3]) / hflw_ref,
                bal, over, itr_list[i].0, itr_list[i].1, r.conv && itr_list[i].2);
        }

        // pass/fail criteria
        // expected: 2nd order in L1/L2, >= 1st order in Linf and boundary heat flow, fitted over meshes with h < dl
        // if fewer than two meshes resolve dl, only a monotonic error decrease is required (pre-asymptotic)
        // global balance |sum Q| / Q < 1e-8; internal faces cancel exactly, so the imbalance is the unconverged residual
        // iterative error of the study solution below 1 % of the L2 discretization error on every mesh
        let name = format!("{case}_beta{beta}");
        let sel: Vec<usize> = (0..res_list.len()).filter(|&i| res_list[i].h < dl).collect();
        let r = &res_list[res_list.len() - 1];
        let qerr: Vec<f64> = res_list.iter().map(|r| (0..4).map(|i| (r.hflw[i] - hflw_ex[i]).abs()).fold(0.0, f64::max) / hflw_ref).collect();
        if sel.len() >= 2 {
            let h_vec: Vec<f64> = sel.iter().map(|&i| res_list[i].h).collect();
            let o1 = fit_order(&sel.iter().map(|&i| res_list[i].err.0).collect::<Vec<f64>>(), &h_vec);
            let o2 = fit_order(&sel.iter().map(|&i| res_list[i].err.1).collect::<Vec<f64>>(), &h_vec);
            let oi = fit_order(&sel.iter().map(|&i| res_list[i].err.2).collect::<Vec<f64>>(), &h_vec);
            let oq = fit_order(&sel.iter().map(|&i| qerr[i]).collect::<Vec<f64>>(), &h_vec);
            println!("fitted orders over {} meshes with h < {dl:.3e}: L1 {o1:.2}, L2 {o2:.2}, Linf {oi:.2}, max boundary heat flow error {oq:.2}", sel.len());
            check(&format!("{name} L1 order"), o1 >= 1.7, format!("{o1:.2} (expected 2; min 1.7), L1 = {:.3e} K", r.err.0), &mut fail_list);
            check(&format!("{name} L2 order"), o2 >= 1.7, format!("{o2:.2} (expected 2; min 1.7), L2 = {:.3e} K", r.err.1), &mut fail_list);
            check(&format!("{name} Linf order"), oi >= 0.85, format!("{oi:.2} (expected >= 1; min 0.85), Linf = {:.3e} K", r.err.2), &mut fail_list);
            check(&format!("{name} heat flow order"), oq >= 0.85, format!("{oq:.2} (expected >= 1; min 0.85), finest max rel err = {:.3e}", qerr[qerr.len() - 1]), &mut fail_list);
        } else {
            let mono = (1..res_list.len()).all(|i| res_list[i].err.1 < res_list[i - 1].err.1);
            println!("pre-asymptotic: no mesh resolves the solution length {dl:.3e} m (finest h = {:.3e} m); order criteria not applicable", r.h);
            check(&format!("{name} monotonic L2 decrease (pre-asymptotic)"), mono, format!("L2 = {:.3e} K on the finest mesh", r.err.1), &mut fail_list);
        }
        let bal = res_list.iter().map(|r| r.hflw.iter().sum::<f64>().abs() / hflw_ref).fold(0.0, f64::max);
        let over = res_list.iter().map(|r| (r.temp_max - TEMP_HOT).max(tmin_ex - r.temp_min).max(0.0)).fold(0.0, f64::max);
        let conv = res_list.iter().all(|r| r.conv);
        let itr_ok = res_list.iter().zip(itr_list.iter()).all(|(r, it)| it.0 < 0.01 * r.err.1);
        check(&format!("{name} global balance"), bal < 1e-8, format!("max |sum Q| / Q = {bal:.3e} (tol 1e-8)"), &mut fail_list);
        check(&format!("{name} boundedness"), over < 1e3 * EPS * TEMP_HOT, format!("max overshoot = {over:.3e} K (tol 1e3 eps T_hot)"), &mut fail_list);
        check(&format!("{name} solver convergence"), conv, format!("picard reached tol {tol_study:e} on all meshes: {conv}"), &mut fail_list);
        check(&format!("{name} iterative vs discretization error"), itr_ok, format!("L2(T_tol1e-10 - T_tol1e-12) < 0.01 L2 error on all meshes: {itr_ok}"), &mut fail_list);
    }}

    println!("\nbench_heat: {} failed criteria", fail_list.len());
    for f in fail_list.iter() { println!("  FAIL {f}"); }
    if !fail_list.is_empty() { std::process::exit(1); }
}
