#[path = "common.rs"]
mod common;
use common::*;
use vfvchem_vfvm2d::*;

/// Benchmark: composite-wall heat conduction with discontinuous conductivity.
/// Run with: `cargo run --release --example bench_vert`
///
/// Geometry (examples/obj/obj_vert.obj; L = H = 0.1 m, interface at x = a = 0.05 m)
///
///         bnd_tl     bnd_tr
///       +--------+--------+
///       |        |        |
/// bnd_l | dom_l itf_v dom_r| bnd_r
///       |        |        |
///       +--------+--------+
///         bnd_bl     bnd_br
///
/// Equations (steady, no source)
/// 0 = div( k * grad(T) ) in dom_l and dom_r
/// T = T_l at x = 0, T = T_r at x = L, zero flux on top and bottom
/// T and k * grad(T) . n continuous on itf_v
///
/// Exact solution
/// k = k_0 * (1 + alp * x / L) in each layer (alp = 0 -> constant layer conductivity)
/// R(x) = int_0^x dx / k = L / (k_0 * alp) * ln(1 + alp * x / L) (or x / k_0 if alp = 0)
/// q = (T_l - T_r) / (R_l(a) + R_r(L) - R_r(a))
/// T = T_l - q * R_l(x) in dom_l; T = T_a - q * (R_r(x) - R_r(a)) in dom_r

const LEN: f64 = 0.1;
const POS_A: f64 = 0.05;
const TEMP_L: f64 = 400.0;
const TEMP_R: f64 = 300.0;

struct Case {
    name: &'static str,
    cond_l: f64,  // k_0 of dom_l
    cond_r: f64,  // k_0 of dom_r
    alp: f64,     // linear variation of k in x
}

struct Result {
    h: f64,
    num_vert: usize,
    err: (f64, f64, f64),  // L1, L2, Linf
    pos_max: (f64, f64),   // location of the max error
    err_far: f64,          // Linf excluding vertices within 6 h of a geometric corner (chamfered corners)
    hflw_l: f64,           // heat flow into domain at bnd_l
    hflw_r: f64,           // heat flow into domain at bnd_r
    hflw_i: f64,           // heat flow from dom_l to dom_r at itf_v
    bal: f64,              // global imbalance divided by its round-off scale
    itf_mis: f64,          // max nodal interface flux mismatch divided by its round-off scale
    eq_res: f64,           // max interior equation residual divided by its round-off scale
    temp_min: f64,
    temp_max: f64,
    conv: bool,            // true if the solver reached its residual tolerance
}

fn res_func(x: f64, cond_0: f64, alp: f64) -> f64 {
    // thermal resistance from 0 to x
    if alp == 0.0 { x / cond_0 } else { LEN / (cond_0 * alp) * (1.0 + alp * x / LEN).ln() }
}

fn exact(case: &Case) -> (f64, Box<dyn Fn(f64) -> f64>) {
    let (kl, kr, alp) = (case.cond_l, case.cond_r, case.alp);
    let q = (TEMP_L - TEMP_R) / (res_func(POS_A, kl, alp) + res_func(LEN, kr, alp) - res_func(POS_A, kr, alp));
    let temp_a = TEMP_L - q * res_func(POS_A, kl, alp);
    let func = Box::new(move |x: f64| {
        if x <= POS_A { TEMP_L - q * res_func(x, kl, alp) } else { temp_a - q * (res_func(x, kr, alp) - res_func(POS_A, kr, alp)) }
    });
    (q, func)
}

fn run(case: &Case, h: f64) -> Result {
    // variables
    let mut vars = Variables::new_obj("examples/obj/obj_vert.obj", h, h, "target/bench_output/vert").unwrap();

    // geometry
    let dom_l = vars.add_dom("dom_l").unwrap();
    let dom_r = vars.add_dom("dom_r").unwrap();
    let bnd_l = vars.add_bnd(dom_l, "bnd_l").unwrap();
    let bnd_r = vars.add_bnd(dom_r, "bnd_r").unwrap();
    let bnd_bl = vars.add_bnd(dom_l, "bnd_bl").unwrap();
    let bnd_tl = vars.add_bnd(dom_l, "bnd_tl").unwrap();
    let bnd_br = vars.add_bnd(dom_r, "bnd_br").unwrap();
    let bnd_tr = vars.add_bnd(dom_r, "bnd_tr").unwrap();
    let itf_v = vars.add_itf(dom_l, dom_r, "itf_v").unwrap();

    // unknowns
    let temp_l = vars.add_scldom_unkn(dom_l, 350.0, "temp").unwrap();
    let temp_r = vars.add_scldom_unkn(dom_r, 350.0, "temp").unwrap();

    // domain values
    let (kl, kr, alp) = (case.cond_l, case.cond_r, case.alp);
    let cond_l = vars.add_scldom_func(dom_l, Box::new(move |_t, (x, _y), _u| kl * (1.0 + alp * x / LEN)), vec![], "cond").unwrap();
    let cond_r = vars.add_scldom_func(dom_r, Box::new(move |_t, (x, _y), _u| kr * (1.0 + alp * x / LEN)), vec![], "cond").unwrap();
    let hlin_l = vars.add_scldom_cnst(dom_l, 0.0, "hlin").unwrap();
    let hsrc_l = vars.add_scldom_cnst(dom_l, 0.0, "hsrc").unwrap();
    let hlin_r = vars.add_scldom_cnst(dom_r, 0.0, "hlin").unwrap();
    let hsrc_r = vars.add_scldom_cnst(dom_r, 0.0, "hsrc").unwrap();

    // boundary values
    let tset_l = vars.add_sclbnd_cnst(bnd_l, TEMP_L).unwrap();
    let tset_r = vars.add_sclbnd_cnst(bnd_r, TEMP_R).unwrap();
    let hflx_bl = vars.add_sclbnd_cnst(bnd_bl, 0.0).unwrap();
    let hflx_tl = vars.add_sclbnd_cnst(bnd_tl, 0.0).unwrap();
    let hflx_br = vars.add_sclbnd_cnst(bnd_br, 0.0).unwrap();
    let hflx_tr = vars.add_sclbnd_cnst(bnd_tr, 0.0).unwrap();

    // physics
    let mut phys = SteadyHeat::new().unwrap();
    phys.add_heat_dom(dom_l, temp_l, cond_l, hlin_l, hsrc_l).unwrap();
    phys.add_heat_dom(dom_r, temp_r, cond_r, hlin_r, hsrc_r).unwrap();
    phys.add_temp_bnd(bnd_l, tset_l).unwrap();
    phys.add_temp_bnd(bnd_r, tset_r).unwrap();
    phys.add_hflx_bnd(bnd_bl, hflx_bl).unwrap();
    phys.add_hflx_bnd(bnd_tl, hflx_tl).unwrap();
    phys.add_hflx_bnd(bnd_br, hflx_br).unwrap();
    phys.add_hflx_bnd(bnd_tr, hflx_tr).unwrap();
    phys.add_hcnt_itf(itf_v).unwrap();

    // solver; linear problem, so the first iteration is the LU solution
    // the normalized residual of the LU solution stagnates at a round-off floor that grows with k_r / k_l
    // (~1e-13 for 1, ~1e-8 for 1e3, ~1e-5 for 1e6); tol = 1e-4 accepts the first LU solution
    let mut lin = RslabLu::new(4).unwrap();
    let conv = phys.solve(&mut vars, &mut lin, 5, 1e-4, 1.0).is_ok();

    // errors
    let (_, temp_func) = exact(case);
    let mut err = Vec::new();
    let mut area = Vec::new();
    let mut temp_min = f64::INFINITY;
    let mut temp_max = f64::NEG_INFINITY;
    let mut pos_max = (0.0, 0.0);
    let mut err_max = 0.0;
    let mut err_far: f64 = 0.0;
    let crn_list = [(0.0, 0.0), (POS_A, 0.0), (LEN, 0.0), (LEN, LEN), (POS_A, LEN), (0.0, LEN)];
    for (dom_id, temp_id) in [(dom_l, temp_l), (dom_r, temp_r)] {
        let pos = vert_pos(&vars, dom_id);
        for (vid, &(x, y)) in pos.iter().enumerate() {
            let temp = vars.scldom[temp_id].value[vid];
            let e = temp - temp_func(x);
            if e.abs() > err_max { err_max = e.abs(); pos_max = (x, y); }
            let dist_crn = crn_list.iter().map(|&(xc, yc): &(f64, f64)| ((x - xc).powi(2) + (y - yc).powi(2)).sqrt()).fold(f64::INFINITY, f64::min);
            if dist_crn > 6.0 * h { err_far = err_far.max(e.abs()); }
            err.push(e);
            area.push(vars.tri[dom_id].vert_area[vid]);
            temp_min = temp_min.min(temp);
            temp_max = temp_max.max(temp);
        }
    }

    // boundary and interface heat flows from the residual of the domain operators
    let diff_l = HeatDomDiffusion::new(&vars, dom_l, temp_l, cond_l).unwrap();
    let diff_r = HeatDomDiffusion::new(&vars, dom_r, temp_r, cond_r).unwrap();
    let (res_l, scl_l) = get_res(&vars, &[&diff_l]);
    let (res_r, scl_r) = get_res(&vars, &[&diff_r]);
    let start_l = scl_start(&vars, temp_l);
    let start_r = scl_start(&vars, temp_r);
    let vid_l = bnd_vid(&vars, bnd_l);
    let vid_r = bnd_vid(&vars, bnd_r);
    let hflw_l: f64 = vid_l.iter().map(|vid| -res_l[start_l + vid]).sum();
    let hflw_r: f64 = vid_r.iter().map(|vid| -res_r[start_r + vid]).sum();
    let bal_scl: f64 = vid_l.iter().map(|vid| scl_l[start_l + vid]).sum::<f64>() + vid_r.iter().map(|vid| scl_r[start_r + vid]).sum::<f64>();
    let bal = (hflw_l + hflw_r).abs() / (EPS * bal_scl);
    let vid_il = itf_vid(&vars, itf_v, dom_l);
    let vid_ir = itf_vid(&vars, itf_v, dom_r);
    let mut hflw_i = 0.0;
    let mut itf_mis: f64 = 0.0;
    for (vl, vr) in vid_il.iter().zip(vid_ir.iter()) {
        let in_l = -res_l[start_l + vl];
        let in_r = -res_r[start_r + vr];
        hflw_i += in_r;
        itf_mis = itf_mis.max((in_l + in_r).abs() / (EPS * (scl_l[start_l + vl] + scl_r[start_r + vr])));
    }
    let mut eq_res: f64 = 0.0;
    for vid in 0..vars.dom[dom_l].num_pts() { eq_res = eq_res.max(res_l[start_l + vid].abs() / (EPS * scl_l[start_l + vid])); }
    for vid in 0..vars.dom[dom_r].num_pts() { eq_res = eq_res.max(res_r[start_r + vid].abs() / (EPS * scl_r[start_r + vid])); }

    let num_vert = vars.dom[dom_l].num_vert() + vars.dom[dom_r].num_vert();
    Result { h, num_vert, err: err_norm(&err, &area), pos_max, err_far, hflw_l, hflw_r, hflw_i, bal, itf_mis, eq_res, temp_min, temp_max, conv }
}

fn main() {
    let case_list = [
        Case { name: "ratio_1", cond_l: 10.0, cond_r: 10.0, alp: 0.0 },
        Case { name: "ratio_1e3", cond_l: 1.0, cond_r: 1000.0, alp: 0.0 },
        Case { name: "ratio_1e-3", cond_l: 1000.0, cond_r: 1.0, alp: 0.0 },
        Case { name: "ratio_1e6", cond_l: 1.0e-3, cond_r: 1.0e3, alp: 0.0 },
        Case { name: "ratio_1e3_kx", cond_l: 1.0, cond_r: 1000.0, alp: 2.0 },
    ];
    let h_list = [0.01, 0.005, 0.0025, 0.00125, 0.000625];
    let mut fail_list = Vec::new();

    for case in case_list.iter() {
        let (q, _) = exact(case);
        let hflw_ex = q * LEN;  // exact heat flow through the full wall height
        let dtemp = TEMP_L - TEMP_R;
        let mut res_list = Vec::new();
        for &h in h_list.iter() {
            res_list.push(run(case, h));
        }

        // table
        // Qerr_* = relative error of heat flow at bnd_l, bnd_r, itf_v; bal, itf_mis, eq_res in units of round-off (eps * sum |terms|)
        println!("\n=== bench_vert case {} (k_l = {}, k_r = {}, alp = {}; q_exact = {:.6e} W m-2) ===", case.name, case.cond_l, case.cond_r, case.alp, q);
        println!("{:>9} {:>6} {:>10} {:>10} {:>10} {:>6} {:>6} {:>6} {:>10} {:>10} {:>10} {:>8} {:>8} {:>8} {:>9} {:>8} {:>8} {:>5} {:>10}", "h", "N", "L1", "L2", "Linf", "o_L1", "o_L2", "o_Linf", "Qerr_l", "Qerr_r", "Qerr_itf", "bal", "itf_mis", "eq_res", "overshoot", "x_max", "y_max", "conv", "Linf_far");
        for (i, r) in res_list.iter().enumerate() {
            let (o1, o2, oi) = if i == 0 { (f64::NAN, f64::NAN, f64::NAN) } else {
                let p = &res_list[i - 1];
                (get_order(p.err.0, r.err.0, p.h, r.h), get_order(p.err.1, r.err.1, p.h, r.h), get_order(p.err.2, r.err.2, p.h, r.h))
            };
            let over = (r.temp_max - TEMP_L).max(TEMP_R - r.temp_min).max(0.0);
            println!("{:>9.3e} {:>6} {:>10.3e} {:>10.3e} {:>10.3e} {:>6.2} {:>6.2} {:>6.2} {:>10.3e} {:>10.3e} {:>10.3e} {:>8.1} {:>8.1} {:>8.1} {:>9.2e} {:>8.5} {:>8.5} {:>5} {:>10.3e}",
                r.h, r.num_vert, r.err.0, r.err.1, r.err.2, o1, o2, oi,
                (r.hflw_l - hflw_ex) / hflw_ex, (-r.hflw_r - hflw_ex) / hflw_ex, (r.hflw_i - hflw_ex) / hflw_ex,
                r.bal, r.itf_mis, r.eq_res, over, r.pos_max.0, r.pos_max.1, r.conv, r.err_far);
        }

        // pass/fail criteria
        // expected: 2nd order in L1/L2 (exact in the interior for alp = 0), >= 1st order in Linf and heat flow
        // orders are least-squares fits over all meshes because the point sampling is random on each mesh
        // conservation: interface flux mismatch, global imbalance, and interior residuals at round-off level (< 1e3 eps * sum |terms|)
        let h_vec: Vec<f64> = res_list.iter().map(|r| r.h).collect();
        let o1 = fit_order(&res_list.iter().map(|r| r.err.0).collect::<Vec<f64>>(), &h_vec);
        let o2 = fit_order(&res_list.iter().map(|r| r.err.1).collect::<Vec<f64>>(), &h_vec);
        let oi = fit_order(&res_list.iter().map(|r| r.err.2).collect::<Vec<f64>>(), &h_vec);
        let qerr: Vec<f64> = res_list.iter().map(|r| ((r.hflw_l - hflw_ex) / hflw_ex).abs()).collect();
        let oq = fit_order(&qerr, &h_vec);
        let r = &res_list[res_list.len() - 1];
        println!("fitted orders: L1 {o1:.2}, L2 {o2:.2}, Linf {oi:.2}, heat flow {oq:.2}");
        check(&format!("{} L1 order", case.name), o1 >= 1.7 || r.err.0 < 1e-10 * dtemp, format!("{o1:.2} (expected 2; min 1.7), L1 = {:.3e} K", r.err.0), &mut fail_list);
        check(&format!("{} L2 order", case.name), o2 >= 1.7 || r.err.1 < 1e-10 * dtemp, format!("{o2:.2} (expected 2; min 1.7), L2 = {:.3e} K", r.err.1), &mut fail_list);
        check(&format!("{} Linf order", case.name), oi >= 0.85 || r.err.2 < 1e-10 * dtemp, format!("{oi:.2} (expected >= 1; min 0.85), Linf = {:.3e} K", r.err.2), &mut fail_list);
        check(&format!("{} heat flow order", case.name), oq >= 0.85, format!("{oq:.2} (expected >= 1; min 0.85), finest rel err = {:.3e}", qerr[qerr.len() - 1]), &mut fail_list);
        let mis = res_list.iter().map(|r| r.itf_mis).fold(0.0, f64::max);
        let bal = res_list.iter().map(|r| r.bal).fold(0.0, f64::max);
        let eqr = res_list.iter().map(|r| r.eq_res).fold(0.0, f64::max);
        let over = res_list.iter().map(|r| (r.temp_max - TEMP_L).max(TEMP_R - r.temp_min).max(0.0)).fold(0.0, f64::max);
        let conv = res_list.iter().all(|r| r.conv);
        check(&format!("{} interface flux continuity", case.name), mis < 1e3, format!("max nodal |q_l + q_r| = {mis:.1} eps * sum |terms| (tol 1e3)"), &mut fail_list);
        check(&format!("{} global balance", case.name), bal < 1e3, format!("|Q_l + Q_r| = {bal:.1} eps * sum |terms| (tol 1e3)"), &mut fail_list);
        check(&format!("{} interior residual", case.name), eqr < 1e3, format!("max |r| = {eqr:.1} eps * sum |terms| (tol 1e3)"), &mut fail_list);
        check(&format!("{} boundedness", case.name), over < 1e3 * EPS * TEMP_L, format!("max overshoot = {over:.3e} K (tol 1e3 eps T_l)"), &mut fail_list);
        check(&format!("{} solver convergence", case.name), conv, format!("{conv}"), &mut fail_list);
    }

    println!("\nbench_vert: {} failed criteria", fail_list.len());
    for f in fail_list.iter() { println!("  FAIL {f}"); }
    if !fail_list.is_empty() { std::process::exit(1); }
}
