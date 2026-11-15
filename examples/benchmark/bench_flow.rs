#[path = "common.rs"]
mod common;
use common::*;
use vfvchem_vfvm2d::*;

/// Benchmark: planar Poiseuille flow with viscosity mu(y) and optionally density rho(y).
/// Run with: `cargo run --release --example bench_flow`
///
/// Geometry: channel 0 <= x <= L, 0 <= y <= H (L = 0.1 m, H = 0.05 m) from Variables::new_bounds
/// bnd_l - inlet, bnd_r - outlet, bnd_b and bnd_t - walls
///
/// Equations (steady)
/// rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) - 2/3 mu div(vel) I )
/// div(rho * vel) = 0
///
/// Fully developed solution with vel = (u(y), 0) and p = p_0 - G x
/// advection and div(rho * vel) vanish for any rho(y); the y-momentum gives dp/dy = 0
/// d/dy( mu * du/dy ) = -G with u(0) = u(H) = 0
/// mu * du/dy = C - G y
/// mu = mu_0 * (1 + a y / H); xi = 1 + a y / H
/// u = H / (a mu_0) * ( (C + G H / a) ln(xi) - (G H / a) (xi - 1) )
/// C = G H / ln(1 + a) - G H / a (wall shear mu du/dy at y = 0; C - G H at y = H)
/// Q = int_0^H u dy
///
/// Cases
/// dir - exact velocity on all four boundaries; pressure reference at one node
/// out - exact velocity at the inlet, no slip at the walls, p = 0 at the outlet (natural outflow)
/// each with rho = rho_0 (constant) and rho = rho_0 * (1 + b y / H)

const LEN: f64 = 0.1;
const HGT: f64 = 0.05;
const GRAD: f64 = 8.0;      // G = -dp/dx (Pa m-1)
const VISC_0: f64 = 0.1;    // mu_0 (Pa s)
const VISC_A: f64 = 3.0;    // a; mu(H) / mu(0) = 4
const DEN_0: f64 = 1000.0;  // rho_0 (kg m-3)

fn shear_0() -> f64 {
    GRAD * HGT / (1.0 + VISC_A).ln() - GRAD * HGT / VISC_A
}

fn vel_exact(y: f64) -> f64 {
    let xi = 1.0 + VISC_A * y / HGT;
    let gha = GRAD * HGT / VISC_A;
    HGT / (VISC_A * VISC_0) * ((shear_0() + gha) * xi.ln() - gha * (xi - 1.0))
}

struct Result {
    h: f64,
    num_vert: usize,
    err_u: (f64, f64, f64),  // L1, L2, Linf of |vel - vel_exact|
    err_p: (f64, f64, f64),  // L1, L2, Linf of p - p_exact with the mean difference removed
    vmax: f64,               // max |v|
    umax: f64,               // max u
    umin: f64,               // min u
    dpdx: f64,               // least-squares pressure gradient
    dpdy: f64,
    flow: f64,               // area-averaged u times H
    mflw_in: f64,            // inlet mass flow from the boundary rule
    mflw_out: f64,           // outlet mass flow from the boundary rule
    tau_b: (f64, f64),       // bottom wall shear (mean, max rel err) from the momentum residual, middle half of the wall
    tau_t: (f64, f64),       // top wall shear
    mass_res: f64,           // max continuity residual at non-dirichlet rows / (rho U h)
    conv: bool,
    vel_vec: Vec<(f64, f64)>,
}

fn run(case: &str, den_b: f64, h: f64, tol: f64) -> Result {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, LEN, HGT, h, h, "target/bench_output/flow").unwrap();

    // geometry
    let dom = vars.add_dom("dom").unwrap();
    let bnd_l = vars.add_bnd(dom, "bnd_l").unwrap();
    let bnd_r = vars.add_bnd(dom, "bnd_r").unwrap();
    let bnd_b = vars.add_bnd(dom, "bnd_b").unwrap();
    let bnd_t = vars.add_bnd(dom, "bnd_t").unwrap();

    // unknowns
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "vel").unwrap();
    let pres = vars.add_scldom_unkn(dom, 0.0, "pres").unwrap();

    // domain values
    let den = vars.add_scldom_func(dom, Box::new(move |_t, (_x, y), _u| DEN_0 * (1.0 + den_b * y / HGT)), vec![], "den").unwrap();
    let visc = vars.add_scldom_func(dom, Box::new(|_t, (_x, y), _u| VISC_0 * (1.0 + VISC_A * y / HGT)), vec![], "visc").unwrap();
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "fce").unwrap();

    // boundary values
    let vel_l = vars.add_vecbnd_func(bnd_l, Box::new(|_t, (_x, y), _u| (vel_exact(y), 0.0)), vec![]).unwrap();
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0)).unwrap();
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0)).unwrap();

    // physics
    let mut phys = SteadyFlow::new().unwrap();
    phys.add_flow_dom(dom, vel, pres, den, visc, fce).unwrap();
    phys.add_vel_bnd(bnd_l, vel_l).unwrap();
    phys.add_vel_bnd(bnd_b, vel_b).unwrap();
    phys.add_vel_bnd(bnd_t, vel_t).unwrap();
    if case == "dir" {
        let vel_r = vars.add_vecbnd_func(bnd_r, Box::new(|_t, (_x, y), _u| (vel_exact(y), 0.0)), vec![]).unwrap();
        phys.add_vel_bnd(bnd_r, vel_r).unwrap();
        let (x0, _) = vert_pos(&vars, dom)[vars.dom[dom].num_pts()];  // reference row is boundary node 0
        phys.set_pres_ref(dom, 0, GRAD * (LEN - x0)).unwrap();
    } else {
        let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0).unwrap();
        phys.add_pres_bnd(bnd_r, pres_r).unwrap();
    }

    // solver
    let mut lin = RslabLu::new(4).unwrap();
    let conv = phys.solve(&mut vars, &mut lin, 200, tol, 1.0).is_ok();

    // velocity and pressure errors
    let pos = vert_pos(&vars, dom);
    let area = vars.tri[dom].vert_area.clone();
    let area_sum: f64 = area.iter().sum();
    let vel_vec = vars.vecdom[vel].value.clone();
    let mut err_u = Vec::new();
    let mut err_p = Vec::new();
    let (mut vmax, mut umax, mut umin): (f64, f64, f64) = (0.0, f64::NEG_INFINITY, f64::INFINITY);
    let mut flow = 0.0;
    for (vid, &(x, y)) in pos.iter().enumerate() {
        let (u, v) = vel_vec[vid];
        err_u.push(((u - vel_exact(y)).powi(2) + v * v).sqrt());
        err_p.push(vars.scldom[pres].value[vid] - GRAD * (LEN - x));
        vmax = vmax.max(v.abs());
        umax = umax.max(u);
        umin = umin.min(u);
        flow += u * area[vid];
    }
    let flow = flow / area_sum * HGT;
    let mean_p: f64 = err_p.iter().zip(area.iter()).map(|(e, a)| e * a).sum::<f64>() / area_sum;
    let err_p: Vec<f64> = err_p.iter().map(|e| e - mean_p).collect();

    // least-squares pressure plane p = c0 + c1 x + c2 y
    let mut m = [[0.0; 3]; 3];
    let mut r = [0.0; 3];
    for (vid, &(x, y)) in pos.iter().enumerate() {
        let f = [1.0, x, y];
        for i in 0..3 { for j in 0..3 { m[i][j] += f[i] * f[j]; } r[i] += f[i] * vars.scldom[pres].value[vid]; }
    }
    let det = |a: [[f64; 3]; 3]| a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    let d = det(m);
    let mut c = [0.0; 3];
    for k in 0..3 {
        let mut mk = m;
        for i in 0..3 { mk[i][k] = r[i]; }
        c[k] = det(mk) / d;
    }

    // boundary mass flows from the 0.75 / 0.25 half-edge rule (as in FlowBndDivergence)
    let tri = &vars.tri[dom];
    let den_dom = &vars.scldom[den];
    let mut mflw = [0.0; 2];
    for (i, &bnd) in [bnd_l, bnd_r].iter().enumerate() {
        for vid in bnd_vid(&vars, bnd) {
            for (blc, &vid_n) in tri.bnd_vid[vid].iter().enumerate() {
                let (nx, ny) = tri.bnd_norm[vid][blc];
                let (u_v, v_v) = vel_vec[vid];
                let (u_n, v_n) = vel_vec[vid_n];
                let fv = den_dom.value[vid] * (u_v * nx + v_v * ny);
                let fn_ = den_dom.value[vid_n] * (u_n * nx + v_n * ny);
                mflw[i] += (0.75 * fv + 0.25 * fn_) * tri.bnd_half[vid][blc];
            }
        }
    }

    // wall shear from the x-momentum residual of the domain operators at wall nodes in the middle half of the wall
    // residual = -(wall force on the fluid); bottom: force = -mu du/dy, top: force = +mu du/dy (per unit length)
    let oper_adv = FlowDomAdvection::new(&vars, dom, vel, den).unwrap();
    let oper_grad = FlowDomGradient::new(&vars, dom, vel, pres).unwrap();
    let oper_diff = FlowDomDiffusion::new(&vars, dom, vel, visc).unwrap();
    let oper_dtrn = FlowDomTranspose::new(&vars, dom, vel, visc).unwrap();
    let oper_ddil = FlowDomDilatatory::new(&vars, dom, vel, visc).unwrap();
    let oper_src = VecDomSource::new(&vars, dom, vel, fce).unwrap();
    let (res_m, _) = get_res(&vars, &[&oper_adv, &oper_grad, &oper_diff, &oper_dtrn, &oper_ddil, &oper_src]);
    let start_v = vec_start(&vars, vel);
    let num_vert = vars.dom[dom].num_vert();
    let mut tau = [(0.0, 0.0); 2];
    for (i, &(bnd, sign, tau_ex)) in [(bnd_b, 1.0, shear_0()), (bnd_t, -1.0, shear_0() - GRAD * HGT)].iter().enumerate() {
        let mut sum = 0.0;
        let mut cnt = 0.0;
        let mut emax: f64 = 0.0;
        for vid in bnd_vid(&vars, bnd) {
            let (x, _) = pos[vid];
            if x < 0.25 * LEN || x > 0.75 * LEN { continue; }
            let len: f64 = tri.bnd_half[vid].iter().sum();
            let tau_v = sign * res_m[start_v + vid] / len;
            sum += tau_v;
            cnt += 1.0;
            emax = emax.max(((tau_v - tau_ex) / tau_ex).abs());
        }
        tau[i] = (sum / cnt, emax);
    }

    // continuity residual at rows that keep the continuity equation
    let oper_div = FlowDomDivergence::new(&vars, dom, pres, vel, den).unwrap();
    let oper_rc = FlowDomRhieChow::new(&vars, dom, pres, vel, den).unwrap();
    let mut oper_list: Vec<Box<dyn OperatorBase>> = vec![Box::new(oper_adv), Box::new(oper_grad), Box::new(oper_diff), Box::new(oper_dtrn), Box::new(oper_ddil), Box::new(oper_src)];
    let mut bnd_list = vec![bnd_l, bnd_b, bnd_t];
    if case == "dir" { bnd_list.push(bnd_r); }
    for &bnd in bnd_list.iter() {
        oper_list.push(Box::new(FlowBndDivergence::new(&vars, bnd, vel, den, pres).unwrap()));
    }
    if case == "out" {
        oper_list.push(Box::new(FlowBndTranspose::new(&vars, bnd_r, vel, visc).unwrap()));
        oper_list.push(Box::new(FlowBndDilatatory::new(&vars, bnd_r, vel, visc).unwrap()));
    }
    oper_list.push(Box::new(oper_div));
    oper_list.push(Box::new(oper_rc));  // needs the momentum diagonal assembled before it
    let oper_ref: Vec<&dyn OperatorBase> = oper_list.iter().map(|o| o.as_ref()).collect();
    let (res_c, _) = get_res(&vars, &oper_ref);
    let start_p = scl_start(&vars, pres);
    let vid_r = bnd_vid(&vars, bnd_r);
    let mut mass_res: f64 = 0.0;
    for vid in 0..num_vert {
        if case == "out" && vid_r.contains(&vid) { continue; }  // pressure dirichlet rows
        if case == "dir" && vid == vars.dom[dom].num_pts() { continue; }  // reference row
        mass_res = mass_res.max(res_c[start_p + vid].abs());
    }
    let umax_ex = vel_exact(0.5 * HGT);
    let mass_res = mass_res / (DEN_0 * umax_ex * h);

    Result { h, num_vert, err_u: err_norm(&err_u, &area), err_p: err_norm(&err_p, &area), vmax, umax, umin, dpdx: c[1], dpdy: c[2], flow, mflw_in: -mflw[0], mflw_out: mflw[1], tau_b: tau[0], tau_t: tau[1], mass_res, conv, vel_vec }
}

fn main() {
    let h_list = [0.005, 0.0025, 0.00125, 0.000625, 0.0003125];
    let tol_study = 1e-9;
    let tol_ref = 1e-11;
    let mut fail_list = Vec::new();

    // exact integrated quantities
    let num_q = 20000;
    let flow_ex: f64 = (0..num_q).map(|i| vel_exact((i as f64 + 0.5) * HGT / num_q as f64)).sum::<f64>() * HGT / num_q as f64;
    let umax_ex = (0..num_q).map(|i| vel_exact(i as f64 * HGT / num_q as f64)).fold(0.0, f64::max);
    let tau_b_ex = shear_0();
    let tau_t_ex = shear_0() - GRAD * HGT;
    println!("exact: Q = {flow_ex:.6e} m2 s-1, u_max = {umax_ex:.6e} m s-1, tau_b = {tau_b_ex:.6e} Pa, tau_t = {tau_t_ex:.6e} Pa, dp/dx = {:.3e} Pa m-1", -GRAD);

    for case in ["dir", "out"] {
    for den_b in [0.0, 1.0] {
        let mut res_list = Vec::new();
        let mut itr_list = Vec::new();
        for &h in h_list.iter() {
            let r = run(case, den_b, h, tol_study);
            let r_ref = run(case, den_b, h, tol_ref);
            let dif: Vec<f64> = r.vel_vec.iter().zip(r_ref.vel_vec.iter()).map(|(a, b)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()).collect();
            itr_list.push((err_norm(&dif, &vec![1.0; dif.len()]).1, r_ref.conv));
            res_list.push(r);
        }

        // table; velocity errors relative to u_max, pressure errors relative to G L
        println!("\n=== bench_flow case {case} (rho = rho_0 * (1 + {den_b} y / H); mu = mu_0 * (1 + {VISC_A} y / H)) ===");
        println!("{:>9} {:>6} {:>9} {:>9} {:>9} {:>5} {:>5} {:>5} {:>9} {:>9} {:>9} {:>5} {:>5} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>5}",
            "h", "N", "u_L1", "u_L2", "u_Linf", "o_L1", "o_L2", "o_inf", "p_L1", "p_L2", "p_Linf", "op_L2", "opinf", "dpdx_err", "dpdy/G", "Q_err", "min_err", "mout_err", "taub_err", "taut_err", "vmax/U", "umax_err", "mass_res", "it_L2", "conv");
        for (i, r) in res_list.iter().enumerate() {
            let (o1, o2, oi, op2, opi) = if i == 0 { (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN) } else {
                let p = &res_list[i - 1];
                (get_order(p.err_u.0, r.err_u.0, p.h, r.h), get_order(p.err_u.1, r.err_u.1, p.h, r.h), get_order(p.err_u.2, r.err_u.2, p.h, r.h), get_order(p.err_p.1, r.err_p.1, p.h, r.h), get_order(p.err_p.2, r.err_p.2, p.h, r.h))
            };
            let mflw_ex = DEN_0 * flow_ex + DEN_0 * den_b / HGT * (0..num_q).map(|j| { let y = (j as f64 + 0.5) * HGT / num_q as f64; y * vel_exact(y) }).sum::<f64>() * HGT / num_q as f64;
            println!("{:>9.3e} {:>6} {:>9.2e} {:>9.2e} {:>9.2e} {:>5.2} {:>5.2} {:>5.2} {:>9.2e} {:>9.2e} {:>9.2e} {:>5.2} {:>5.2} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>5}",
                r.h, r.num_vert, r.err_u.0 / umax_ex, r.err_u.1 / umax_ex, r.err_u.2 / umax_ex, o1, o2, oi,
                r.err_p.0 / (GRAD * LEN), r.err_p.1 / (GRAD * LEN), r.err_p.2 / (GRAD * LEN), op2, opi,
                (r.dpdx + GRAD) / GRAD, r.dpdy / GRAD, (r.flow - flow_ex) / flow_ex, (r.mflw_in - mflw_ex) / mflw_ex, (r.mflw_out - mflw_ex) / mflw_ex,
                (r.tau_b.0 - tau_b_ex) / tau_b_ex, (r.tau_t.0 - tau_t_ex) / tau_t_ex, r.vmax / umax_ex, (r.umax - umax_ex) / umax_ex, r.mass_res, itr_list[i].0 / umax_ex, r.conv && itr_list[i].1);
        }
        let _ = res_list.iter().map(|r| r.umin).fold(0.0, f64::min);

        // pass/fail criteria
        // expected: velocity 2nd order in L1/L2 and >= 1st in Linf; pressure >= 1st order in L2
        // dp/dx, flow rate, and mid-wall shear errors >= 1st order; continuity residual and iterative error small
        let name = format!("{case}_rho{den_b}");
        let h_vec: Vec<f64> = res_list.iter().map(|r| r.h).collect();
        let fit = |f: &dyn Fn(&Result) -> f64| fit_order(&res_list.iter().map(|r| f(r).abs().max(1e-300)).collect::<Vec<f64>>(), &h_vec);
        let o1 = fit(&|r| r.err_u.0);
        let o2 = fit(&|r| r.err_u.1);
        let oi = fit(&|r| r.err_u.2);
        let op = fit(&|r| r.err_p.1);
        let og = fit(&|r| r.dpdx + GRAD);
        let oq = fit(&|r| r.flow - flow_ex);
        let otb = fit(&|r| r.tau_b.0 - tau_b_ex);
        let ott = fit(&|r| r.tau_t.0 - tau_t_ex);
        println!("fitted orders: u L1 {o1:.2}, u L2 {o2:.2}, u Linf {oi:.2}, p L2 {op:.2}, dp/dx {og:.2}, Q {oq:.2}, tau_b {otb:.2}, tau_t {ott:.2}");
        let r = &res_list[res_list.len() - 1];
        check(&format!("{name} velocity L1 order"), o1 >= 1.7, format!("{o1:.2} (expected 2; min 1.7)"), &mut fail_list);
        check(&format!("{name} velocity L2 order"), o2 >= 1.7, format!("{o2:.2} (expected 2; min 1.7)"), &mut fail_list);
        check(&format!("{name} velocity Linf order"), oi >= 0.85, format!("{oi:.2} (expected >= 1; min 0.85)"), &mut fail_list);
        check(&format!("{name} pressure L2 order"), op >= 0.85, format!("{op:.2} (expected >= 1; min 0.85)"), &mut fail_list);
        check(&format!("{name} pressure gradient order"), og >= 0.85, format!("{og:.2} (expected >= 1; min 0.85), finest rel err {:.2e}", (r.dpdx + GRAD) / GRAD), &mut fail_list);
        check(&format!("{name} flow rate order"), oq >= 0.85, format!("{oq:.2} (expected >= 1; min 0.85), finest rel err {:.2e}", (r.flow - flow_ex) / flow_ex), &mut fail_list);
        check(&format!("{name} bottom wall shear order"), otb >= 0.85, format!("{otb:.2} (expected >= 1; min 0.85), finest rel err {:.2e}", (r.tau_b.0 - tau_b_ex) / tau_b_ex), &mut fail_list);
        check(&format!("{name} top wall shear order"), ott >= 0.85, format!("{ott:.2} (expected >= 1; min 0.85), finest rel err {:.2e}", (r.tau_t.0 - tau_t_ex) / tau_t_ex), &mut fail_list);
        let mres = res_list.iter().map(|r| r.mass_res).fold(0.0, f64::max);
        let conv = res_list.iter().all(|r| r.conv);
        let itr_ok = res_list.iter().zip(itr_list.iter()).all(|(r, it)| it.0 < 0.01 * r.err_u.1);
        check(&format!("{name} continuity residual"), mres < 1e-6, format!("max |r_mass| / (rho U h) = {mres:.2e} (tol 1e-6)"), &mut fail_list);
        check(&format!("{name} solver convergence"), conv, format!("picard reached tol {tol_study:e}: {conv}"), &mut fail_list);
        check(&format!("{name} iterative vs discretization error"), itr_ok, format!("L2(vel_tol1e-9 - vel_tol1e-11) < 0.01 L2 error on all meshes: {itr_ok}"), &mut fail_list);
    }}

    println!("\nbench_flow: {} failed criteria", fail_list.len());
    for f in fail_list.iter() { println!("  FAIL {f}"); }
    if !fail_list.is_empty() { std::process::exit(1); }
}
