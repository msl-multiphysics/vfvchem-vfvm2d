#[path = "common.rs"]
mod common;
use common::*;
use std::collections::HashMap;
use std::f64::consts::PI;
use vfvchem_vfvm2d::*;

/// Benchmark: manufactured advection-diffusion-reaction solution for heat (SteadyHeatFlow) and species (SteadyMoleFlow)
/// Run with: `cargo run --release --example bench_heatflow`
///
/// The velocity is prescribed (vecdom function, no flow domain) and is not solenoidal.
/// All properties vary in space.
///
/// Equations solved by the code (non-conservative form)
/// heat:    rho * cp * vel . grad(T) = div( k * grad(T) ) + L * T + S
/// species: rho * vel . grad(c) = div( rho * D * grad(c) ) + L * c + S
/// generic: a * vel . grad(u) = div( g * grad(u) ) + L * u + S with a = rho * cp or rho, g = k or rho * D
///
/// Manufactured solution (s = x / X, t = y / X, X = 0.1 m)
/// phi = sin(c0 s + c1) * cos(c2 t + c3)
/// ext:  (c0, c1, c2, c3) = (1.5, 0.3, 1.2, -0.4); smooth maximum inside the domain (exercises the slope limiter)
/// mono: (c0, c1, c2, c3) = (0.6, 0.2, 0.5, -0.6); monotone in x and y, no interior extremum
/// heat:    T = 300 + 50 * phi; rho = 1 * (1 + 0.3 s); cp = 1000 * (1 + 0.2 t); k = 5e-3 * (1 + 0.5 s + 0.3 t); L = -1.0 * (1 + s t)
/// species: c = 0.5 + 0.3 * phi; rho = 40 * (1 + 0.3 s); D = 5e-6 * (1 + 0.5 s + 0.3 t); L = -0.04 * (1 + s t)
/// S = a * vel . grad(u) - g * lap(u) - grad(g) . grad(u) - L * u
///
/// Velocity (U = 1e-3 m s-1; Peclet number U X a / g = 20 for both equations at the origin)
/// v1   = U * ( (1 + 0.5 t) * (1 + 0.3 s), 0.4 * sin(pi s) * (1 - t) ); div(v1) != 0; v1 . n = 0 on the top boundary
/// vpar = U * (1 + 0.5 s + 0.3 t) * (1, 1) / sqrt(2); parallel to the diagonal x = y
///
/// Cases
/// dir     - square (new_bounds), u exact on all boundaries, v1
/// neu     - square, exact outward diffusive flux on the top boundary (v1 . n = 0 there), u exact elsewhere
/// out     - square, exact outward diffusive flux on the right boundary where v1 . n > 0 (outflow), u exact elsewhere
/// itf_par - obj_diag (two triangles, interface on x = y), continuity interface, u exact on all boundaries, vpar (no flow across the interface)
/// itf_x   - obj_diag, as itf_par but with v1, which crosses the interface

const LEN: f64 = 0.1;
const VEL: f64 = 1e-3;

// coefficients of phi = sin(c0 s + c1) * cos(c2 t + c3)
const SOL_EXT: [f64; 4] = [1.5, 0.3, 1.2, -0.4];   // maximum inside the domain at s = 0.85, t = 0.33
const SOL_MONO: [f64; 4] = [0.6, 0.2, 0.5, -0.6];  // monotone in s and t (no interior extremum)

fn phi(cf: [f64; 4], x: f64, y: f64) -> (f64, (f64, f64), f64) {
    // value, gradient, laplacian
    let (s, t) = (x / LEN, y / LEN);
    let (sa, ca) = (cf[0] * s + cf[1]).sin_cos();
    let (sb, cb) = (cf[2] * t + cf[3]).sin_cos();
    let val = sa * cb;
    let grad = (cf[0] / LEN * ca * cb, -cf[2] / LEN * sa * sb);
    let lap = -(cf[0] * cf[0] + cf[2] * cf[2]) / (LEN * LEN) * val;
    (val, grad, lap)
}

fn vel_func(case: &str, x: f64, y: f64) -> (f64, f64) {
    let (s, t) = (x / LEN, y / LEN);
    if case == "itf_par" {
        let m = VEL * (1.0 + 0.5 * s + 0.3 * t) / 2.0f64.sqrt();
        (m, m)
    } else {
        (VEL * (1.0 + 0.5 * t) * (1.0 + 0.3 * s), VEL * 0.4 * (PI * s).sin() * (1.0 - t))
    }
}

// properties: (u0, du, rho_0, cp_0, g_0, lin_0); species use cp = 1 and g = rho * D
fn param(eqn: &str) -> (f64, f64, f64, f64, f64, f64) {
    if eqn == "heat" { (300.0, 50.0, 1.0, 1000.0, 5e-3, -1.0) } else { (0.5, 0.3, 40.0, 1.0, 5e-6, -0.04) }
}

fn den_func(eqn: &str, x: f64, _y: f64) -> f64 { param(eqn).2 * (1.0 + 0.3 * x / LEN) }
fn hcap_func(eqn: &str, _x: f64, y: f64) -> f64 { if eqn == "heat" { param(eqn).3 * (1.0 + 0.2 * y / LEN) } else { 1.0 } }
fn coef_func(eqn: &str, x: f64, y: f64) -> f64 { param(eqn).4 * (1.0 + 0.5 * x / LEN + 0.3 * y / LEN) }  // k or D
fn lin_func(eqn: &str, x: f64, y: f64) -> f64 { param(eqn).5 * (1.0 + x * y / (LEN * LEN)) }

fn exact(eqn: &str, cf: [f64; 4], x: f64, y: f64) -> f64 {
    let (u0, du, _, _, _, _) = param(eqn);
    u0 + du * phi(cf, x, y).0
}

fn gam_func(eqn: &str, x: f64, y: f64) -> (f64, (f64, f64)) {
    // diffusive coefficient g and its gradient
    let g0 = param(eqn).4;
    let kd = coef_func(eqn, x, y);
    let dkd = (0.5 * g0 / LEN, 0.3 * g0 / LEN);
    if eqn == "heat" { return (kd, dkd); }
    let den = den_func(eqn, x, y);
    let dden = (0.3 * param(eqn).2 / LEN, 0.0);
    (den * kd, (dden.0 * kd + den * dkd.0, dden.1 * kd + den * dkd.1))
}

fn src_func(eqn: &str, cf: [f64; 4], case: &str, x: f64, y: f64) -> f64 {
    let du = param(eqn).1;
    let (val, (px, py), lap) = phi(cf, x, y);
    let (ux, uy, ulap) = (du * px, du * py, du * lap);
    let u = param(eqn).0 + du * val;
    let a = den_func(eqn, x, y) * hcap_func(eqn, x, y);
    let (vx, vy) = vel_func(case, x, y);
    let (g, (gx, gy)) = gam_func(eqn, x, y);
    a * (vx * ux + vy * uy) - g * ulap - (gx * ux + gy * uy) - lin_func(eqn, x, y) * u
}

fn flux_out(eqn: &str, cf: [f64; 4], x: f64, y: f64, n: (f64, f64)) -> f64 {
    // exact outward diffusive flux -g grad(u) . n
    let du = param(eqn).1;
    let (_, (px, py), _) = phi(cf, x, y);
    -gam_func(eqn, x, y).0 * du * (px * n.0 + py * n.1)
}

fn int_side(eqn: &str, cf: [f64; 4], p0: (f64, f64), p1: (f64, f64), n: (f64, f64)) -> f64 {
    // exact inward diffusive flow g grad(u) . n_out over a straight side (midpoint rule)
    let num = 4000;
    let len = ((p1.0 - p0.0).powi(2) + (p1.1 - p0.1).powi(2)).sqrt();
    let mut sum = 0.0;
    for i in 0..num {
        let r = (i as f64 + 0.5) / num as f64;
        sum -= flux_out(eqn, cf, p0.0 + r * (p1.0 - p0.0), p0.1 + r * (p1.1 - p0.1), n);
    }
    sum * len / num as f64
}

struct Result {
    h: f64,
    num_vert: usize,
    err: (f64, f64, f64),   // L1, L2, Linf
    err_far: f64,           // Linf over vertices farther than 6 h from boundaries and the interface
    loc_max: (f64, f64),    // location of the largest error
    flow_side: Vec<f64>,    // inward diffusive flow per boundary side (bnd_l, bnd_r, bnd_b, bnd_t) from the residual
    flow_itf: (f64, f64),   // inward diffusive flow into dom_br and dom_tl through the interface
    itf_mis: f64,           // max |r_a + r_b| / sum |terms| over interface nodes (r = residual of each domain's operators)
    over: f64,              // overshoot beyond the exact range
    conv: bool,
    unk_vec: Vec<f64>,      // solution of all domains for the iterative-error comparison
}

fn run(eqn: &'static str, cf: [f64; 4], case: &'static str, h: f64, tol: f64) -> Result {
    // variables and geometry
    let two = case.starts_with("itf");
    let mut vars = if two {
        Variables::new_obj("examples/obj/obj_diag.obj", h, h, "target/bench_output/heatflow").unwrap()
    } else {
        Variables::new_bounds(0.0, 0.0, LEN, LEN, h, h, "target/bench_output/heatflow").unwrap()
    };
    let mut dom_list = Vec::new();
    let mut bnd_list = Vec::new();  // (domain index, bnd_id) ordered l, r, b, t
    let mut itf = 0;
    if two {
        let dom_br = vars.add_dom("dom_br").unwrap();
        let dom_tl = vars.add_dom("dom_tl").unwrap();
        bnd_list.push((1, vars.add_bnd(dom_tl, "bnd_l").unwrap()));
        bnd_list.push((0, vars.add_bnd(dom_br, "bnd_r").unwrap()));
        bnd_list.push((0, vars.add_bnd(dom_br, "bnd_b").unwrap()));
        bnd_list.push((1, vars.add_bnd(dom_tl, "bnd_t").unwrap()));
        itf = vars.add_itf(dom_br, dom_tl, "itf_d").unwrap();
        dom_list.push(dom_br);
        dom_list.push(dom_tl);
    } else {
        let dom = vars.add_dom("dom").unwrap();
        for name in ["bnd_l", "bnd_r", "bnd_b", "bnd_t"] { bnd_list.push((0, vars.add_bnd(dom, name).unwrap())); }
        dom_list.push(dom);
    }

    // unknowns and properties per domain
    let (u0, _, _, _, _, _) = param(eqn);
    let mut unk_list = Vec::new();
    let mut prop_list = Vec::new();  // (den, hcap, vel, coef, lin, src)
    for &dom in dom_list.iter() {
        let unk = vars.add_scldom_unkn(dom, u0, "unk").unwrap();
        let den = vars.add_scldom_func(dom, Box::new(move |_t, (x, y), _u| den_func(eqn, x, y)), vec![], "den").unwrap();
        let hcap = vars.add_scldom_func(dom, Box::new(move |_t, (x, y), _u| hcap_func(eqn, x, y)), vec![], "hcap").unwrap();
        let vel = vars.add_vecdom_func(dom, Box::new(move |_t, (x, y), _u| vel_func(case, x, y)), vec![], "vel").unwrap();
        let coef = vars.add_scldom_func(dom, Box::new(move |_t, (x, y), _u| coef_func(eqn, x, y)), vec![], "coef").unwrap();
        let lin = vars.add_scldom_func(dom, Box::new(move |_t, (x, y), _u| lin_func(eqn, x, y)), vec![], "lin").unwrap();
        let src = vars.add_scldom_func(dom, Box::new(move |_t, (x, y), _u| src_func(eqn, cf, case, x, y)), vec![], "src").unwrap();
        unk_list.push(unk);
        prop_list.push((den, hcap, vel, coef, lin, src));
    }

    // boundary conditions; flux side: neu -> top, out -> right
    let side_flx = if case == "neu" { 3 } else if case == "out" { 1 } else { 9 };
    let norm_list = [(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0)];
    let mut bc_list = Vec::new();  // (domain index, bnd_id, value id, is flux)
    for (side, &(di, bnd)) in bnd_list.iter().enumerate() {
        if side == side_flx {
            let n = norm_list[side];
            let flx = vars.add_sclbnd_func(bnd, Box::new(move |_t, (x, y), _u| flux_out(eqn, cf, x, y, n)), vec![]).unwrap();
            bc_list.push((di, bnd, flx, true));
        } else {
            let set = vars.add_sclbnd_func(bnd, Box::new(move |_t, (x, y), _u| exact(eqn, cf, x, y)), vec![]).unwrap();
            bc_list.push((di, bnd, set, false));
        }
    }

    // physics and solver
    let mut lin_sol = RslabLu::new(4).unwrap();
    let conv = if eqn == "heat" {
        let mut phys = SteadyHeatFlow::new().unwrap();
        for (i, &dom) in dom_list.iter().enumerate() {
            let (den, hcap, vel, coef, lin, src) = prop_list[i];
            phys.add_heat_dom(dom, unk_list[i], den, hcap, vel, coef, lin, src).unwrap();
        }
        for &(_, bnd, val, is_flx) in bc_list.iter() {
            if is_flx { phys.add_hflx_bnd(bnd, val).unwrap(); } else { phys.add_temp_bnd(bnd, val).unwrap(); }
        }
        if two { phys.add_hcnt_itf(itf).unwrap(); }
        phys.solve(&mut vars, &mut lin_sol, 200, tol, 1.0).is_ok()
    } else {
        let mut phys = SteadyMoleFlow::new(1).unwrap();
        for (i, &dom) in dom_list.iter().enumerate() {
            let (den, _, vel, coef, lin, src) = prop_list[i];
            phys.add_mole_dom(0, dom, unk_list[i], den, vel, HashMap::from([(0, coef)]), lin, src).unwrap();
        }
        for &(_, bnd, val, is_flx) in bc_list.iter() {
            if is_flx { phys.add_nflx_bnd(0, bnd, val).unwrap(); } else { phys.add_conc_bnd(0, bnd, val).unwrap(); }
        }
        if two { phys.add_ncnt_itf(0, itf).unwrap(); }
        phys.solve(&mut vars, &mut lin_sol, 200, tol, 1.0).is_ok()
    };

    // errors and boundedness
    let mut err = Vec::new();
    let mut area = Vec::new();
    let mut err_far: f64 = 0.0;
    let mut loc_max = (0.0, 0.0);
    let mut emax: f64 = 0.0;
    let mut unk_vec = Vec::new();
    let (mut umin, mut umax) = (f64::INFINITY, f64::NEG_INFINITY);
    for (i, &dom) in dom_list.iter().enumerate() {
        for (vid, &(x, y)) in vert_pos(&vars, dom).iter().enumerate() {
            let u = vars.scldom[unk_list[i]].value[vid];
            let e = u - exact(eqn, cf, x, y);
            let dist = x.min(y).min(LEN - x).min(LEN - y).min(if two { (x - y).abs() / 2.0f64.sqrt() } else { f64::INFINITY });
            if dist > 6.0 * h { err_far = err_far.max(e.abs()); }
            if e.abs() > emax { emax = e.abs(); loc_max = (x, y); }
            umin = umin.min(u);
            umax = umax.max(u);
            err.push(e);
            area.push(vars.tri[dom].vert_area[vid]);
            unk_vec.push(u);
        }
    }
    let num = 400;
    let (mut emin, mut emax_ex) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..=num { for j in 0..=num {
        let u = exact(eqn, cf, LEN * i as f64 / num as f64, LEN * j as f64 / num as f64);
        emin = emin.min(u);
        emax_ex = emax_ex.max(u);
    }}
    let over = (umax - emax_ex).max(emin - umin).max(0.0);

    // residual of the domain operators (before boundary and interface rows) per domain
    let mut res_dom = Vec::new();
    for (i, &dom) in dom_list.iter().enumerate() {
        let (den, hcap, vel, coef, lin, src) = prop_list[i];
        let unk = unk_list[i];
        let (res, scl) = if eqn == "heat" {
            let o1 = HeatDomAdvection::new(&vars, dom, unk, vel, den, hcap).unwrap();
            let o2 = HeatDomDiffusion::new(&vars, dom, unk, coef).unwrap();
            let o3 = SclDomLinear::new(&vars, dom, unk, lin).unwrap();
            let o4 = SclDomSource::new(&vars, dom, unk, src).unwrap();
            get_res(&vars, &[&o1, &o2, &o3, &o4])
        } else {
            let o1 = MoleDomAdvection::new(&vars, dom, unk, vel, den).unwrap();
            let o2 = MoleDomDiffusion::new(&vars, dom, unk, den, coef, unk).unwrap();
            let o3 = SclDomLinear::new(&vars, dom, unk, lin).unwrap();
            let o4 = SclDomSource::new(&vars, dom, unk, src).unwrap();
            get_res(&vars, &[&o1, &o2, &o3, &o4])
        };
        let start = scl_start(&vars, unk);
        res_dom.push((res[start..start + vars.dom[dom].num_vert()].to_vec(), scl[start..start + vars.dom[dom].num_vert()].to_vec()));
    }

    // inward flow per side; a node shared by two sides counts for the first
    // includes the advective term omitted at the boundary half-edges (see NOTES)
    let mut flow_side = vec![0.0; 4];
    let mut seen: Vec<Vec<bool>> = dom_list.iter().map(|&d| vec![false; vars.dom[d].num_vert()]).collect();
    for (side, &(di, bnd)) in bnd_list.iter().enumerate() {
        for vid in bnd_vid(&vars, bnd) {
            if seen[di][vid] { continue; }
            seen[di][vid] = true;
            flow_side[side] -= res_dom[di].0[vid];
        }
    }

    // interface flow; excludes interface nodes on boundaries
    let mut flow_itf = (0.0, 0.0);
    let mut itf_mis: f64 = 0.0;
    if two {
        let vid_a = itf_vid(&vars, itf, dom_list[0]);
        let vid_b = itf_vid(&vars, itf, dom_list[1]);
        for k in 0..vid_a.len() {
            if seen[0][vid_a[k]] || seen[1][vid_b[k]] { continue; }
            let (ra, rb) = (res_dom[0].0[vid_a[k]], res_dom[1].0[vid_b[k]]);
            flow_itf.0 -= ra;
            flow_itf.1 -= rb;
            itf_mis = itf_mis.max((ra + rb).abs() / (res_dom[0].1[vid_a[k]] + res_dom[1].1[vid_b[k]]));
        }
    }

    Result { h, num_vert: unk_vec.len(), err: err_norm(&err, &area), err_far, loc_max, flow_side, flow_itf, itf_mis, over, conv, unk_vec }
}

fn main() {
    let h_list = [0.01, 0.005, 0.0025, 0.00125, 0.000625, 0.0003125];
    let tol_study = 1e-10;
    let tol_ref = 1e-12;
    let mut fail_list = Vec::new();

    for (sol, cf) in [("ext", SOL_EXT), ("mono", SOL_MONO)] {
    for eqn in ["heat", "mole"] {
        // exact inward diffusive flows for the full square sides and the diagonal (into dom_br, n_out = (-1, 1) / sqrt 2)
        let flow_ex = [
            int_side(eqn, cf, (0.0, 0.0), (0.0, LEN), (-1.0, 0.0)), int_side(eqn, cf, (LEN, 0.0), (LEN, LEN), (1.0, 0.0)),
            int_side(eqn, cf, (0.0, 0.0), (LEN, 0.0), (0.0, -1.0)), int_side(eqn, cf, (0.0, LEN), (LEN, LEN), (0.0, 1.0)),
        ];
        let s2 = 0.5f64.sqrt();
        let fitf_ex = int_side(eqn, cf, (0.0, 0.0), (LEN, LEN), (-s2, s2));
        let flow_ref = flow_ex.iter().map(|f| f.abs()).fold(0.0, f64::max);
        let du = param(eqn).1;
        println!("\n##### bench_heatflow solution {sol} {cf:?}, equation {eqn}: exact inward diffusive flows (per unit depth) l {:.4e}, r {:.4e}, b {:.4e}, t {:.4e}, interface into dom_br {:.4e}",
            flow_ex[0], flow_ex[1], flow_ex[2], flow_ex[3], fitf_ex);

        for case in ["dir", "neu", "out", "itf_par", "itf_x"] {
            let mut res_list = Vec::new();
            let mut itr_list = Vec::new();
            for &h in h_list.iter() {
                let r = run(eqn, cf, case, h, tol_study);
                let r_ref = run(eqn, cf, case, h, tol_ref);
                let dif: Vec<f64> = r.unk_vec.iter().zip(r_ref.unk_vec.iter()).map(|(a, b)| a - b).collect();
                itr_list.push((err_norm(&dif, &vec![1.0; dif.len()]).1, r_ref.conv, r_ref.itf_mis));
                res_list.push(r);
            }

            // table; errors relative to the amplitude du; flows relative to the largest side flow
            println!("\n=== bench_heatflow {sol} {eqn} case {case} (errors / {du}; flow errors / {flow_ref:.3e}) ===");
            println!("{:>9} {:>6} {:>9} {:>9} {:>9} {:>5} {:>5} {:>5} {:>9} {:>17} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>5}",
                "h", "N", "L1", "L2", "Linf", "o_L1", "o_L2", "o_inf", "Linf_far", "loc_max", "Qerr_l", "Qerr_r", "Qerr_b", "Qerr_t", "Qerr_itf", "itf_mis", "overshoot", "it_L2", "Q_bal", "conv");
            // itf_mis is for the study solution (tol 1e-10); the continuity criterion uses the tol 1e-12 solution
            for (i, r) in res_list.iter().enumerate() {
                let (o1, o2, oi) = if i == 0 { (f64::NAN, f64::NAN, f64::NAN) } else {
                    let p = &res_list[i - 1];
                    (get_order(p.err.0, r.err.0, p.h, r.h), get_order(p.err.1, r.err.1, p.h, r.h), get_order(p.err.2, r.err.2, p.h, r.h))
                };
                let qitf = if case.starts_with("itf") { (r.flow_itf.0 - fitf_ex) / flow_ref } else { f64::NAN };
                let qbal = (r.flow_side.iter().sum::<f64>() - flow_ex.iter().sum::<f64>()) / flow_ref;
                println!("{:>9.3e} {:>6} {:>9.2e} {:>9.2e} {:>9.2e} {:>5.2} {:>5.2} {:>5.2} {:>9.2e} {:>8.5},{:>8.5} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>5}",
                    r.h, r.num_vert, r.err.0 / du, r.err.1 / du, r.err.2 / du, o1, o2, oi, r.err_far / du, r.loc_max.0, r.loc_max.1,
                    (r.flow_side[0] - flow_ex[0]) / flow_ref, (r.flow_side[1] - flow_ex[1]) / flow_ref, (r.flow_side[2] - flow_ex[2]) / flow_ref, (r.flow_side[3] - flow_ex[3]) / flow_ref,
                    qitf, r.itf_mis, r.over / du, itr_list[i].0 / du, qbal, r.conv && itr_list[i].1);
            }

            // pass/fail criteria
            // expected: 2nd order in L1/L2 and >= 1st order in Linf; total boundary diffusive flow (balance) >= 1st order
            // interface continuity of the diffusive flow at round-off; overshoot within the Linf error; iterative error < 1 % of the L2 error
            let name = format!("{sol}_{eqn}_{case}");
            // orders are fitted over the meshes with h <= 5e-3 (the 0.01 mesh has 10 cells per side and is pre-asymptotic)
            let sel: Vec<&Result> = res_list.iter().filter(|r| r.h <= 0.005).collect();
            let h_vec: Vec<f64> = sel.iter().map(|r| r.h).collect();
            let fit = |f: &dyn Fn(&Result) -> f64| fit_order(&sel.iter().map(|r| f(r).abs().max(1e-300)).collect::<Vec<f64>>(), &h_vec);
            let o1 = fit(&|r| r.err.0);
            let o2 = fit(&|r| r.err.1);
            let oi = fit(&|r| r.err.2);
            let sel_far: Vec<&Result> = res_list.iter().filter(|r| r.err_far > 0.0 && r.h <= 0.005).collect();
            let of = fit_order(&sel_far.iter().map(|r| r.err_far).collect::<Vec<f64>>(), &sel_far.iter().map(|r| r.h).collect::<Vec<f64>>());
            let ob = fit(&|r| r.flow_side.iter().sum::<f64>() - flow_ex.iter().sum::<f64>());
            println!("fitted orders (h <= 5e-3): L1 {o1:.2}, L2 {o2:.2}, Linf {oi:.2}, Linf_far {of:.2}, total boundary flow {ob:.2}");
            check(&format!("{name} L1 order"), o1 >= 1.7, format!("{o1:.2} (expected 2; min 1.7)"), &mut fail_list);
            check(&format!("{name} L2 order"), o2 >= 1.7, format!("{o2:.2} (expected 2; min 1.7)"), &mut fail_list);
            check(&format!("{name} Linf order"), oi >= 0.85, format!("{oi:.2} (expected >= 1; min 0.85)"), &mut fail_list);
            check(&format!("{name} total boundary flow order"), ob >= 0.85, format!("{ob:.2} (expected >= 1; min 0.85)"), &mut fail_list);
            if case.starts_with("itf") {
                // the interface row is the sum of both domain rows, so r_a + r_b is the unconverged residual of that row
                let mis = itr_list.iter().map(|it| it.2).fold(0.0, f64::max);
                check(&format!("{name} interface flux continuity"), mis < 1e-8, format!("max |r_a + r_b| / sum |terms| = {mis:.2e} at tol {tol_ref:e} (tol 1e-8)"), &mut fail_list);
                let oq = fit(&|r| r.flow_itf.0 - fitf_ex);
                check(&format!("{name} interface flow order"), oq >= 0.85, format!("{oq:.2} (expected >= 1; min 0.85)"), &mut fail_list);
            }
            let bnd_ok = res_list.iter().all(|r| r.over <= r.err.2);
            check(&format!("{name} boundedness"), bnd_ok, format!("overshoot beyond the exact range <= Linf error on all meshes: {bnd_ok}"), &mut fail_list);
            let conv = res_list.iter().all(|r| r.conv);
            check(&format!("{name} solver convergence"), conv, format!("picard reached tol {tol_study:e}: {conv}"), &mut fail_list);
            let itr_ok = res_list.iter().zip(itr_list.iter()).all(|(r, it)| it.0 < 0.01 * r.err.1 && it.1);
            check(&format!("{name} iterative vs discretization error"), itr_ok, format!("L2(u_tol1e-10 - u_tol1e-12) < 0.01 L2 error on all meshes: {itr_ok}"), &mut fail_list);
        }
    }

    }

    println!("\nbench_heatflow: {} failed criteria", fail_list.len());
    for f in fail_list.iter() { println!("  FAIL {f}"); }
    if !fail_list.is_empty() { std::process::exit(1); }
}
