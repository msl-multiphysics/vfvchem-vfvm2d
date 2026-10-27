use vfvchem_vfvm2d::*;

// local thermal non-equilibrium heat transfer in porous media (SteadyHeatPore, TransientHeatPore)
// fluid: epsl * (rho cp)_f * dT_f/dt + div( (rho cp)_f * vel * T_f ) = div( (epsl * k_f * I + K_d) * grad(T_f) ) + h * a_v * (T_s - T_f) + epsl * (s_f * T_f + q_f)
// solid: (1 - epsl) * (rho c)_s * dT_s/dt = div( (1 - epsl) * k_s * grad(T_s) ) - h * a_v * (T_s - T_f) + (1 - epsl) * (s_s * T_s + q_s) + a_v * (s_a * T_s + q_a)
// run with: cargo test --release --test heatpore -- --nocapture --test-threads 1

// positions of all vertices in a domain; vid -> (x, y)
fn positions(vars: &Variables, dom: usize) -> Vec<(f64, f64)> {
    let mut pos: Vec<(f64, f64)> = vars.dom[dom].pts_gpt.iter().map(|&g| vars.pts.pts_pos[g]).collect();
    pos.extend(vars.dom[dom].nds_gnd.iter().map(|&g| vars.nds.nds_pos[g]));
    pos
}

// channel [0, l] x [0, h] with four boundaries
struct Mesh { vars: Variables, dom: usize, bnd: [usize; 4] }

fn mesh(l: f64, h: f64, dx: f64) -> Result<Mesh, VFVChemError> {
    let mut vars = Variables::new_bounds(0.0, 0.0, l, h, dx, dx, dx, 100, 42)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;
    Ok(Mesh { vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t] })
}

// dense 4 x 4 solve with partial pivoting
fn solve4(mut a: [[f64; 4]; 4], mut b: [f64; 4]) -> [f64; 4] {
    for k in 0..4 {
        let p = (k..4).max_by(|&i, &j| a[i][k].abs().partial_cmp(&a[j][k].abs()).unwrap()).unwrap();
        a.swap(k, p); b.swap(k, p);
        for i in k + 1..4 {
            let f = a[i][k] / a[k][k];
            for j in k..4 { a[i][j] -= f * a[k][j]; }
            b[i] -= f * b[k];
        }
    }
    let mut x = [0.0; 4];
    for k in (0..4).rev() {
        x[k] = (b[k] - (k + 1..4).map(|j| a[k][j] * x[j]).sum::<f64>()) / a[k][k];
    }
    x
}

// ---------------------------------------------------------------- epsl = 1 regression

// Poiseuille channel (Re = 50) with conduction, a production coefficient, and a non-uniform source
// heatflow solver and heatpore solver with epsl = 1, K = inf, F = 0, alpha_L = alpha_T = 0
// h and a_v are non-zero and the solid has its own data and BCs; the solid mask must remove all of it
fn poiseuille(pore: bool, transient: bool) -> Result<Vec<f64>, VFVChemError> {
    let Mesh { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t] } = mesh(0.02, 0.005, 5e-4)?;
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 0.001, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    let epsl = vars.add_scldom_cnst(dom, 1.0, "")?;
    let perm = vars.add_scldom_cnst(dom, f64::INFINITY, "")?;  // free fluid
    let forch = vars.add_scldom_cnst(dom, 0.0, "")?;
    let temp = vars.add_scldom_unkn(dom, 300.0, "")?;
    let vlcp = vars.add_scldom_cnst(dom, 4.18e6, "")?;
    // spatially varying conductivity exercises the skew correction
    let cond = vars.add_scldom_func(dom, Box::new(|_t, (x, _y), _u| 0.6 * (1.0 + 20.0 * x)), vec![], "")?;
    let hprd = vars.add_scldom_cnst(dom, -2e3, "")?;
    let hsrc = vars.add_scldom_func(dom, Box::new(|_t, (x, y), _u| 2e6 * (x / 0.02) * (y / 0.005) + 2e3 * 300.0), vec![], "")?;
    // the heatflow solver must not see an extra unknown
    let stemp = if pore { vars.add_scldom_unkn(dom, 350.0, "")? } else { vars.add_scldom_cnst(dom, 350.0, "")? };
    let svlcp = vars.add_scldom_cnst(dom, 2e6, "")?;
    let scond = vars.add_scldom_cnst(dom, 20.0, "")?;
    let sprd = vars.add_scldom_cnst(dom, -7.0, "")?;
    let ssrc = vars.add_scldom_cnst(dom, 3e5, "")?;
    let htc = vars.add_scldom_cnst(dom, 500.0, "")?;
    let area = vars.add_scldom_cnst(dom, 600.0, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _u| (4.0 * 0.01 * y * (0.005 - y) / 0.005f64.powi(2), 0.0));
    let vel_l = vars.add_vecbnd_func(bnd_l, vel_func, vec![], "")?;
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;
    let temp_l = vars.add_sclbnd_func(bnd_l, Box::new(|_t, (_x, y), _u| 300.0 + 20.0 * y / 0.005), vec![], "")?;
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, -5e3, "")?;
    let stemp_l = vars.add_sclbnd_cnst(bnd_l, 400.0, "")?;
    let shflx_r = vars.add_sclbnd_cnst(bnd_r, 1e3, "")?;
    let shflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let shflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;
    let mut lin = FaerLu::new(1)?;
    macro_rules! setup {
        ($phys:ident) => {
            $phys.add_vel_bnd(bnd_l, vel_l)?; $phys.add_pres_bnd(bnd_r, pres_r)?; $phys.add_vel_bnd(bnd_b, vel_b)?; $phys.add_vel_bnd(bnd_t, vel_t)?;
        };
    }
    macro_rules! setup_pore {
        ($phys:ident) => {
            $phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, forch)?;
            $phys.add_heatl_dom(dom, temp, vlcp, cond, hprd, hsrc)?;
            $phys.add_heats_dom(dom, stemp, svlcp, scond, sprd, ssrc, zero, ssrc)?;
            $phys.add_matrix_dom(dom, htc, area, zero, zero)?;
            $phys.add_templ_bnd(bnd_l, temp_l)?; $phys.add_houtl_bnd(bnd_r)?; $phys.add_hflxl_bnd(bnd_b, hflx_b)?; $phys.add_hflxl_bnd(bnd_t, hflx_t)?;
            $phys.add_temps_bnd(bnd_l, stemp_l)?; $phys.add_hflxs_bnd(bnd_r, shflx_r)?; $phys.add_hflxs_bnd(bnd_b, shflx_b)?; $phys.add_hflxs_bnd(bnd_t, shflx_t)?;
        };
    }
    match (pore, transient) {
        (false, false) => {
            let mut phys = SteadyHeatFlow::new()?;
            phys.add_heat_dom(dom, temp, vlcp, cond, hprd, hsrc)?;
            phys.add_temp_bnd(bnd_l, temp_l)?; phys.add_hout_bnd(bnd_r)?; phys.add_hflx_bnd(bnd_b, hflx_b)?; phys.add_hflx_bnd(bnd_t, hflx_t)?;
            phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
            setup!(phys);
            SteadyBase::solve(&mut phys, &mut vars, &mut lin, 50, 1e-10, 1.0)?;
        }
        (true, false) => {
            let mut phys = SteadyHeatPore::new()?;
            setup_pore!(phys);
            setup!(phys);
            SteadyBase::solve(&mut phys, &mut vars, &mut lin, 50, 1e-10, 1.0)?;
            let dts = vars.scldom[temp].value.iter().zip(&vars.scldom[stemp].value).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            println!("RESULT regression steady: max|T_s - T_f| = {dts:.3e}");
            assert!(dts == 0.0);
        }
        (false, true) => {
            let mut phys = TransientHeatFlow::new()?;
            phys.add_heat_dom(dom, temp, vlcp, cond, hprd, hsrc)?;
            phys.add_temp_bnd(bnd_l, temp_l)?; phys.add_hout_bnd(bnd_r)?; phys.add_hflx_bnd(bnd_b, hflx_b)?; phys.add_hflx_bnd(bnd_t, hflx_t)?;
            phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
            setup!(phys);
            TransientBase::solve(&mut phys, &mut vars, &mut lin, 0.05, 4, 100, 50, 1e-10, 1.0)?;
        }
        (true, true) => {
            let mut phys = TransientHeatPore::new()?;
            setup_pore!(phys);
            setup!(phys);
            TransientBase::solve(&mut phys, &mut vars, &mut lin, 0.05, 4, 100, 50, 1e-10, 1.0)?;
        }
    }
    Ok(vars.scldom[temp].value.clone())
}

fn regression(transient: bool) {
    let temp_f = poiseuille(false, transient).unwrap();
    let temp_p = poiseuille(true, transient).unwrap();
    let tmax = temp_f.iter().map(|a| a.abs()).fold(0.0, f64::max);
    let dt = temp_f.iter().zip(&temp_p).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    let trng = temp_f.iter().cloned().fold(f64::MIN, f64::max) - temp_f.iter().cloned().fold(f64::MAX, f64::min);
    println!("RESULT regression transient={transient}: T_max = {tmax:.4e}; T range = {trng:.4e}; max|dT|/T_max = {:.3e}", dt / tmax);
    assert!(dt / tmax < 1e-12);
}

#[test]
fn regression_steady() { regression(false); }

#[test]
fn regression_transient() { regression(true); }

// ---------------------------------------------------------------- plug flow through a packed column

// plug flow (superficial u_in) on [0, l] x [0, h]; vel BC on left, top, bottom; p = 0 on the right
struct Column { vars: Variables, dom: usize, temp: usize, stemp: usize }

// fluid: T_f(0) = 0, outflow on the right; solid: adiabatic everywhere; top and bottom adiabatic for both phases
// (rho cp)_f = 1 and u = 1 so the advective coefficient is 1
fn column(dx: f64, epsl_c: f64, kf: f64, ks: f64, h: f64, av: f64, ssrc_c: f64, sprd_c: f64) -> Result<Column, VFVChemError> {
    let Mesh { mut vars, dom, bnd } = mesh(1.0, 0.2, dx)?;
    let [bnd_l, bnd_r, bnd_b, bnd_t] = bnd;
    let u_in = 1.0;
    let vel = vars.add_vecdom_unkn(dom, (u_in, 0.0), "")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom, 1.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 1.0, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    let epsl = vars.add_scldom_cnst(dom, epsl_c, "")?;
    let perm = vars.add_scldom_cnst(dom, 1e-2, "")?;
    let forch = vars.add_scldom_cnst(dom, 0.5, "")?;
    let temp = vars.add_scldom_unkn(dom, 0.0, "")?;
    let stemp = vars.add_scldom_unkn(dom, 0.0, "")?;
    let one = vars.add_scldom_cnst(dom, 1.0, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let cond = vars.add_scldom_cnst(dom, kf, "")?;
    let scond = vars.add_scldom_cnst(dom, ks, "")?;
    let ssrc = vars.add_scldom_cnst(dom, ssrc_c, "")?;
    let sprd = vars.add_scldom_cnst(dom, sprd_c, "")?;
    let htc = vars.add_scldom_cnst(dom, h, "")?;
    let area = vars.add_scldom_cnst(dom, av, "")?;
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (u_in, 0.0), "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (u_in, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (u_in, 0.0), "")?;
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;
    let flx0: Vec<usize> = bnd.iter().map(|&b| vars.add_sclbnd_cnst(b, 0.0, "")).collect::<Result<_, _>>()?;
    let mut phys = SteadyHeatPore::new()?;
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, forch)?;
    phys.add_heatl_dom(dom, temp, one, cond, zero, zero)?;
    phys.add_heats_dom(dom, stemp, one, scond, sprd, ssrc, zero, zero)?;
    phys.add_matrix_dom(dom, htc, area, zero, zero)?;
    phys.add_templ_bnd(bnd_l, temp_l)?; phys.add_houtl_bnd(bnd_r)?; phys.add_hflxl_bnd(bnd_b, flx0[2])?; phys.add_hflxl_bnd(bnd_t, flx0[3])?;
    for i in 0..4 { phys.add_hflxs_bnd(bnd[i], flx0[i])?; }
    phys.add_vel_bnd(bnd_l, vel_l)?; phys.add_pres_bnd(bnd_r, pres_r)?; phys.add_vel_bnd(bnd_b, vel_b)?; phys.add_vel_bnd(bnd_t, vel_t)?;
    let mut lin = FaerLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-10, 1.0)?;
    Ok(Column { vars, dom, temp, stemp })
}

// rms and max error of a vertex field against an exact function of x
fn errors(vars: &Variables, dom: usize, val: usize, exact: &dyn Fn(f64) -> f64) -> (f64, f64) {
    let pos = positions(vars, dom);
    let v = &vars.scldom[val].value;
    let rms = (pos.iter().zip(v).map(|(p, v)| (v - exact(p.0)).powi(2)).sum::<f64>() / pos.len() as f64).sqrt();
    let emax = pos.iter().zip(v).map(|(p, v)| (v - exact(p.0)).abs()).fold(0.0, f64::max);
    (rms, emax)
}

// ---------------------------------------------------------------- 1D LTNE column with a heated solid

// steady 1D with B = (rho cp)_f u, A = epsl k_f, C = (1 - epsl) k_s, H = h a_v, Q = (1 - epsl) q_s
// fluid: B T_f' = A T_f'' + H (T_s - T_f); solid: 0 = C T_s'' - H (T_s - T_f) + Q
// T_f(0) = 0, T_f'(1) = 0 (outflow), T_s'(0) = T_s'(1) = 0 (adiabatic)
// particular: T_f = Q x / B, T_s = Q x / B + Q / H; homogeneous: e^{lam x} with lam = 0 or a root of
// A C lam^3 - B C lam^2 - H (A + C) lam + B H = 0; T_s / T_f = H / (H - C lam^2)
fn ltne_exact(a: f64, b: f64, c: f64, hh: f64, q: f64) -> (impl Fn(f64) -> f64, impl Fn(f64) -> f64) {
    // real roots of the cubic by bisection on sign changes
    let f = move |l: f64| a * c * l.powi(3) - b * c * l * l - hh * (a + c) * l + b * hh;
    let mut roots = vec![];
    let (lo, hi, n) = (-1e3, 1e3, 2_000_000);
    let mut x0 = lo;
    for i in 1..=n {
        let x1 = lo + (hi - lo) * i as f64 / n as f64;
        if f(x0) * f(x1) < 0.0 {
            let (mut p, mut r) = (x0, x1);
            for _ in 0..200 { let m = 0.5 * (p + r); if f(p) * f(m) <= 0.0 { r = m } else { p = m } }
            roots.push(0.5 * (p + r));
        }
        x0 = x1;
    }
    assert!(roots.len() == 3);
    let lam = [0.0, roots[0], roots[1], roots[2]];
    let rat: Vec<f64> = lam.iter().map(|&l| hh / (hh - c * l * l)).collect();
    // scaled modes e^{lam (x - x_ref)} with x_ref = 1 for positive lam to avoid overflow
    let xr: Vec<f64> = lam.iter().map(|&l| if l > 0.0 { 1.0 } else { 0.0 }).collect();
    let e = |k: usize, x: f64| (lam[k] * (x - xr[k])).exp();
    let mut m = [[0.0; 4]; 4];
    for k in 0..4 {
        m[0][k] = e(k, 0.0);                          // T_f(0)
        m[1][k] = lam[k] * e(k, 1.0);                 // T_f'(1)
        m[2][k] = rat[k] * lam[k] * e(k, 0.0);        // T_s'(0)
        m[3][k] = rat[k] * lam[k] * e(k, 1.0);        // T_s'(1)
    }
    let g = q / b;
    let coef = solve4(m, [0.0, -g, -g, -g]);
    let (lam2, rat2, xr2, coef2) = (lam, rat.clone(), xr.clone(), coef);
    let tf = move |x: f64| g * x + (0..4).map(|k| coef[k] * (lam[k] * (x - xr[k])).exp()).sum::<f64>();
    let ts = move |x: f64| g * x + q / hh + (0..4).map(|k| coef2[k] * rat2[k] * (lam2[k] * (x - xr2[k])).exp()).sum::<f64>();
    (tf, ts)
}

#[test]
fn ltne_column() {
    let (epsl, kf, ks, h, av, qs) = (0.4, 0.05, 0.5, 2.0, 2.5, 10.0);
    let (a, b, c, hh, q) = (epsl * kf, 1.0, (1.0 - epsl) * ks, h * av, (1.0 - epsl) * qs);
    let (tf_ex, ts_ex) = ltne_exact(a, b, c, hh, q);
    println!("RESULT ltne_column: T_f(1) = {:.6}; T_s(0) = {:.6}; T_s(1) = {:.6}", tf_ex(1.0), ts_ex(0.0), ts_ex(1.0));
    let dxs = [0.05, 0.025, 0.0125];
    let mut errs = vec![];
    for &dx in &dxs {
        let col = column(dx, epsl, kf, ks, h, av, qs, 0.0).unwrap();
        let ef = errors(&col.vars, col.dom, col.temp, &tf_ex);
        let es = errors(&col.vars, col.dom, col.stemp, &ts_ex);
        errs.push((ef, es));
    }
    for i in 0..dxs.len() {
        let pf = if i > 0 { (errs[i - 1].0.0 / errs[i].0.0).ln() / 2f64.ln() } else { f64::NAN };
        let ps = if i > 0 { (errs[i - 1].1.0 / errs[i].1.0).ln() / 2f64.ln() } else { f64::NAN };
        println!("RESULT ltne_column: dx = {:.4}; fluid rms = {:.3e} max = {:.3e} order = {pf:.2}; solid rms = {:.3e} max = {:.3e} order = {ps:.2}", dxs[i], errs[i].0.0, errs[i].0.1, errs[i].1.0, errs[i].1.1);
    }
    let tscale = ts_ex(1.0);
    assert!(errs[2].0.0 < 0.5 * errs[0].0.0 && errs[2].1.0 < 0.5 * errs[0].1.0);
    assert!(errs[2].0.1 < 1e-2 * tscale && errs[2].1.1 < 1e-2 * tscale);
}

// ---------------------------------------------------------------- local thermal equilibrium limit

#[test]
fn lte_limit() {
    // large h a_v forces T_s = T_f; LTE: B T' = (A + C) T'' + (1 - epsl) s_s (T - T_ref) form with s_s T + q_s
    // (1 - epsl) (s_s T + q_s) with s_s = -k_r, q_s = k_r: B T' = D T'' + Kr (1 - T); T(0) = 0, T'(1) = 0
    let (epsl, kf, ks, av, kr): (f64, f64, f64, f64, f64) = (0.4, 0.05, 0.5, 2.5, 2.0);
    let (a, b, c) = (epsl * kf, 1.0, (1.0 - epsl) * ks);
    let (d, kk) = (a + c, (1.0 - epsl) * kr);
    let disc = (b * b + 4.0 * d * kk).sqrt();
    let (r1, r2) = ((b + disc) / (2.0 * d), (b - disc) / (2.0 * d));
    // T = 1 + P e^{r1 (x - 1)} + R e^{r2 x}; T(0) = 0; T'(1) = 0
    let m = [[(-r1).exp(), 1.0], [r1, r2 * r2.exp()]];
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let (p, r) = ((-1.0 * m[1][1]) / det, (m[1][0] * 1.0) / det);
    let exact = move |x: f64| 1.0 + p * (r1 * (x - 1.0)).exp() + r * (r2 * x).exp();
    let mut prev = f64::NAN;
    for &h in &[1e2, 1e4, 1e6] {
        let col = column(0.025, epsl, kf, ks, h, av, kr, -kr).unwrap();
        let tf = &col.vars.scldom[col.temp].value;
        let ts = &col.vars.scldom[col.stemp].value;
        let dts = tf.iter().zip(ts).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        let (rms, emax) = errors(&col.vars, col.dom, col.temp, &exact);
        println!("RESULT lte_limit: h a_v = {:.1e}; max|T_s - T_f| = {dts:.3e}; T_f vs LTE exact rms = {rms:.3e} max = {emax:.3e}", h * av);
        if !prev.is_nan() { assert!(dts < 0.1 * prev); }
        prev = dts;
        if h == 1e6 { assert!(dts < 1e-4 && emax < 1e-2); }
    }
}

// ---------------------------------------------------------------- separate BCs for each phase

#[test]
fn per_phase_bcs() {
    // no flow; fluid: T_f(0) = 1, insulated at x = 1; solid: heat flux q_in into x = 0, T_s(1) = 0
    // A T_f'' + H (T_s - T_f) = 0; C T_s'' - H (T_s - T_f) = 0
    // S = A T_f + C T_s = al + be x; th = T_f - T_s = p cosh(m x) + r sinh(m x); m^2 = H (1 / A + 1 / C)
    let (epsl, kf, ks, h, av, qin): (f64, f64, f64, f64, f64, f64) = (0.4, 0.5, 2.0, 1.0, 2.0, 0.5);
    let (a, c, hh) = (epsl * kf, (1.0 - epsl) * ks, h * av);
    let m = (hh * (1.0 / a + 1.0 / c)).sqrt();
    let s = a + c;
    // unknowns (al, be, p, r); T_f = (S + C th) / s; T_s = (S - A th) / s
    let mat = [
        [1.0 / s, 0.0, c / s, 0.0],                                          // T_f(0) = 1
        [0.0, 1.0 / s, c / s * m * m.sinh(), c / s * m * m.cosh()],          // T_f'(1) = 0
        [0.0, c / s, 0.0, -a * c / s * m],                                   // outward flux C T_s'(0) = -q_in
        [1.0 / s, 1.0 / s, -a / s * m.cosh(), -a / s * m.sinh()],            // T_s(1) = 0
    ];
    let co = solve4(mat, [1.0, 0.0, -qin, 0.0]);
    let tf_ex = move |x: f64| (co[0] + co[1] * x + c * (co[2] * (m * x).cosh() + co[3] * (m * x).sinh())) / s;
    let ts_ex = move |x: f64| (co[0] + co[1] * x - a * (co[2] * (m * x).cosh() + co[3] * (m * x).sinh())) / s;
    println!("RESULT per_phase_bcs: T_f(1) = {:.6}; T_s(0) = {:.6}", tf_ex(1.0), ts_ex(0.0));

    let mut errs = vec![];
    let dxs = [0.05, 0.025];
    for &dx in &dxs {
        let Mesh { mut vars, dom, bnd } = mesh(1.0, 0.2, dx).unwrap();
        let [bnd_l, bnd_r, bnd_b, bnd_t] = bnd;
        let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "").unwrap();
        let pres = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
        let one = vars.add_scldom_cnst(dom, 1.0, "").unwrap();
        let zero = vars.add_scldom_cnst(dom, 0.0, "").unwrap();
        let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "").unwrap();
        let epsl_id = vars.add_scldom_cnst(dom, epsl, "").unwrap();
        let perm = vars.add_scldom_cnst(dom, 1e-2, "").unwrap();
        let temp = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
        let stemp = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
        let cond = vars.add_scldom_cnst(dom, kf, "").unwrap();
        let scond = vars.add_scldom_cnst(dom, ks, "").unwrap();
        let htc = vars.add_scldom_cnst(dom, h, "").unwrap();
        let area = vars.add_scldom_cnst(dom, av, "").unwrap();
        let vel0: Vec<usize> = bnd.iter().map(|&b| vars.add_vecbnd_cnst(b, (0.0, 0.0), "")).collect::<Result<_, _>>().unwrap();
        let flx0: Vec<usize> = bnd.iter().map(|&b| vars.add_sclbnd_cnst(b, 0.0, "")).collect::<Result<_, _>>().unwrap();
        let temp_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "").unwrap();
        let stemp_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "").unwrap();
        let shflx_l = vars.add_sclbnd_cnst(bnd_l, -qin, "").unwrap();  // outward flux; negative for inflow
        let mut phys = SteadyHeatPore::new().unwrap();
        phys.add_pore_dom(dom, vel, pres, one, one, fce, epsl_id, perm, zero).unwrap();
        phys.add_heatl_dom(dom, temp, one, cond, zero, zero).unwrap();
        phys.add_heats_dom(dom, stemp, one, scond, zero, zero, zero, zero).unwrap();
        phys.add_matrix_dom(dom, htc, area, zero, zero).unwrap();
        phys.add_templ_bnd(bnd_l, temp_l).unwrap(); phys.add_hflxl_bnd(bnd_r, flx0[1]).unwrap(); phys.add_hflxl_bnd(bnd_b, flx0[2]).unwrap(); phys.add_hflxl_bnd(bnd_t, flx0[3]).unwrap();
        phys.add_hflxs_bnd(bnd_l, shflx_l).unwrap(); phys.add_temps_bnd(bnd_r, stemp_r).unwrap(); phys.add_hflxs_bnd(bnd_b, flx0[2]).unwrap(); phys.add_hflxs_bnd(bnd_t, flx0[3]).unwrap();
        for i in 0..4 { phys.add_vel_bnd(bnd[i], vel0[i]).unwrap(); }
        phys.set_pres_ref(dom, 0, 0.0).unwrap();
        let mut lin = FaerLu::new(1).unwrap();
        phys.solve(&mut vars, &mut lin, 20, 1e-10, 1.0).unwrap();
        let ef = errors(&vars, dom, temp, &tf_ex);
        let es = errors(&vars, dom, stemp, &ts_ex);
        println!("RESULT per_phase_bcs: dx = {dx:.4}; fluid rms = {:.3e} max = {:.3e}; solid rms = {:.3e} max = {:.3e}", ef.0, ef.1, es.0, es.1);
        errs.push((ef, es));
    }
    assert!(errs[1].0.1 < 5e-3 && errs[1].1.1 < 5e-3);
    assert!(errs[1].0.0 < errs[0].0.0 && errs[1].1.0 < errs[0].1.0);
}

// ---------------------------------------------------------------- closed transient energy balance

// closed box with no flow and adiabatic walls for both phases
// E = sum_v [ epsl (rho cp)_f T_f + (1 - epsl) (rho c)_s T_s ] A_v is conserved by the exchange and diffusion
// run 1: variable epsl (with an epsl = 1 corner), non-uniform initial temperatures, no sources; E is constant
// run 2: constant epsl with uniform sources in both phases and on the surface; E grows by t * (epsl q_f + (1 - epsl) q_s + a_v q_a) * area
fn closed_box(sources: bool) -> (f64, f64, f64, f64) {
    let Mesh { mut vars, dom, bnd } = mesh(0.1, 0.05, 4e-3).unwrap();
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "").unwrap();
    let pres = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
    let den = vars.add_scldom_cnst(dom, 1.2, "").unwrap();
    let visc = vars.add_scldom_cnst(dom, 1.8e-5, "").unwrap();
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "").unwrap();
    let epsl = if sources { vars.add_scldom_cnst(dom, 0.4, "").unwrap() } else {
        vars.add_scldom_func(dom, Box::new(|_t, (x, y), _u| (0.35 + 0.5 * (x / 0.1) * (y / 0.05) + 0.3 * (x / 0.1).powi(2)).min(1.0)), vec![], "").unwrap()
    };
    let perm = vars.add_scldom_cnst(dom, 1e-8, "").unwrap();
    let zero = vars.add_scldom_cnst(dom, 0.0, "").unwrap();
    let temp = vars.add_scldom_unkn(dom, 300.0, "").unwrap();
    let stemp = vars.add_scldom_unkn(dom, 300.0, "").unwrap();
    let vlcp = vars.add_scldom_cnst(dom, 1.2e3, "").unwrap();
    let cond = vars.add_scldom_cnst(dom, 0.026, "").unwrap();
    let svlcp = vars.add_scldom_cnst(dom, 2.0e6, "").unwrap();
    let scond = vars.add_scldom_func(dom, Box::new(|_t, (x, _y), _u| 1.0 + 10.0 * x), vec![], "").unwrap();
    let htc = vars.add_scldom_cnst(dom, 50.0, "").unwrap();
    let area = vars.add_scldom_cnst(dom, 600.0, "").unwrap();
    let (qf, qs, qa) = if sources { (1e5, 2e5, 50.0) } else { (0.0, 0.0, 0.0) };
    let hsrc = vars.add_scldom_cnst(dom, qf, "").unwrap();
    let ssrc = vars.add_scldom_cnst(dom, qs, "").unwrap();
    let asrc = vars.add_scldom_cnst(dom, qa, "").unwrap();
    let vel0: Vec<usize> = bnd.iter().map(|&b| vars.add_vecbnd_cnst(b, (0.0, 0.0), "")).collect::<Result<_, _>>().unwrap();
    let flx0: Vec<usize> = bnd.iter().map(|&b| vars.add_sclbnd_cnst(b, 0.0, "")).collect::<Result<_, _>>().unwrap();

    // non-uniform initial temperatures
    if !sources {
        let pos = positions(&vars, dom);
        for (v, p) in pos.iter().enumerate() {
            vars.scldom[temp].value[v] = 300.0 + 100.0 * (p.0 / 0.1);
            vars.scldom[stemp].value[v] = 350.0 - 50.0 * (p.1 / 0.05);
        }
    }

    let mut phys = TransientHeatPore::new().unwrap();
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, zero).unwrap();
    phys.add_heatl_dom(dom, temp, vlcp, cond, zero, hsrc).unwrap();
    phys.add_heats_dom(dom, stemp, svlcp, scond, zero, ssrc, zero, asrc).unwrap();
    phys.add_matrix_dom(dom, htc, area, zero, zero).unwrap();
    for i in 0..4 {
        phys.add_hflxl_bnd(bnd[i], flx0[i]).unwrap();
        phys.add_hflxs_bnd(bnd[i], flx0[i]).unwrap();
        phys.add_vel_bnd(bnd[i], vel0[i]).unwrap();
    }
    phys.set_pres_ref(dom, 0, 0.0).unwrap();

    // E from the dual control volumes
    let energy = |vars: &Variables| -> f64 {
        let tri = &vars.tri[dom];
        let e = &vars.scldom[epsl].value;
        let (tf, ts) = (&vars.scldom[temp].value, &vars.scldom[stemp].value);
        (0..tf.len()).map(|v| (e[v] * 1.2e3 * tf[v] + (1.0 - e[v]) * 2.0e6 * ts[v]) * tri.cvol_area[v]).sum()
    };
    let (dt, nts) = (5.0, 20);
    let mut lin = FaerLu::new(1).unwrap();
    TransientBase::solve(&mut phys, &mut vars, &mut lin, dt, 1, 100, 20, 1e-12, 1.0).unwrap();
    // E needs the control volume areas and the epsl values, which only exist after a solve
    // the first one-step solve builds them; value_prev holds the initial state, which is restored before measuring E0
    let tf0 = vars.scldom[temp].value_prev.clone();
    let ts0 = vars.scldom[stemp].value_prev.clone();
    vars.scldom[temp].value = tf0;
    vars.scldom[stemp].value = ts0;
    let e0 = energy(&vars);
    TransientBase::solve(&mut phys, &mut vars, &mut lin, dt, nts, 100, 20, 1e-12, 1.0).unwrap();
    let e_end = energy(&vars);
    let tot_area: f64 = vars.tri[dom].cvol_area.iter().sum();
    let inj = (nts as f64) * dt * (0.4 * qf + 0.6 * qs + 600.0 * qa) * tot_area;
    let dts = vars.scldom[temp].value.iter().zip(&vars.scldom[stemp].value).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    (e0, e_end, inj, dts)
}

#[test]
fn closed_energy() {
    let (e0, e1, _, dts) = closed_box(false);
    println!("RESULT closed_energy (no sources): E0 = {e0:.10e}; E_end = {e1:.10e}; rel change = {:.3e}; max|T_s - T_f| at end = {dts:.3e}", (e1 - e0) / e0);
    assert!(((e1 - e0) / e0).abs() < 1e-10);
    let (e0, e1, inj, _) = closed_box(true);
    println!("RESULT closed_energy (sources): dE = {:.10e}; injected = {inj:.10e}; rel residual = {:.3e}", e1 - e0, (e1 - e0 - inj) / inj);
    assert!(((e1 - e0 - inj) / inj).abs() < 1e-8);
}

// ---------------------------------------------------------------- two layers joined by an interface

// [0, 0.2] x [0, 1] split at y = 1/2 into dom_p (bottom) and dom_f (top); plug flow U in y
// the interface is created from dom_p to dom_f, or from dom_f to dom_p if swap
struct Layered { vars: Variables, dom: [usize; 2], bnd: [usize; 6], itf: usize }

fn layered(dx: f64, swap: bool) -> Result<Layered, VFVChemError> {
    let mut vars = Variables::new_obj("tests/obj/obj_layered.obj", dx, dx, dx, 100, 42)?;
    let dom_p = vars.add_dom("dom_p")?;
    let dom_f = vars.add_dom("dom_f")?;
    let bnd_b = vars.add_bnd(dom_p, "bnd_b")?;
    let bnd_lp = vars.add_bnd(dom_p, "bnd_lp")?;
    let bnd_rp = vars.add_bnd(dom_p, "bnd_rp")?;
    let bnd_t = vars.add_bnd(dom_f, "bnd_t")?;
    let bnd_lf = vars.add_bnd(dom_f, "bnd_lf")?;
    let bnd_rf = vars.add_bnd(dom_f, "bnd_rf")?;
    let itf = if swap { vars.add_itf(dom_f, dom_p, "itf_m")? } else { vars.add_itf(dom_p, dom_f, "itf_m")? };
    Ok(Layered { vars, dom: [dom_p, dom_f], bnd: [bnd_b, bnd_lp, bnd_rp, bnd_t, bnd_lf, bnd_rf], itf })
}

// flow data of one layer: plug flow U in y; vel BC on the inlet and the side walls, p = 0 on the outlet
// returns (vel, pres, den, visc, fce, epsl, perm, forch)
fn layered_flow(vars: &mut Variables, dom: usize, u: f64, epsl: f64, perm: f64, forch: f64) -> Result<[usize; 8], VFVChemError> {
    Ok([
        vars.add_vecdom_unkn(dom, (0.0, u), "")?,
        vars.add_scldom_unkn(dom, 0.0, "")?,
        vars.add_scldom_cnst(dom, 1.0, "")?,
        vars.add_scldom_cnst(dom, 1.0, "")?,
        vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?,
        vars.add_scldom_cnst(dom, epsl, "")?,
        vars.add_scldom_cnst(dom, perm, "")?,
        vars.add_scldom_cnst(dom, forch, "")?,
    ])
}

#[test]
fn interface_split() {
    // the 1D LTNE column of ltne_column along y, split into two identical layers joined by fluid and solid interfaces
    let (u, epsl, kf, ks, h, av, qs) = (1.0, 0.4, 0.05, 0.5, 2.0, 2.5, 10.0);
    let (a, b, c, hh, q) = (epsl * kf, 1.0, (1.0 - epsl) * ks, h * av, (1.0 - epsl) * qs);
    let (tf_ex, ts_ex) = ltne_exact(a, b, c, hh, q);
    let dxs = [0.05, 0.025, 0.0125];
    let mut errs = vec![];
    for &dx in &dxs {
        let Layered { mut vars, dom, bnd, itf } = layered(dx, false).unwrap();
        let [bnd_b, bnd_lp, bnd_rp, bnd_t, bnd_lf, bnd_rf] = bnd;
        let mut phys = SteadyHeatPore::new().unwrap();
        let mut temp = vec![];
        for d in dom {
            let [vel, pres, den, visc, fce, epsl_id, perm, forch] = layered_flow(&mut vars, d, u, epsl, 1e-2, 0.5).unwrap();
            let templ = vars.add_scldom_unkn(d, 0.0, "").unwrap();
            let temps = vars.add_scldom_unkn(d, 0.0, "").unwrap();
            let one = vars.add_scldom_cnst(d, 1.0, "").unwrap();
            let zero = vars.add_scldom_cnst(d, 0.0, "").unwrap();
            let condl = vars.add_scldom_cnst(d, kf, "").unwrap();
            let conds = vars.add_scldom_cnst(d, ks, "").unwrap();
            let hsrcs = vars.add_scldom_cnst(d, qs, "").unwrap();
            let htc = vars.add_scldom_cnst(d, h, "").unwrap();
            let area = vars.add_scldom_cnst(d, av, "").unwrap();
            phys.add_pore_dom(d, vel, pres, den, visc, fce, epsl_id, perm, forch).unwrap();
            phys.add_heatl_dom(d, templ, one, condl, zero, zero).unwrap();
            phys.add_heats_dom(d, temps, one, conds, zero, hsrcs, zero, zero).unwrap();
            phys.add_matrix_dom(d, htc, area, zero, zero).unwrap();
            temp.push((templ, temps));
        }
        for bnd in [bnd_b, bnd_lp, bnd_rp, bnd_lf, bnd_rf] {
            let vel = vars.add_vecbnd_cnst(bnd, (0.0, u), "").unwrap();
            phys.add_vel_bnd(bnd, vel).unwrap();
        }
        let pres_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "").unwrap();
        phys.add_pres_bnd(bnd_t, pres_t).unwrap();
        let templ_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "").unwrap();
        phys.add_templ_bnd(bnd_b, templ_b).unwrap();
        phys.add_houtl_bnd(bnd_t).unwrap();
        for bnd in [bnd_lp, bnd_rp, bnd_lf, bnd_rf] {
            let hflxl = vars.add_sclbnd_cnst(bnd, 0.0, "").unwrap();
            phys.add_hflxl_bnd(bnd, hflxl).unwrap();
        }
        for bnd in bnd {
            let hflxs = vars.add_sclbnd_cnst(bnd, 0.0, "").unwrap();
            phys.add_hflxs_bnd(bnd, hflxs).unwrap();
        }
        phys.add_fcnt_itf(itf).unwrap();
        phys.add_hcntl_itf(itf).unwrap();
        phys.add_hcnts_itf(itf).unwrap();
        let mut lin = FaerLu::new(1).unwrap();
        phys.solve(&mut vars, &mut lin, 100, 1e-10, 1.0).unwrap();

        // errors over both layers against the exact solution in y
        let (mut sf, mut ss, mut mf, mut ms, mut num): (f64, f64, f64, f64, usize) = (0.0, 0.0, 0.0, 0.0, 0);
        for (d, &(templ, temps)) in dom.iter().zip(&temp) {
            for (v, p) in positions(&vars, *d).iter().enumerate() {
                let (ef, es) = (vars.scldom[templ].value[v] - tf_ex(p.1), vars.scldom[temps].value[v] - ts_ex(p.1));
                sf += ef * ef; ss += es * es;
                mf = mf.max(ef.abs()); ms = ms.max(es.abs());
                num += 1;
            }
        }
        errs.push((((sf / num as f64).sqrt(), mf), ((ss / num as f64).sqrt(), ms)));
    }
    for i in 0..dxs.len() {
        println!("RESULT interface_split: dx = {:.4}; fluid rms = {:.3e} max = {:.3e}; solid rms = {:.3e} max = {:.3e}", dxs[i], errs[i].0.0, errs[i].0.1, errs[i].1.0, errs[i].1.1);
    }
    let tscale = ts_ex(1.0);
    assert!(errs[2].0.0 < 0.5 * errs[0].0.0 && errs[2].1.0 < 0.5 * errs[0].1.0);
    assert!(errs[2].0.1 < 1e-2 * tscale && errs[2].1.1 < 1e-2 * tscale);
}

// packed bed with a heated solid (dom_p) under a clear fluid without a solid (dom_f)
// prescribed non-advective inflow q_in at the inlet, outflow at the outlet; adiabatic walls; adiabatic solid surface on the interface
// fluid sink s_f T in both layers; without it, T is only defined up to a constant (no Dirichlet BC)
// returns the global balance residual relative to the inflow plus generation, and the fluid temperatures of both layers
fn bed_fluid(swap: bool) -> (f64, Vec<f64>, Vec<f64>) {
    let (dx, u, epsl, kf, ks, h, av, qs, q_in, sf) = (0.025, 1.0, 0.4, 0.05, 0.5, 2.0, 2.5, 10.0, 0.5, -0.5);
    let Layered { mut vars, dom: [dom_p, dom_f], bnd, itf } = layered(dx, swap).unwrap();
    let [bnd_b, bnd_lp, bnd_rp, bnd_t, bnd_lf, bnd_rf] = bnd;
    let mut phys = SteadyHeatPore::new().unwrap();
    let mut templs = vec![];
    for (d, epsl_d, perm, forch) in [(dom_p, epsl, 1e-2, 0.5), (dom_f, 1.0, f64::INFINITY, 0.0)] {
        let [vel, pres, den, visc, fce, epsl_id, perm, forch] = layered_flow(&mut vars, d, u, epsl_d, perm, forch).unwrap();
        vars.set_scldom_scale(pres, 1.0).unwrap();  // rho U^2; the clear fluid pressure is zero (no flow resistance)
        let templ = vars.add_scldom_unkn(d, 0.0, "").unwrap();
        let one = vars.add_scldom_cnst(d, 1.0, "").unwrap();
        let zero = vars.add_scldom_cnst(d, 0.0, "").unwrap();
        let condl = vars.add_scldom_cnst(d, kf, "").unwrap();
        let hprdl = vars.add_scldom_cnst(d, sf, "").unwrap();
        let htc = vars.add_scldom_cnst(d, h, "").unwrap();
        let area = vars.add_scldom_cnst(d, if d == dom_p { av } else { 0.0 }, "").unwrap();
        phys.add_pore_dom(d, vel, pres, den, visc, fce, epsl_id, perm, forch).unwrap();
        phys.add_heatl_dom(d, templ, one, condl, hprdl, zero).unwrap();
        phys.add_matrix_dom(d, htc, area, zero, zero).unwrap();
        if d == dom_p {
            let temps = vars.add_scldom_unkn(d, 0.0, "").unwrap();
            let conds = vars.add_scldom_cnst(d, ks, "").unwrap();
            let hsrcs = vars.add_scldom_cnst(d, qs, "").unwrap();
            phys.add_heats_dom(d, temps, one, conds, zero, hsrcs, zero, zero).unwrap();
        }
        templs.push(templ);
    }
    for bnd in [bnd_b, bnd_lp, bnd_rp, bnd_lf, bnd_rf] {
        let vel = vars.add_vecbnd_cnst(bnd, (0.0, u), "").unwrap();
        phys.add_vel_bnd(bnd, vel).unwrap();
    }
    let pres_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "").unwrap();
    phys.add_pres_bnd(bnd_t, pres_t).unwrap();
    let hflxl_b = vars.add_sclbnd_cnst(bnd_b, -q_in, "").unwrap();  // outward flux; negative for inflow
    phys.add_hflxl_bnd(bnd_b, hflxl_b).unwrap();
    phys.add_houtl_bnd(bnd_t).unwrap();
    for bnd in [bnd_lp, bnd_rp, bnd_lf, bnd_rf] {
        let hflxl = vars.add_sclbnd_cnst(bnd, 0.0, "").unwrap();
        phys.add_hflxl_bnd(bnd, hflxl).unwrap();
    }
    for bnd in [bnd_b, bnd_lp, bnd_rp] {
        let hflxs = vars.add_sclbnd_cnst(bnd, 0.0, "").unwrap();
        phys.add_hflxs_bnd(bnd, hflxs).unwrap();
    }
    phys.add_fcnt_itf(itf).unwrap();
    phys.add_hcntl_itf(itf).unwrap();
    let mut lin = FaerLu::new(1).unwrap();
    phys.solve(&mut vars, &mut lin, 100, 1e-10, 1.0).unwrap();

    // outward advective flux of each control volume; T_v times the net outward face mass flux ((rho cp)_f = rho = 1)
    // the interior faces and the interface halves cancel in pairs, so the sum is the advective outflow minus inflow
    // sink as assembled: epsl s_f (T_v + grad T . skew) dA over all control volumes
    let (mut adv_net, mut sink) = (0.0, 0.0);
    for ((d, &templ), epsl_d) in [dom_p, dom_f].iter().zip(&templs).zip([epsl, 1.0]) {
        let flow = &vars.sclflw[phys.itr_flw[d]];
        let (t, grad) = (&vars.scldom[templ].value, &vars.scldom[templ].grad);
        let tri = &vars.tri[*d];
        adv_net += (0..t.len()).map(|v| -t[v] * flow.flux[v].iter().sum::<f64>()).sum::<f64>();
        sink += (0..t.len()).map(|v| epsl_d * sf * (t[v] + grad[v].0 * tri.cvol_skew[v].0 + grad[v].1 * tri.cvol_skew[v].1) * tri.cvol_area[v]).sum::<f64>();
    }
    let tri = &vars.tri[dom_p];
    let prod = (1.0 - epsl) * qs * tri.cvol_area.iter().sum::<f64>();
    let hflx_in = q_in * vars.bnd[bnd_b].nds_gnd.iter().map(|g| tri.bnd_half[vars.dom[dom_p].num_pts() + vars.dom[dom_p].nds_nid[g]].iter().sum::<f64>()).sum::<f64>();
    let resid = (hflx_in + prod + sink - adv_net) / (hflx_in + prod);
    println!("RESULT interface_adiabatic (swap = {swap}): inflow = {hflx_in:.6e}; generation = {prod:.6e}; sink = {sink:.6e}; advective net outflow = {adv_net:.6e}; rel residual = {resid:.3e}");
    (resid, vars.scldom[templs[0]].value.clone(), vars.scldom[templs[1]].value.clone())
}

#[test]
fn interface_adiabatic() {
    let (resid, tp, tf) = bed_fluid(false);
    let (resid_s, tp_s, tf_s) = bed_fluid(true);
    let tmax = tf.iter().cloned().fold(0.0, f64::max);
    let dt = tp.iter().zip(&tp_s).chain(tf.iter().zip(&tf_s)).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    println!("RESULT interface_adiabatic: T_f max = {tmax:.6e}; max|dT_f| after swapping the sides = {:.3e}", dt / tmax);
    assert!(resid.abs() < 1e-8 && resid_s.abs() < 1e-8);
    assert!(dt / tmax < 1e-8);
}
