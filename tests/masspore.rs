use std::collections::HashMap;
use vfvchem_vfvm2d::*;

// mass transfer in porous media (SteadyMassPore, TransientMassPore)
// epsl * rho * dc/dt + rho * vel . grad(c) = div( rho * ((epsl / tau) * D + D_disp) * grad(c) ) + epsl * R_v + a_v * R_s
// run with: cargo test --release --test masspore -- --nocapture --test-threads 1

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

// erf; Abramowitz and Stegun 7.1.26 is too coarse, so use the series / continued fraction via erfc
fn erf(x: f64) -> f64 {
    // W. J. Cody rational approximation is overkill here; Taylor series for |x| < 3, 1 beyond
    if x.abs() > 4.0 { return x.signum(); }
    let mut sum = x;
    let mut term = x;
    for n in 1..200 {
        term *= -x * x / n as f64;
        sum += term / (2 * n + 1) as f64;
        if term.abs() < 1e-17 { break; }
    }
    2.0 / std::f64::consts::PI.sqrt() * sum
}

// ---------------------------------------------------------------- epsl = 1 regression

// Poiseuille channel (Re = 50) with two components, cross diffusion, first-order production, and a source
// massflow solver and masspore solver with epsl = 1, tau = 1, alpha_L = alpha_T = 0, a_v = 0
fn poiseuille(pore: bool, transient: bool) -> Result<Vec<Vec<f64>>, VFVChemError> {
    let Mesh { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t] } = mesh(0.02, 0.005, 5e-4)?;
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 0.001, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    let epsl = vars.add_scldom_cnst(dom, 1.0, "")?;
    let perm = vars.add_scldom_cnst(dom, f64::INFINITY, "")?;  // free fluid
    let forch = vars.add_scldom_cnst(dom, 0.0, "")?;
    let tort = vars.add_scldom_cnst(dom, 1.0, "")?;
    let area = vars.add_scldom_cnst(dom, 0.0, "")?;
    let displ = vars.add_scldom_cnst(dom, 0.0, "")?;
    let dispt = vars.add_scldom_cnst(dom, 0.0, "")?;
    let conc: Vec<usize> = (0..2).map(|_| vars.add_scldom_unkn(dom, 0.0, "")).collect::<Result<_, _>>()?;
    // rho * D; D_00 = 2e-7, D_01 = 5e-8, D_11 = 1e-7; spatially varying D_00 exercises the skew correction
    let diff_00 = vars.add_scldom_func(dom, Box::new(|_t, (x, _y), _u| 1000.0 * 2e-7 * (1.0 + 20.0 * x)), vec![], "")?;
    let diff_01 = vars.add_scldom_cnst(dom, 1000.0 * 5e-8, "")?;
    let diff_11 = vars.add_scldom_cnst(dom, 1000.0 * 1e-7, "")?;
    let mprd = [vars.add_scldom_cnst(dom, -1000.0 * 0.5, "")?, vars.add_scldom_cnst(dom, 0.0, "")?];
    let msrc = [vars.add_scldom_cnst(dom, 0.0, "")?, vars.add_scldom_func(dom, Box::new(|_t, (x, y), _u| 1000.0 * 0.1 * (x / 0.02) * (y / 0.005)), vec![], "")?];
    let sprd = vars.add_scldom_cnst(dom, -7.0, "")?;  // ignored at a_v = 0
    let ssrc = vars.add_scldom_cnst(dom, 3.0, "")?;   // ignored at a_v = 0
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _u| (4.0 * 0.01 * y * (0.005 - y) / 0.005f64.powi(2), 0.0));
    let vel_l = vars.add_vecbnd_func(bnd_l, vel_func, vec![], "")?;
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;
    let conc_l = [vars.add_sclbnd_cnst(bnd_l, 1.0, "")?, vars.add_sclbnd_func(bnd_l, Box::new(|_t, (_x, y), _u| y / 0.005), vec![], "")?];
    let mflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let mflx_t = vars.add_sclbnd_cnst(bnd_t, 1e-4, "")?;
    let diffs = [HashMap::from([(0, diff_00), (1, diff_01)]), HashMap::from([(1, diff_11)])];
    let mut lin = FaerLu::new(1)?;
    macro_rules! setup {
        ($phys:ident) => {
            for i in 0..2 {
                $phys.add_conc_bnd(i, bnd_l, conc_l[i])?;
                $phys.add_mout_bnd(i, bnd_r)?;
                $phys.add_mflx_bnd(i, bnd_b, mflx_b)?;
                $phys.add_mflx_bnd(i, bnd_t, mflx_t)?;
            }
            $phys.add_vel_bnd(bnd_l, vel_l)?; $phys.add_pres_bnd(bnd_r, pres_r)?; $phys.add_vel_bnd(bnd_b, vel_b)?; $phys.add_vel_bnd(bnd_t, vel_t)?;
        };
    }
    match (pore, transient) {
        (false, false) => {
            let mut phys = SteadyMassFlow::new(2)?;
            for i in 0..2 { phys.add_mass_dom(i, dom, conc[i], diffs[i].clone(), mprd[i], msrc[i])?; }
            phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
            setup!(phys);
            SteadyBase::solve(&mut phys, &mut vars, &mut lin, 50, 1e-10, 1.0)?;
        }
        (true, false) => {
            let mut phys = SteadyMassPore::new(2)?;
            phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, forch)?;
            phys.add_matrix_dom(dom, tort, area, displ, dispt)?;
            for i in 0..2 { phys.add_mass_dom(i, dom, conc[i], diffs[i].clone(), mprd[i], msrc[i], sprd, ssrc)?; }
            setup!(phys);
            SteadyBase::solve(&mut phys, &mut vars, &mut lin, 50, 1e-10, 1.0)?;
        }
        (false, true) => {
            let mut phys = TransientMassFlow::new(2)?;
            for i in 0..2 { phys.add_mass_dom(i, dom, conc[i], diffs[i].clone(), mprd[i], msrc[i])?; }
            phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
            setup!(phys);
            TransientBase::solve(&mut phys, &mut vars, &mut lin, 0.05, 4, 100, 50, 1e-10, 1.0)?;
        }
        (true, true) => {
            let mut phys = TransientMassPore::new(2)?;
            phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, forch)?;
            phys.add_matrix_dom(dom, tort, area, displ, dispt)?;
            for i in 0..2 { phys.add_mass_dom(i, dom, conc[i], diffs[i].clone(), mprd[i], msrc[i], sprd, ssrc)?; }
            setup!(phys);
            TransientBase::solve(&mut phys, &mut vars, &mut lin, 0.05, 4, 100, 50, 1e-10, 1.0)?;
        }
    }
    Ok(conc.iter().map(|&c| vars.scldom[c].value.clone()).collect())
}

fn regression(transient: bool) {
    let conc_f = poiseuille(false, transient).unwrap();
    let conc_p = poiseuille(true, transient).unwrap();
    for i in 0..2 {
        let cmax = conc_f[i].iter().map(|a| a.abs()).fold(0.0, f64::max);
        let dc = conc_f[i].iter().zip(&conc_p[i]).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        println!("RESULT regression transient={transient} comp {i}: c_max = {cmax:.4e}; max|dc|/c_max = {:.3e}", dc / cmax);
        assert!(dc / cmax < 1e-10);
    }
}

#[test]
fn regression_steady() { regression(false); }

#[test]
fn regression_transient() { regression(true); }

// ---------------------------------------------------------------- plug flow through a packed column

// plug flow (superficial u_in) on [0, l] x [0, h]; vel BC on left, top, bottom; p = 0 on the right
// the mass BCs are set by the caller through the closure
struct Column { vars: Variables, dom: usize, bnd: [usize; 4], vel: usize, pres: usize, den: usize, visc: usize, fce: usize, epsl: usize, perm: usize, forch: usize, vel_bnd: [usize; 3], pres_r: usize }

fn column(l: f64, h: f64, dx: f64, u_in: f64, epsl_c: f64) -> Result<Column, VFVChemError> {
    let Mesh { mut vars, dom, bnd } = mesh(l, h, dx)?;
    let vel = vars.add_vecdom_unkn(dom, (u_in, 0.0), "")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom, 1.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 1.0, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    let epsl = vars.add_scldom_cnst(dom, epsl_c, "")?;
    let perm = vars.add_scldom_cnst(dom, 1e-2, "")?;
    let forch = vars.add_scldom_cnst(dom, 0.5, "")?;
    let vel_l = vars.add_vecbnd_cnst(bnd[0], (u_in, 0.0), "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd[2], (u_in, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd[3], (u_in, 0.0), "")?;
    let pres_r = vars.add_sclbnd_cnst(bnd[1], 0.0, "")?;
    Ok(Column { vars, dom, bnd, vel, pres, den, visc, fce, epsl, perm, forch, vel_bnd: [vel_l, vel_b, vel_t], pres_r })
}

// ---------------------------------------------------------------- steady 1D convection-dispersion-reaction

#[test]
fn packed_column_cdr() {
    // U c' = D_ax c'' - a_v k_s c; D_ax = (epsl / tau) D + alpha_L U / epsl; c(0) = 1, c'(L) = 0 (Danckwerts outlet)
    // epsl = 0.4, tau = 2, D = 0.1, alpha_L = 0.02, alpha_T = 0.005; a_v = 6 (1 - epsl) / d_p with d_p = 0.1; k_s = 0.05
    // the heterogeneous production coefficient is s_s = -rho k_s; alpha_T and the cross term act only through the random mesh
    let (l, h, u, epsl, tort, diff, displ, dispt, dp, ks): (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) = (1.0, 0.2, 1.0, 0.4, 2.0, 0.1, 0.02, 0.005, 0.1, 0.05);
    let av = 6.0 * (1.0 - epsl) / dp;
    let dax = epsl / tort * diff + displ * u / epsl;
    let kr = av * ks;
    let disc = (u * u + 4.0 * dax * kr).sqrt();
    let (r1, r2) = ((u + disc) / (2.0 * dax), (u - disc) / (2.0 * dax));
    // c = A e^{r1 x} + B e^{r2 x}; A + B = 1; A r1 e^{r1 L} + B r2 e^{r2 L} = 0
    let bb = 1.0 / (1.0 - r2 * (r2 * l).exp() / (r1 * (r1 * l).exp()));
    let aa = 1.0 - bb;
    let exact = move |x: f64| aa * (r1 * x).exp() + bb * (r2 * x).exp();
    println!("RESULT packed_column_cdr: D_ax = {dax:.4}; Pe = {:.2}; Da = {:.2}; c_exact(L) = {:.6e}", u * l / dax, kr * l / u, exact(l));
    let dxs = [0.05, 0.025, 0.0125];
    let mut errs = vec![];
    for &dx in &dxs {
        let Column { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t], vel, pres, den, visc, fce, epsl: epsl_id, perm, forch, vel_bnd: [vel_l, vel_b, vel_t], pres_r } = column(l, h, dx, u, epsl).unwrap();
        let conc = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
        let diff_id = vars.add_scldom_cnst(dom, diff, "").unwrap();
        let zero = vars.add_scldom_cnst(dom, 0.0, "").unwrap();
        let sprd = vars.add_scldom_cnst(dom, -ks, "").unwrap();
        let tort_id = vars.add_scldom_cnst(dom, tort, "").unwrap();
        let area_id = vars.add_scldom_cnst(dom, av, "").unwrap();
        let displ_id = vars.add_scldom_cnst(dom, displ, "").unwrap();
        let dispt_id = vars.add_scldom_cnst(dom, dispt, "").unwrap();
        let conc_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "").unwrap();
        let mflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "").unwrap();
        let mflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "").unwrap();
        let mut phys = SteadyMassPore::new(1).unwrap();
        phys.add_mass_dom(0, dom, conc, HashMap::from([(0, diff_id)]), zero, zero, sprd, zero).unwrap();
        phys.add_matrix_dom(dom, tort_id, area_id, displ_id, dispt_id).unwrap();
        phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl_id, perm, forch).unwrap();
        phys.add_conc_bnd(0, bnd_l, conc_l).unwrap(); phys.add_mout_bnd(0, bnd_r).unwrap();
        phys.add_mflx_bnd(0, bnd_b, mflx_b).unwrap(); phys.add_mflx_bnd(0, bnd_t, mflx_t).unwrap();
        phys.add_vel_bnd(bnd_l, vel_l).unwrap(); phys.add_pres_bnd(bnd_r, pres_r).unwrap();
        phys.add_vel_bnd(bnd_b, vel_b).unwrap(); phys.add_vel_bnd(bnd_t, vel_t).unwrap();
        let mut lin = FaerLu::new(1).unwrap();
        phys.solve(&mut vars, &mut lin, 100, 1e-10, 1.0).unwrap();
        let pos = positions(&vars, dom);
        let c = &vars.scldom[conc].value;
        let n = pos.len() as f64;
        let rms = (pos.iter().zip(c).map(|(p, c)| (c - exact(p.0)).powi(2)).sum::<f64>() / n).sqrt();
        let emax = pos.iter().zip(c).map(|(p, c)| (c - exact(p.0)).abs()).fold(0.0, f64::max);
        errs.push((rms, emax));
    }
    for i in 0..dxs.len() {
        let p = if i > 0 { (errs[i - 1].0 / errs[i].0).ln() / (dxs[i - 1] / dxs[i]).ln() } else { f64::NAN };
        println!("RESULT packed_column_cdr: dx = {:.4}; rms = {:.3e}; max = {:.3e}; order(rms) = {:.2}", dxs[i], errs[i].0, errs[i].1, p);
    }
    assert!(errs[2].0 < 0.5 * errs[0].0 && errs[2].1 < 2e-2);
}

// ---------------------------------------------------------------- transient Gaussian pulse with longitudinal dispersion

// moments of c over the domain: (mass, mean x, variance in x)
fn moments(vars: &Variables, dom: usize, conc: usize) -> (f64, f64, f64) {
    let pos = positions(vars, dom);
    let area = &vars.tri[dom].cvol_area;
    let c = &vars.scldom[conc].value;
    let m0: f64 = (0..pos.len()).map(|i| c[i] * area[i]).sum();
    let m1: f64 = (0..pos.len()).map(|i| c[i] * area[i] * pos[i].0).sum::<f64>() / m0;
    let m2: f64 = (0..pos.len()).map(|i| c[i] * area[i] * (pos[i].0 - m1).powi(2)).sum::<f64>() / m0;
    (m0, m1, m2)
}

fn pulse(displ: f64, dt: f64) -> Result<(f64, f64, f64, f64, f64), VFVChemError> {
    // epsl c_t + U c_x = D_ax c_xx; D_ax = (epsl / tau) D + alpha_L U / epsl
    // pulse: mean x0 + (U / epsl) t; variance s0^2 + 2 (D_ax / epsl) t
    let (l, h, dx, u, epsl, tort, diff, x0, s0, t_end) = (1.0f64, 0.05, 0.01, 0.4, 0.4, 1.5, 1e-5, 0.25, 0.05, 0.4);
    let Column { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t], vel, pres, den, visc, fce, epsl: epsl_id, perm, forch, vel_bnd: [vel_l, vel_b, vel_t], pres_r } = column(l, h, dx, u, epsl)?;
    let conc = vars.add_scldom_unkn(dom, 0.0, "")?;
    let pos = positions(&vars, dom);
    for (i, p) in pos.iter().enumerate() { vars.scldom[conc].value[i] = (-(p.0 - x0).powi(2) / (2.0 * s0 * s0)).exp(); }
    let diff_id = vars.add_scldom_cnst(dom, diff, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let tort_id = vars.add_scldom_cnst(dom, tort, "")?;
    let displ_id = vars.add_scldom_cnst(dom, displ, "")?;
    let dispt_id = vars.add_scldom_cnst(dom, 0.1 * displ, "")?;
    let conc_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;
    let mflx = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let mflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;
    let mut phys = TransientMassPore::new(1)?;
    phys.add_mass_dom(0, dom, conc, HashMap::from([(0, diff_id)]), zero, zero, zero, zero)?;
    phys.add_matrix_dom(dom, tort_id, zero, displ_id, dispt_id)?;
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl_id, perm, forch)?;
    phys.add_conc_bnd(0, bnd_l, conc_l)?; phys.add_mout_bnd(0, bnd_r)?;
    phys.add_mflx_bnd(0, bnd_b, mflx)?; phys.add_mflx_bnd(0, bnd_t, mflx_t)?;
    phys.add_vel_bnd(bnd_l, vel_l)?; phys.add_pres_bnd(bnd_r, pres_r)?;
    phys.add_vel_bnd(bnd_b, vel_b)?; phys.add_vel_bnd(bnd_t, vel_t)?;
    let m0_init = s0 * (2.0 * std::f64::consts::PI).sqrt() * h;  // the CV areas are built in solve
    let mut lin = FaerLu::new(1)?;
    let num_ts = (t_end / dt).round() as usize;
    phys.solve(&mut vars, &mut lin, dt, num_ts, num_ts, 50, 1e-10, 1.0)?;
    let (m0, m1, m2) = moments(&vars, dom, conc);
    let dax = epsl / tort * diff + displ * u / epsl;
    let mean = x0 + u / epsl * t_end;
    let var = s0 * s0 + 2.0 * dax / epsl * t_end;
    // backward Euler adds v^2 dt / 2 to D_ax / epsl at first order
    let var_be = var + (u / epsl).powi(2) * dt * t_end;
    println!("RESULT pulse alpha_L = {displ}; dt = {dt}: mass ratio = {:.6}; mean = {m1:.5} (exact {mean:.5}); var = {m2:.5e} (exact {var:.5e}; with BE diffusion {var_be:.5e})", m0 / m0_init);
    Ok((m1 - mean, (m2 - var) / var, (m2 - var_be) / var_be, m0 / m0_init - 1.0, m2))
}

#[test]
fn transient_pulse_dispersion() {
    let (dm_a, dv_a, dvbe_a, dmass_a, _) = pulse(0.005, 0.002).unwrap();
    let (dm_b, dv_b, dvbe_b, _, var_b) = pulse(0.005, 0.001).unwrap();
    let (_, _, _, _, var_0) = pulse(0.0, 0.001).unwrap();
    // without alpha_L the cell Peclet number is ~ 400 and the power-law scheme falls back to upwinding (numerical D ~ v dx / 2)
    let var_up = 2.5053e-3 + 2.0 * (0.4 / 0.4) * 0.01 / 2.0 * 0.4 + (0.4f64 / 0.4).powi(2) * 0.001 * 0.4;
    println!("RESULT pulse: alpha_L = 0 variance = {var_0:.5e} (upwind estimate {var_up:.5e}); alpha_L = 0.005 variance = {var_b:.5e}");
    println!("RESULT pulse: variance rel err = {dv_a:.3e} (dt = 0.002), {dv_b:.3e} (dt = 0.001); after BE correction {dvbe_a:.3e}, {dvbe_b:.3e}; mean err {dm_a:.2e}, {dm_b:.2e}");
    // the mass reference is the exact integral of the initial Gaussian; the time error is first order; the remainder is the spatial error of the power-law scheme
    assert!(dmass_a.abs() < 5e-3 && dm_a.abs() < 5e-3 && dm_b.abs() < 5e-3);
    assert!(dvbe_b.abs() < 0.1 && dv_b.abs() < dv_a.abs());
}

// ---------------------------------------------------------------- steady transverse dispersion of a step

#[test]
fn transverse_dispersion() {
    // inlet step c = 1 for y > h / 2; downstream U c_x = D_T c_yy with D_T = (epsl / tau) D + alpha_T U / epsl
    // c = (1 + erf((y - h/2) / sqrt(4 D_T x / U))) / 2; alpha_L = 5 alpha_T so that a wrong tensor would show
    // the longitudinal term D_L c_xx is neglected in the reference (relative size ~ D_L / (U x))
    let (l, h, dx, u, epsl, tort, diff, displ, dispt) = (1.0f64, 0.6, 0.0125, 1.0, 0.8, 1.0, 1e-6, 0.02, 0.004);
    let Column { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t], vel, pres, den, visc, fce, epsl: epsl_id, perm, forch, vel_bnd: [vel_l, vel_b, vel_t], pres_r } = column(l, h, dx, u, epsl).unwrap();
    let conc = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
    let diff_id = vars.add_scldom_cnst(dom, diff, "").unwrap();
    let zero = vars.add_scldom_cnst(dom, 0.0, "").unwrap();
    let tort_id = vars.add_scldom_cnst(dom, tort, "").unwrap();
    let displ_id = vars.add_scldom_cnst(dom, displ, "").unwrap();
    let dispt_id = vars.add_scldom_cnst(dom, dispt, "").unwrap();
    let conc_l = vars.add_sclbnd_func(bnd_l, Box::new(move |_t, (_x, y), _u| if y > 0.5 * h { 1.0 } else { 0.0 }), vec![], "").unwrap();
    let mflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "").unwrap();
    let mflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "").unwrap();
    let mut phys = SteadyMassPore::new(1).unwrap();
    phys.add_mass_dom(0, dom, conc, HashMap::from([(0, diff_id)]), zero, zero, zero, zero).unwrap();
    phys.add_matrix_dom(dom, tort_id, zero, displ_id, dispt_id).unwrap();
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl_id, perm, forch).unwrap();
    phys.add_conc_bnd(0, bnd_l, conc_l).unwrap(); phys.add_mout_bnd(0, bnd_r).unwrap();
    phys.add_mflx_bnd(0, bnd_b, mflx_b).unwrap(); phys.add_mflx_bnd(0, bnd_t, mflx_t).unwrap();
    phys.add_vel_bnd(bnd_l, vel_l).unwrap(); phys.add_pres_bnd(bnd_r, pres_r).unwrap();
    phys.add_vel_bnd(bnd_b, vel_b).unwrap(); phys.add_vel_bnd(bnd_t, vel_t).unwrap();
    let mut lin = FaerLu::new(1).unwrap();
    phys.solve(&mut vars, &mut lin, 100, 1e-10, 1.0).unwrap();
    let pos = positions(&vars, dom);
    let c = &vars.scldom[conc].value;
    let dtr = epsl / tort * diff + dispt * u / epsl;
    let dlg = epsl / tort * diff + displ * u / epsl;
    // fit D from the slab 0.6 < x < 0.8: dc/dy at the centre line ~ 1 / sqrt(4 pi D x / U); least squares of c against erf over D
    let pts: Vec<(f64, f64, f64)> = pos.iter().zip(c).filter(|(p, _)| p.0 > 0.6 && p.0 < 0.8).map(|(p, &c)| (p.0, p.1, c)).collect();
    let resid = |d: f64| pts.iter().map(|&(x, y, c)| (c - 0.5 * (1.0 + erf((y - 0.5 * h) / (4.0 * d * x / u).sqrt()))).powi(2)).sum::<f64>();
    let (mut lo, mut hi) = (0.1 * dtr, 10.0 * dtr);
    for _ in 0..200 {
        let (m1, m2) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if resid(m1) < resid(m2) { hi = m2 } else { lo = m1 }
    }
    let d_fit = 0.5 * (lo + hi);
    let rms = (resid(dtr) / pts.len() as f64).sqrt();
    println!("RESULT transverse_dispersion: D_T = {dtr:.4e}; D_L = {dlg:.4e}; fitted D = {d_fit:.4e}; rel err = {:.3e}; rms vs erf(D_T) = {rms:.3e}", d_fit / dtr - 1.0);
    // includes the numerical crosswind diffusion of the upwinded faces on the random mesh
    assert!((d_fit / dtr - 1.0).abs() < 0.25);
}

// ---------------------------------------------------------------- global balance with dispersion at the outlet

#[test]
fn outflow_conservation() {
    // packed column with a prescribed non-advective inflow at the inlet, a surface reaction, and dispersion reaching the outlet
    // no node has a Dirichlet condition, so the sum of all node balances is the global balance
    // inflow (prescribed non-advective + advective) = outflow (advective) + consumption
    // the interior faces cancel in pairs only if the outlet control volumes carry the same diffusion and dispersion as their neighbors
    let (l, h, dx, u, epsl, tort, diff, displ, dispt, av, ks, flux_in) = (1.0, 0.2, 0.025, 1.0, 0.4, 1.5, 0.05, 0.1, 0.02, 3.0, 0.5, 1.0);
    let Column { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t], vel, pres, den, visc, fce, epsl: epsl_id, perm, forch, vel_bnd: [vel_l, vel_b, vel_t], pres_r } = column(l, h, dx, u, epsl).unwrap();
    let conc = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
    let diff_id = vars.add_scldom_cnst(dom, diff, "").unwrap();
    let zero = vars.add_scldom_cnst(dom, 0.0, "").unwrap();
    let sprd = vars.add_scldom_cnst(dom, -ks, "").unwrap();
    let tort_id = vars.add_scldom_cnst(dom, tort, "").unwrap();
    let area_id = vars.add_scldom_cnst(dom, av, "").unwrap();
    let displ_id = vars.add_scldom_cnst(dom, displ, "").unwrap();
    let dispt_id = vars.add_scldom_cnst(dom, dispt, "").unwrap();
    let mflx_l = vars.add_sclbnd_cnst(bnd_l, -flux_in, "").unwrap();  // outward flux; negative for inflow
    let mflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "").unwrap();
    let mflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "").unwrap();
    let mut phys = SteadyMassPore::new(1).unwrap();
    phys.add_mass_dom(0, dom, conc, HashMap::from([(0, diff_id)]), zero, zero, sprd, zero).unwrap();
    phys.add_matrix_dom(dom, tort_id, area_id, displ_id, dispt_id).unwrap();
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl_id, perm, forch).unwrap();
    phys.add_mflx_bnd(0, bnd_l, mflx_l).unwrap(); phys.add_mout_bnd(0, bnd_r).unwrap();
    phys.add_mflx_bnd(0, bnd_b, mflx_b).unwrap(); phys.add_mflx_bnd(0, bnd_t, mflx_t).unwrap();
    phys.add_vel_bnd(bnd_l, vel_l).unwrap(); phys.add_pres_bnd(bnd_r, pres_r).unwrap();
    phys.add_vel_bnd(bnd_b, vel_b).unwrap(); phys.add_vel_bnd(bnd_t, vel_t).unwrap();
    let mut lin = FaerLu::new(1).unwrap();
    phys.solve(&mut vars, &mut lin, 100, 1e-12, 1.0).unwrap();

    // outward advective flux of each control volume; c_v times the net outward face mass flux (rho = 1)
    let tri = &vars.tri[dom];
    let flow = &vars.sclflw[phys.itr_flw[&dom]];
    let c = &vars.scldom[conc];
    let adv: Vec<f64> = (0..c.value.len()).map(|v| -c.value[v] * flow.flux[v].iter().sum::<f64>()).collect();
    let nds_vid = |bnd: usize| -> Vec<usize> { vars.bnd[bnd].nds_gnd.iter().map(|g| vars.dom[dom].num_pts() + vars.dom[dom].nds_nid[g]).collect() };
    let adv_in = -nds_vid(bnd_l).iter().map(|&v| adv[v]).sum::<f64>();
    let adv_out = nds_vid(bnd_r).iter().map(|&v| adv[v]).sum::<f64>();
    let adv_net: f64 = adv.iter().sum();

    // prescribed inflow over the boundary half-edges of the inlet nodes; consumption with the centroid skewness correction
    let mflx_in = flux_in * nds_vid(bnd_l).iter().map(|&v| tri.bnd_half[v].iter().sum::<f64>()).sum::<f64>();
    let cons: f64 = (0..c.value.len()).map(|v| av * ks * (c.value[v] + c.grad[v].0 * tri.cvol_skew[v].0 + c.grad[v].1 * tri.cvol_skew[v].1) * tri.cvol_area[v]).sum();
    let resid = mflx_in - adv_net - cons;
    println!("RESULT outflow_conservation: inflow = {:.6e} (prescribed {mflx_in:.6e} + advective {adv_in:.6e}); outflow = {adv_out:.6e}; walls = {:.3e}; consumption = {cons:.6e}; rel residual = {:.3e}", mflx_in + adv_in, adv_net - adv_out + adv_in, resid / (mflx_in + adv_in));
    assert!((resid / (mflx_in + adv_in)).abs() < 1e-8);
}

// ---------------------------------------------------------------- two porous layers in series joined by an interface

// 1D reference: U c' = (D_ax c')' - k_r c on [0, 1]; c(0) = 1, c'(1) = 0; D_ax and k_r jump at y = a
// c and D_ax c' are continuous at y = a
// node-based conservative FD with central advection, mirrored ghost node at the outlet, Thomas algorithm
fn layered_cdr(n: usize, u: f64, a: f64, dax: [f64; 2], kr: [f64; 2]) -> Vec<f64> {
    let h = 1.0 / n as f64;
    let dax_y = |y: f64| if y < a { dax[0] } else { dax[1] };
    let kr_y = |y: f64| if y < a { kr[0] } else { kr[1] };
    let (mut lo, mut di, mut up, mut rhs) = (vec![0.0; n + 1], vec![1.0; n + 1], vec![0.0; n + 1], vec![0.0; n + 1]);
    rhs[0] = 1.0;
    for i in 1..=n {
        let y = i as f64 * h;
        let (dax_lo, dax_up) = (dax_y(y - 0.5 * h), if i < n { dax_y(y + 0.5 * h) } else { dax_y(y - 0.5 * h) });
        lo[i] = -dax_lo / (h * h) - u / (2.0 * h);
        up[i] = -dax_up / (h * h) + u / (2.0 * h);
        di[i] = (dax_lo + dax_up) / (h * h) + 0.5 * (kr_y(y - 0.25 * h) + kr_y(y + 0.25 * h));
    }
    lo[n] += up[n];
    for i in 1..=n {
        let m = lo[i] / di[i - 1];
        di[i] -= m * up[i - 1];
        rhs[i] -= m * rhs[i - 1];
    }
    let mut c = vec![0.0; n + 1];
    c[n] = rhs[n] / di[n];
    for i in (0..n).rev() { c[i] = (rhs[i] - up[i] * c[i + 1]) / di[i]; }
    c
}

#[test]
fn interface_layered() {
    // plug flow U in y through two packed layers in series on [0, 0.2] x [0, 1], joined by an interface at y = 1/2
    // the layers share epsl but differ in tau, a_v, and alpha_L; surface reaction s_s = -rho k_s
    // D_ax = (epsl / tau) D + alpha_L U / epsl, k_r = a_v k_s; c = 1 at the inlet (y = 0), outflow at y = 1
    let (u, epsl, diff, ks) = (1.0f64, 0.4, 0.1, 0.2);
    let (tort, av, displ) = ([1.5, 3.0], [6.0, 2.0], [0.02, 0.05]);
    let dax = [epsl / tort[0] * diff + displ[0] * u / epsl, epsl / tort[1] * diff + displ[1] * u / epsl];
    let ref_1d = layered_cdr(20000, u, 0.5, dax, [av[0] * ks, av[1] * ks]);
    let n_ref = (ref_1d.len() - 1) as f64;
    let exact = |y: f64| { let t = (y * n_ref).clamp(0.0, n_ref - 1e-9); let i = t.floor() as usize; ref_1d[i] + (t - i as f64) * (ref_1d[i + 1] - ref_1d[i]) };
    let dxs = [0.05, 0.025, 0.0125];
    let mut errs = vec![];
    for &dx in &dxs {
        let mut vars = Variables::new_obj("tests/obj/obj_layered.obj", dx, dx, dx, 100, 42).unwrap();
        let dom_p = vars.add_dom("dom_p").unwrap();
        let dom_f = vars.add_dom("dom_f").unwrap();
        let bnd_b = vars.add_bnd(dom_p, "bnd_b").unwrap();
        let bnd_lp = vars.add_bnd(dom_p, "bnd_lp").unwrap();
        let bnd_rp = vars.add_bnd(dom_p, "bnd_rp").unwrap();
        let bnd_t = vars.add_bnd(dom_f, "bnd_t").unwrap();
        let bnd_lf = vars.add_bnd(dom_f, "bnd_lf").unwrap();
        let bnd_rf = vars.add_bnd(dom_f, "bnd_rf").unwrap();
        let itf_m = vars.add_itf(dom_p, dom_f, "itf_m").unwrap();
        let mut phys = SteadyMassPore::new(1).unwrap();
        let mut concs = vec![];
        for (i, dom) in [dom_p, dom_f].into_iter().enumerate() {
            let vel = vars.add_vecdom_unkn(dom, (0.0, u), "").unwrap();
            let pres = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
            let den = vars.add_scldom_cnst(dom, 1.0, "").unwrap();
            let visc = vars.add_scldom_cnst(dom, 1.0, "").unwrap();
            let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "").unwrap();
            let epsl_id = vars.add_scldom_cnst(dom, epsl, "").unwrap();
            let perm = vars.add_scldom_cnst(dom, 1e-2, "").unwrap();
            let forch = vars.add_scldom_cnst(dom, 0.5, "").unwrap();
            let conc = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
            let diff_id = vars.add_scldom_cnst(dom, diff, "").unwrap();
            let zero = vars.add_scldom_cnst(dom, 0.0, "").unwrap();
            let sprd = vars.add_scldom_cnst(dom, -ks, "").unwrap();
            let tort_id = vars.add_scldom_cnst(dom, tort[i], "").unwrap();
            let area_id = vars.add_scldom_cnst(dom, av[i], "").unwrap();
            let displ_id = vars.add_scldom_cnst(dom, displ[i], "").unwrap();
            let dispt_id = vars.add_scldom_cnst(dom, 0.2 * displ[i], "").unwrap();
            phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl_id, perm, forch).unwrap();
            phys.add_matrix_dom(dom, tort_id, area_id, displ_id, dispt_id).unwrap();
            phys.add_mass_dom(0, dom, conc, HashMap::from([(0, diff_id)]), zero, zero, sprd, zero).unwrap();
            concs.push(conc);
        }
        for bnd in [bnd_b, bnd_lp, bnd_rp, bnd_lf, bnd_rf] {
            let vel = vars.add_vecbnd_cnst(bnd, (0.0, u), "").unwrap();
            phys.add_vel_bnd(bnd, vel).unwrap();
        }
        let pres_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "").unwrap();
        phys.add_pres_bnd(bnd_t, pres_t).unwrap();
        let conc_b = vars.add_sclbnd_cnst(bnd_b, 1.0, "").unwrap();
        phys.add_conc_bnd(0, bnd_b, conc_b).unwrap();
        phys.add_mout_bnd(0, bnd_t).unwrap();
        for bnd in [bnd_lp, bnd_rp, bnd_lf, bnd_rf] {
            let mflx = vars.add_sclbnd_cnst(bnd, 0.0, "").unwrap();
            phys.add_mflx_bnd(0, bnd, mflx).unwrap();
        }
        phys.add_fcnt_itf(itf_m).unwrap();
        phys.add_mcnt_itf(0, itf_m).unwrap();
        let mut lin = FaerLu::new(1).unwrap();
        phys.solve(&mut vars, &mut lin, 100, 1e-10, 1.0).unwrap();
        let mut sum = 0.0;
        let mut emax: f64 = 0.0;
        let mut num = 0;
        for (dom, conc) in [(dom_p, concs[0]), (dom_f, concs[1])] {
            for (p, c) in positions(&vars, dom).iter().zip(&vars.scldom[conc].value) {
                let e = c - exact(p.1);
                sum += e * e;
                emax = emax.max(e.abs());
                num += 1;
            }
        }
        errs.push(((sum / num as f64).sqrt(), emax));
    }
    for i in 0..dxs.len() {
        let p = if i > 0 { (errs[i - 1].0 / errs[i].0).ln() / (dxs[i - 1] / dxs[i]).ln() } else { f64::NAN };
        println!("RESULT interface_layered: dx = {:.4}; rms = {:.3e}; max = {:.3e}; order(rms) = {:.2}", dxs[i], errs[i].0, errs[i].1, p);
    }
    println!("RESULT interface_layered: D_ax = {:.4}, {:.4}; c_exact(1/2) = {:.6e}; c_exact(1) = {:.6e}", dax[0], dax[1], exact(0.5), exact(1.0));
    assert!(errs[2].0 < 0.5 * errs[0].0 && errs[2].1 < 2e-2);
}
