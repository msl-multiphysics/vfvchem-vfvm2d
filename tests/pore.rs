use vfvchem_vfvm2d::*;

// Darcy-Brinkman-Forchheimer solvers (SteadyPore, TransientPore)
// K and F come from the Ergun relation of the porosity (epsl) and particle diameter (diam); see PoreOpVecDomDrag::closure
// run with: cargo test --release --test pore -- --nocapture --test-threads 1

// positions of all vertices in a domain; vid -> (x, y)
fn positions(vars: &Variables, dom: usize) -> Vec<(f64, f64)> {
    let mut pos: Vec<(f64, f64)> = vars.dom[dom].pts_gpt.iter().map(|&g| vars.pts.pts_pos[g]).collect();
    pos.extend(vars.dom[dom].nds_gnd.iter().map(|&g| vars.nds.nds_pos[g]));
    pos
}

// int u dy over a vertical boundary from its node values (trapezoid); the walls at y = 0 and y = h have u = 0
fn end_flux(vars: &Variables, bnd: usize, vel: usize, h: f64) -> f64 {
    let dom = &vars.dom[vars.bnd[bnd].dom_id];
    let mut prof: Vec<(f64, f64)> = vars.bnd[bnd].nds_gnd.iter().map(|&g| (vars.nds.nds_pos[g].1, vars.vecdom[vel].value[dom.num_pts() + dom.nds_nid[&g]].0)).collect();
    prof.push((0.0, 0.0));
    prof.push((h, 0.0));
    prof.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    prof.windows(2).map(|w| 0.5 * (w[0].1 + w[1].1) * (w[1].0 - w[0].0)).sum()
}

// Ergun relation; returns (1 / K, F / sqrt(K))
fn ergun(epsl: f64, diam: f64) -> (f64, f64) {
    (150.0 * (1.0 - epsl).powi(2) / (diam * diam * epsl.powi(3)), 1.75 * (1.0 - epsl) / (diam * epsl.powi(3)))
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

// ---------------------------------------------------------------- epsl = 1 regression

// Poiseuille channel with parabolic inlet (Re = 50); flow solver and pore solver with epsl = 1 (no drag)
fn poiseuille(pore: bool, transient: bool) -> Result<(Vec<(f64, f64)>, Vec<f64>), VFVChemError> {
    let Mesh { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t] } = mesh(0.02, 0.005, 5e-4)?;
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 0.001, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    let epsl = vars.add_scldom_cnst(dom, 1.0, "")?;
    let diam = vars.add_scldom_cnst(dom, 1e-3, "")?;
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _u| (4.0 * 0.01 * y * (0.005 - y) / 0.005f64.powi(2), 0.0));
    let vel_l = vars.add_vecbnd_func(bnd_l, vel_func, vec![], "")?;
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;
    let mut lin = LinearLu::new(1)?;
    match (pore, transient) {
        (false, false) => {
            let mut phys = SteadyFlow::new()?;
            phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
            phys.add_vel_bnd(bnd_l, vel_l)?; phys.add_pres_bnd(bnd_r, pres_r)?; phys.add_vel_bnd(bnd_b, vel_b)?; phys.add_vel_bnd(bnd_t, vel_t)?;
            SteadyBase::solve(&mut phys, &mut vars, &mut lin, 50, 1e-10, 1.0)?;
        }
        (true, false) => {
            let mut phys = SteadyPore::new()?;
            phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, diam)?;
            phys.add_vel_bnd(bnd_l, vel_l)?; phys.add_pres_bnd(bnd_r, pres_r)?; phys.add_vel_bnd(bnd_b, vel_b)?; phys.add_vel_bnd(bnd_t, vel_t)?;
            SteadyBase::solve(&mut phys, &mut vars, &mut lin, 50, 1e-10, 1.0)?;
        }
        (false, true) => {
            let mut phys = TransientFlow::new()?;
            phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
            phys.add_vel_bnd(bnd_l, vel_l)?; phys.add_pres_bnd(bnd_r, pres_r)?; phys.add_vel_bnd(bnd_b, vel_b)?; phys.add_vel_bnd(bnd_t, vel_t)?;
            TransientBase::solve(&mut phys, &mut vars, &mut lin, 0.05, 4, 100, 50, 1e-10, 1.0)?;
        }
        (true, true) => {
            let mut phys = TransientPore::new()?;
            phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, diam)?;
            phys.add_vel_bnd(bnd_l, vel_l)?; phys.add_pres_bnd(bnd_r, pres_r)?; phys.add_vel_bnd(bnd_b, vel_b)?; phys.add_vel_bnd(bnd_t, vel_t)?;
            TransientBase::solve(&mut phys, &mut vars, &mut lin, 0.05, 4, 100, 50, 1e-10, 1.0)?;
        }
    }
    Ok((vars.vecdom[vel].value.clone(), vars.scldom[pres].value.clone()))
}

fn regression(transient: bool) {
    let (vel_f, pres_f) = poiseuille(false, transient).unwrap();
    let (vel_p, pres_p) = poiseuille(true, transient).unwrap();
    let du = vel_f.iter().zip(&vel_p).map(|(a, b)| (a.0 - b.0).abs().max((a.1 - b.1).abs())).fold(0.0, f64::max);
    let dp = pres_f.iter().zip(&pres_p).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    let pmax = pres_f.iter().map(|a| a.abs()).fold(0.0, f64::max);
    println!("RESULT regression transient={transient}: max|du|/u_max = {:.3e}; max|dp|/p_max = {:.3e}", du / 0.01, dp / pmax);
    assert!(du / 0.01 < 1e-10 && dp / pmax < 1e-10);
}

#[test]
fn regression_steady() { regression(false); }

#[test]
fn regression_transient() { regression(true); }

// ---------------------------------------------------------------- body-force channel, no-slip walls, open ends (p = 0)

// returns (positions, velocity, int u dy at the left end, int u dy at the right end)
// epsl and diam are functions of y; mu = rho = 1
fn channel(dx: f64, g: f64, epsl: fn(f64) -> f64, diam: fn(f64) -> f64) -> Result<(Vec<(f64, f64)>, Vec<(f64, f64)>, f64, f64), VFVChemError> {
    let Mesh { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t] } = mesh(0.2, 1.0, dx)?;
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom, 1.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 1.0, "")?;
    let fce = vars.add_vecdom_cnst(dom, (g, 0.0), "")?;
    let epsl = vars.add_scldom_func(dom, Box::new(move |_t, (_x, y), _u| epsl(y)), vec![], "")?;
    let diam = vars.add_scldom_func(dom, Box::new(move |_t, (_x, y), _u| diam(y)), vec![], "")?;
    let pres_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;
    let mut phys = SteadyPore::new()?;
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, diam)?;
    phys.add_pres_bnd(bnd_l, pres_l)?; phys.add_pres_bnd(bnd_r, pres_r)?;
    phys.add_vel_bnd(bnd_b, vel_b)?; phys.add_vel_bnd(bnd_t, vel_t)?;
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-10, 1.0)?;
    let flux_l = end_flux(&vars, bnd_l, vel, 1.0);
    let flux_r = end_flux(&vars, bnd_r, vel, 1.0);
    Ok((positions(&vars, dom), vars.vecdom[vel].value.clone(), flux_l, flux_r))
}

// (rms error, max error) normalised by max |u_exact|; v error is included in max
fn errors(pos: &[(f64, f64)], vel: &[(f64, f64)], exact: &dyn Fn(f64) -> f64) -> (f64, f64) {
    let umax = pos.iter().map(|p| exact(p.1).abs()).fold(0.0, f64::max);
    let mut sum = 0.0;
    let mut emax: f64 = 0.0;
    for (p, v) in pos.iter().zip(vel) {
        let e = v.0 - exact(p.1);
        sum += e * e;
        emax = emax.max(e.abs()).max(v.1.abs());
    }
    ((sum / pos.len() as f64).sqrt() / umax, emax / umax)
}

fn order(name: &str, dxs: &[f64], errs: &[(f64, f64)]) {
    for i in 0..dxs.len() {
        let p = if i > 0 { (errs[i - 1].0 / errs[i].0).ln() / (dxs[i - 1] / dxs[i]).ln() } else { f64::NAN };
        println!("RESULT {name}: dx = {:.4}; rms = {:.3e}; max = {:.3e}; order(rms) = {:.2}", dxs[i], errs[i].0, errs[i].1, p);
    }
}

// 1D reference for the channel: d/dy(visc * du/dy) - frac * (mu / K + rho * F / sqrt(K) * |u|) * u + G = 0 with mu = rho = 1
// porous layer y < a with uniform K and F; visc = mu_eff is visc_p in the porous layer and visc_f above
// u and visc * du/dy are continuous at y = a
// node-based conservative FD, Picard on |u|, Thomas algorithm
fn layered_1d(n: usize, g: f64, a: f64, visc_p: f64, visc_f: f64, kinv: f64, fsk: f64) -> Vec<f64> {
    let h = 1.0 / n as f64;
    let visc = |y: f64| if y < a { visc_p } else { visc_f };
    let mut u: Vec<f64> = vec![0.0; n + 1];
    for _ in 0..500 {
        let (mut lo, mut di, mut up, mut rhs) = (vec![0.0; n + 1], vec![1.0; n + 1], vec![0.0; n + 1], vec![0.0; n + 1]);
        for i in 1..n {
            let y = i as f64 * h;
            // fraction of the node cell [y - h/2, y + h/2] inside the porous layer
            let frac = ((a - (y - 0.5 * h)) / h).clamp(0.0, 1.0);
            let drag = frac * (kinv + fsk * u[i].abs());
            let (visc_lo, visc_up) = (visc(y - 0.5 * h), visc(y + 0.5 * h));
            lo[i] = -visc_lo / (h * h);
            up[i] = -visc_up / (h * h);
            di[i] = (visc_lo + visc_up) / (h * h) + drag;
            rhs[i] = g;
        }
        for i in 1..=n {
            let m = lo[i] / di[i - 1];
            di[i] -= m * up[i - 1];
            rhs[i] -= m * rhs[i - 1];
        }
        let mut new = vec![0.0; n + 1];
        new[n] = rhs[n] / di[n];
        for i in (0..n).rev() { new[i] = (rhs[i] - up[i] * new[i + 1]) / di[i]; }
        let diff = new.iter().zip(&u).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
        u = new;
        if diff < 1e-14 { break; }
    }
    u
}

fn interp(u: &[f64], y: f64) -> f64 {
    let n = u.len() - 1;
    let t = (y * n as f64).clamp(0.0, n as f64 - 1e-12);
    let i = (t.floor() as usize).min(n - 1);
    u[i] + (t - i as f64) * (u[i + 1] - u[i])
}

// int u dy of the 1D reference (trapezoid)
fn integral(u: &[f64]) -> f64 {
    let h = 1.0 / (u.len() - 1) as f64;
    u.windows(2).map(|w| 0.5 * (w[0] + w[1]) * h).sum()
}

#[test]
fn brinkman_channel() {
    // uniform epsl = 0.5 and diam = sqrt(3), so that K = 0.01 and F / sqrt(K) = 1.75 / (0.5^2 * sqrt(3)); G = 1
    // mu_eff = mu / epsl; checks the boundary half-edge viscosity of PoreOpVecBndTranspose and PoreOpVecBndDilatatory at the open ends
    let (kinv, fsk) = ergun(0.5, 3f64.sqrt());
    let ref_1d = layered_1d(20000, 1.0, 1.0, 2.0, 2.0, kinv, fsk);
    let flux_1d = integral(&ref_1d);
    let exact = move |y: f64| interp(&ref_1d, y);
    let dxs = [0.1, 0.05, 0.025, 0.0125];
    let mut errs = vec![];
    let mut flux_errs = vec![];
    for &dx in &dxs {
        let (pos, vel, flux_l, flux_r) = channel(dx, 1.0, |_| 0.5, |_| 3f64.sqrt()).unwrap();
        errs.push(errors(&pos, &vel, &exact));
        let flux_err = ((flux_l - flux_1d).abs().max((flux_r - flux_1d).abs())) / flux_1d;
        println!("RESULT brinkman_channel: dx = {dx:.4}; flux_l = {flux_l:.6e}; flux_r = {flux_r:.6e}; flux_1d = {flux_1d:.6e}; max rel err = {flux_err:.3e}");
        flux_errs.push(flux_err);
    }
    order("brinkman_channel", &dxs, &errs);
    // the error keeps decreasing under refinement and the open ends carry the inlet mass flux
    assert!(errs[3].0 < 1e-3 && errs[3].1 < 0.25 * errs[1].1);
    assert!(flux_errs[3] < 0.5 * flux_errs[1] && flux_errs[3] < 1e-2);
}

#[test]
fn vafai_kim_layered() {
    // Vafai-Kim (1990) configuration: porous layer (epsl = 0.5, diam = sqrt(3), so K = 0.01) on the wall y < 1/2, clear fluid (epsl = 1) above
    // mu_eff = mu / epsl jumps from 2 in the porous layer to 1 in the clear fluid
    // driven by G = 1 (Brinkman-dominated) and G = 100 (Forchheimer drag significant)
    // the interface is not aligned with the random mesh, and epsl is evaluated at the vertices
    let epsl: fn(f64) -> f64 = |y| if y < 0.5 { 0.5 } else { 1.0 };
    let diam: fn(f64) -> f64 = |_| 3f64.sqrt();
    let (kinv, fsk) = ergun(0.5, 3f64.sqrt());
    let a = 0.5f64;
    let visc_p = 2.0;

    // check the 1D reference against the analytic Darcy-Brinkman solution (F = 0); s = 1 / sqrt(K * visc_p)
    // porous: u = K (1 - cosh(s y)) + B sinh(s y); fluid: u = -(y^2 - 1) / 2 + C (y - 1)
    let k = 1.0 / kinv;
    let s = (kinv / visc_p).sqrt();
    // u(a): K (1 - ch) + B sh = -(a^2 - 1)/2 + C (a - 1); visc * u'(a): visc_p (-K s sh + B s ch) = -a + C
    let (sh, ch) = ((s * a).sinh(), (s * a).cosh());
    let (a11, a12, b1) = (sh, -(a - 1.0), -(a * a - 1.0) / 2.0 - k * (1.0 - ch));
    let (a21, a22, b2) = (visc_p * s * ch, -1.0, -a + visc_p * k * s * sh);
    let det = a11 * a22 - a12 * a21;
    let bb = (b1 * a22 - a12 * b2) / det;
    let cc = (a11 * b2 - b1 * a21) / det;
    let analytic = move |y: f64| if y < a { k * (1.0 - (s * y).cosh()) + bb * (s * y).sinh() } else { -(y * y - 1.0) / 2.0 + cc * (y - 1.0) };
    let ref_1d = layered_1d(20000, 1.0, a, visc_p, 1.0, kinv, 0.0);
    let check = (0..=100).map(|i| (interp(&ref_1d, i as f64 / 100.0) - analytic(i as f64 / 100.0)).abs()).fold(0.0, f64::max);
    println!("RESULT vafai_kim 1D reference vs analytic (F = 0): max = {:.3e}", check / analytic(0.75));
    assert!(check / analytic(0.75) < 1e-6);

    let dxs = [0.05, 0.025, 0.0125];
    for g in [1.0, 100.0] {
        let ref_1d = layered_1d(20000, g, a, visc_p, 1.0, kinv, fsk);
        let exact = move |y: f64| interp(&ref_1d, y);
        let mut errs = vec![];
        for &dx in &dxs {
            let (pos, vel, _, _) = channel(dx, g, epsl, diam).unwrap();
            errs.push(errors(&pos, &vel, &exact));
        }
        order(&format!("vafai_kim G={g}"), &dxs, &errs);
        assert!(errs[2].0 < 2e-2 && errs[2].0 < errs[0].0);
    }
}

// ---------------------------------------------------------------- packed column: Ergun pressure drop

// water through a bed of diam = 1 mm, epsl = 0.4; uniform inlet, top and bottom held at the inlet velocity (no wall layer)
// returns (dp/dx from a least-squares fit, max |vel - U| / U)
fn column(u_in: f64) -> Result<(f64, f64), VFVChemError> {
    let Mesh { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t] } = mesh(0.1, 0.02, 0.002)?;
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 0.001, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    let epsl = vars.add_scldom_cnst(dom, 0.4, "")?;
    let diam = vars.add_scldom_cnst(dom, 1e-3, "")?;
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (u_in, 0.0), "")?;
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (u_in, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (u_in, 0.0), "")?;
    let mut phys = SteadyPore::new()?;
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, diam)?;
    phys.add_vel_bnd(bnd_l, vel_l)?; phys.add_pres_bnd(bnd_r, pres_r)?;
    phys.add_vel_bnd(bnd_b, vel_b)?; phys.add_vel_bnd(bnd_t, vel_t)?;
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-10, 1.0)?;
    let pos = positions(&vars, dom);
    let p = &vars.scldom[pres].value;
    let n = vars.dom[dom].num_pts();
    let (mx, mp) = ((0..n).map(|i| pos[i].0).sum::<f64>() / n as f64, (0..n).map(|i| p[i]).sum::<f64>() / n as f64);
    let slope = (0..n).map(|i| (pos[i].0 - mx) * (p[i] - mp)).sum::<f64>() / (0..n).map(|i| (pos[i].0 - mx).powi(2)).sum::<f64>();
    let dev = vars.vecdom[vel].value.iter().map(|v| (v.0 - u_in).abs().max(v.1.abs())).fold(0.0, f64::max) / u_in;
    Ok((slope, dev))
}

#[test]
fn darcy_ergun_column() {
    let (mu, rho, diam): (f64, f64, f64) = (1e-3, 1000.0, 1e-3);
    let (kinv, fsk) = ergun(0.4, diam);
    for u_in in [1e-4, 1e-2, 5e-2] {
        let ergun = mu * u_in * kinv + rho * fsk * u_in * u_in;
        let darcy = mu * u_in * kinv;
        let (slope, dev) = column(u_in).unwrap();
        let rel = (-slope - ergun) / ergun;
        println!("RESULT column U = {u_in:.0e}: Re_p = {:.2}; -dp/dx = {:.6e}; Ergun = {:.6e}; Darcy-only = {:.6e}; rel err = {:.3e}; max|vel - U|/U = {:.3e}", rho * u_in * diam / mu, -slope, ergun, darcy, rel, dev);
        assert!(rel.abs() < 1e-6 && dev < 1e-6);
    }
}

// ---------------------------------------------------------------- transient Darcy-Forchheimer spin-up

#[test]
fn transient_spinup() {
    // epsl = 0.5 and diam = sqrt(0.3), so that K = 1e-3; mu = rho = 1, G = 1
    // channel core (Brinkman layer sqrt(K mu_eff / mu) = 0.045 << 1/2): (rho / epsl) du/dt = G - (mu / K + rho F / sqrt(K) |u|) u
    // backward Euler: (F / sqrt(K)) dt u_n^2 + (1 / epsl + dt / K) u_n - (u_{n-1} / epsl + dt G) = 0; tau = rho K / (epsl mu) = 2e-3
    let (kinv, fsk) = ergun(0.5, 0.3f64.sqrt());
    let Mesh { mut vars, dom, bnd: [bnd_l, bnd_r, bnd_b, bnd_t] } = mesh(0.2, 1.0, 0.05).unwrap();
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "").unwrap();
    let pres = vars.add_scldom_unkn(dom, 0.0, "").unwrap();
    let den = vars.add_scldom_cnst(dom, 1.0, "").unwrap();
    let visc = vars.add_scldom_cnst(dom, 1.0, "").unwrap();
    let fce = vars.add_vecdom_cnst(dom, (1.0, 0.0), "").unwrap();
    let epsl = vars.add_scldom_cnst(dom, 0.5, "").unwrap();
    let diam = vars.add_scldom_cnst(dom, 0.3f64.sqrt(), "").unwrap();
    let pres_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "").unwrap();
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "").unwrap();
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "").unwrap();
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "").unwrap();
    let mut phys = TransientPore::new().unwrap();
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, diam).unwrap();
    phys.add_pres_bnd(bnd_l, pres_l).unwrap(); phys.add_pres_bnd(bnd_r, pres_r).unwrap();
    phys.add_vel_bnd(bnd_b, vel_b).unwrap(); phys.add_vel_bnd(bnd_t, vel_t).unwrap();
    let mut lin = LinearLu::new(1).unwrap();
    let dt = 5e-4;
    let mut u_be = 0.0;
    let pos = positions(&vars, dom);
    for step in 1..=3 {
        // one step per call; previous values are carried in vars
        phys.solve(&mut vars, &mut lin, dt, 1, 100, 50, 1e-12, 1.0).unwrap();
        let (qa, qb, qc) = (fsk * dt, 2.0 + kinv * dt, 2.0 * u_be + dt);
        u_be = (-qb + (qb * qb + 4.0 * qa * qc).sqrt()) / (2.0 * qa);
        let core: Vec<f64> = pos.iter().zip(&vars.vecdom[vel].value).filter(|(p, _)| (p.1 - 0.5).abs() < 0.2 && (p.0 - 0.1).abs() < 0.05).map(|(_, v)| v.0).collect();
        let err = core.iter().map(|u| (u - u_be).abs()).fold(0.0, f64::max) / u_be;
        println!("RESULT transient_spinup step {step}: u_BE = {:.6e}; core max rel err = {:.3e}", u_be, err);
        assert!(err < 2e-3);
    }
}
