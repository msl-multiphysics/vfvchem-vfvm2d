use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Transient local thermal non-equilibrium heat and momentum transfer in a metal foam heated by a localized surface heat source.
/// Run with: `cargo run --release --example poreheat_surface`
///
/// Geometry:
/// Aluminum foam with a denser layer near the bottom; air flows from the left inlet to the right outlet.
/// The struts are heated over a Gaussian band (induction heating) that switches on at t = 0.
/// Prints the peak temperatures and the mean outlet fluid temperature at t = 15 s.
///
///                  bnd_t
///       +-----------------------+
/// bnd_l |          dom          | bnd_r
///       +-----------------------+
///                  bnd_b
///
/// Domains:
/// dom - metal foam (0.1 m x 0.02 m)
///
/// Boundaries:
/// bnd_l - inlet (0.02 m)
/// bnd_r - outlet (0.02 m)
/// bnd_b - base plate (0.1 m)
/// bnd_t - top wall (0.1 m)
///
/// Heat Transfer
/// T_f is the fluid temperature and T_s is the solid temperature.
///
/// Properties:
/// [non-constant properties in SI units]
/// dom - fluid volumetric heat capacity ((rho cp)_f = 1.2 * 1005.0 J m-3 K-1)
/// dom - fluid conductivity (k_f = 0.026 W m-1 K-1)
/// dom - solid volumetric heat capacity ((rho c)_s = 2.43e6 J m-3 K-1); aluminum
/// dom - solid conductivity (k_s = 20.0 W m-1 K-1); effective strut conductivity
/// dom - heat transfer coefficient (h = 60.0 W m-2 K-1)
/// dom - specific surface area (a_v = 1500.0 m-1)
/// dom - longitudinal dispersivity (alpha_L = 1.0e-3 m)
/// dom - transverse dispersivity (alpha_T = 2.0e-4 m)
/// dom - fluid production coefficient and heat source are zero
/// dom - solid production coefficient (s_s = -500.0 W m-3 K-1); linearized loss to the surroundings at 300 K
/// dom - solid heat source (q_s = 500.0 * 300.0 W m-3)
/// dom - surface heat source (q_a = 2000.0 * exp(-((x - 0.03) / 0.01)^2) W m-2)
///
/// Initial conditions:
/// dom - fluid and solid temperature (T_f = T_s = 300.0 K)
///
/// Boundary conditions:
/// bnd_l - fluid temperature (T_f = 300.0 K); solid no flux (q_s = 0.0 W m-2)
/// bnd_r - fluid outflow; solid no flux (q_s = 0.0 W m-2)
/// bnd_b - fluid no flux (q_f = 0.0 W m-2); solid temperature (T_s = 300.0 K)
/// bnd_t - fluid and solid no flux (q_f = q_s = 0.0 W m-2)
///
/// Momentum Transfer
///
/// Properties:
/// dom - density (rho = 1.2 kg m-3)
/// dom - viscosity (mu = 1.8e-5 Pa s)
/// dom - body force (f = <0.0, 0.0> N m-3)
/// dom - porosity (epsl = 0.80 + 0.12 * y / 0.02)
/// dom - permeability (K = 1.0e-7 * (epsl^3 / (1 - epsl)^2) / (0.86^3 / 0.14^2) m2); Kozeny-Carman scaled to 1.0e-7 m2 at the mean porosity
/// dom - Forchheimer coefficient (F = 0.1)
///
/// Initial conditions:
/// dom - superficial velocity (v = <1.0, 0.0> m s-1)
/// dom - pressure (p = 0.0 Pa)
///
/// Boundary conditions:
/// bnd_l - superficial velocity (v = <1.0, 0.0> m s-1)
/// bnd_r - pressure (p = 0.0 Pa)
/// bnd_b, bnd_t - no-slip (v = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    create_dir_all("examples/output_poreheat_surface").unwrap();
    // 0.5 mm spacing resolves the Brinkman layer at the top wall, where K is largest; 1 mm does not converge
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.1, 0.02, 5e-4, 5e-4, 5e-4, 100, 42)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let templ = vars.add_scldom_unkn(dom, 300.0, "examples/output_poreheat_surface/templ.vtu")?;
    let temps = vars.add_scldom_unkn(dom, 300.0, "examples/output_poreheat_surface/temps.vtu")?;
    let vel = vars.add_vecdom_unkn(dom, (1.0, 0.0), "examples/output_poreheat_surface/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_poreheat_surface/pres.vtu")?;

    // flow properties; Kozeny-Carman permeability scaled to 1e-7 m2 at the mean porosity
    let den = vars.add_scldom_cnst(dom, 1.2, "")?;
    let visc = vars.add_scldom_cnst(dom, 1.8e-5, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    let epsl = vars.add_scldom_func(dom, Box::new(|_t, (_x, y), _u| 0.80 + 0.12 * y / 0.02), vec![], "")?;
    let perm = vars.add_scldom_func(dom, Box::new(|_t, (_x, y), _u| {
        let e: f64 = 0.80 + 0.12 * y / 0.02;
        1e-7 * (e.powi(3) / (1.0 - e).powi(2)) / (0.86f64.powi(3) / 0.14f64.powi(2))
    }), vec![], "")?;  // permeability K
    let forch = vars.add_scldom_cnst(dom, 0.1, "")?;  // Forchheimer coefficient F

    // heat properties
    let vlcpl = vars.add_scldom_cnst(dom, 1.2 * 1005.0, "")?;  // (rho cp)_f
    let condl = vars.add_scldom_cnst(dom, 0.026, "")?;  // k_f
    let vlcps = vars.add_scldom_cnst(dom, 2.43e6, "")?;  // (rho c)_s
    let conds = vars.add_scldom_cnst(dom, 20.0, "")?;  // k_s
    let htc = vars.add_scldom_cnst(dom, 60.0, "")?;  // h
    let area = vars.add_scldom_cnst(dom, 1500.0, "")?;  // a_v
    let displ = vars.add_scldom_cnst(dom, 1e-3, "")?;  // alpha_L
    let dispt = vars.add_scldom_cnst(dom, 2e-4, "")?;  // alpha_T
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let hprds = vars.add_scldom_cnst(dom, -500.0, "")?;  // s_s
    let hsrcs = vars.add_scldom_cnst(dom, 500.0 * 300.0, "")?;  // q_s
    let asrcs = vars.add_scldom_func(dom, Box::new(|_t, (x, _y), _u| 2000.0 * (-((x - 0.03) / 0.01f64).powi(2)).exp()), vec![], "")?;  // q_a

    // boundary values
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (1.0, 0.0), "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let templ_l = vars.add_sclbnd_cnst(bnd_l, 300.0, "")?;
    let temps_b = vars.add_sclbnd_cnst(bnd_b, 300.0, "")?;
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;
    let hflx_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;

    // transient heat transfer in porous media
    let mut phys = TransientHeatPore::new()?;
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, forch)?;  // arguments: domain, v, p, rho, mu, f, epsl, K, F
    phys.add_heatl_dom(dom, templ, vlcpl, condl, zero, zero)?;  // arguments: domain, T_f, (rho cp)_f, k_f, s_f, q_f
    phys.add_heats_dom(dom, temps, vlcps, conds, hprds, hsrcs, zero, asrcs)?;  // arguments: domain, T_s, (rho c)_s, k_s, s_s, q_s, s_a, q_a
    phys.add_matrix_dom(dom, htc, area, displ, dispt)?;  // arguments: domain, h, a_v, alpha_L, alpha_T
    phys.add_templ_bnd(bnd_l, templ_l)?;
    phys.add_houtl_bnd(bnd_r)?;
    phys.add_hflxl_bnd(bnd_b, hflx_b)?;
    phys.add_hflxl_bnd(bnd_t, hflx_t)?;
    phys.add_hflxs_bnd(bnd_l, hflx_l)?;
    phys.add_hflxs_bnd(bnd_r, hflx_r)?;
    phys.add_temps_bnd(bnd_b, temps_b)?;
    phys.add_hflxs_bnd(bnd_t, hflx_t)?;
    phys.add_vel_bnd(bnd_l, vel_l)?;
    phys.add_pres_bnd(bnd_r, pres_r)?;
    phys.add_vel_bnd(bnd_b, vel_b)?;
    phys.add_vel_bnd(bnd_t, vel_t)?;

    // 15 s; dt = 1 s
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 1.0, 15, 5, 50, 1e-8, 1.0)?;

    // post-processing; mean outlet T_f over the points within 2 mm of the outlet
    let pos: Vec<(f64, f64)> = vars.dom[dom].pts_gpt.iter().map(|&g| vars.pts.pts_pos[g]).collect();
    let tf = &vars.scldom[templ].value;
    let ts = &vars.scldom[temps].value;
    let ts_max = ts.iter().cloned().fold(f64::MIN, f64::max);
    let tf_max = tf.iter().cloned().fold(f64::MIN, f64::max);
    let out: Vec<f64> = pos.iter().zip(tf).filter(|(p, _)| p.0 > 0.098).map(|(_, t)| *t).collect();
    let tf_out = out.iter().sum::<f64>() / out.len() as f64;
    println!("t = 15 s; peak T_s = {ts_max:.2} K; peak T_f = {tf_max:.2} K; mean outlet T_f = {tf_out:.2} K");
    
    Ok(())
}
