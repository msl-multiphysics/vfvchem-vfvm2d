use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Steady-state local thermal non-equilibrium heat and momentum transfer in a catalytic packed bed with an exothermic surface reaction.
/// Run with: `cargo run --release --example poreheat_packedbed`
///
/// Geometry:
/// Packed bed of 5 mm catalyst pellets; air flows from the left inlet to the right outlet; the bottom wall is cooled.
/// Prints the peak solid temperature, the largest solid-fluid temperature difference, and the global heat balance
/// (heat released, net advective outflow, and heat removed through the inlet and the cooled wall).
///
///                  bnd_t
///       +-----------------------+
/// bnd_l |          dom          | bnd_r
///       +-----------------------+
///                  bnd_b
///
/// Domains:
/// dom - packed bed (0.2 m x 0.05 m)
///
/// Boundaries:
/// bnd_l - inlet (0.05 m)
/// bnd_r - outlet (0.05 m)
/// bnd_b - cooled bottom wall (0.2 m)
/// bnd_t - top wall (0.2 m)
///
/// Heat Transfer
/// T_f is the fluid temperature and T_s is the solid temperature.
///
/// Properties:
/// [non-constant properties in SI units]
/// dom - fluid volumetric heat capacity ((rho cp)_f = 1.0 * 1010.0 J m-3 K-1)
/// dom - fluid conductivity (k_f = 0.03 W m-1 K-1)
/// dom - solid volumetric heat capacity ((rho c)_s = 2.3e6 J m-3 K-1); not used in steady state
/// dom - solid conductivity (k_s = 1.0 W m-1 K-1); alumina pellets
/// dom - heat transfer coefficient (h = 110.0 W m-2 K-1); Wakao-Kaguei with Re_p = 109
/// dom - specific surface area (a_v = 6 * (1 - epsl) / diam = 720.0 m-1)
/// dom - longitudinal dispersivity (alpha_L = diam / 2 = 2.5e-3 m)
/// dom - transverse dispersivity (alpha_T = diam / 10 = 5.0e-4 m)
/// dom - volumetric production coefficients and heat sources are zero
/// dom - surface heat release (q_a = 40.0 * exp(-3000 * (1 / T_s - 1 / 400)) W m-2); Arrhenius
///
/// Boundary conditions:
/// bnd_l - fluid temperature (T_f = 400.0 K); solid no flux (q_s = 0.0 W m-2)
/// bnd_r - fluid outflow; solid no flux (q_s = 0.0 W m-2)
/// bnd_b - fluid and solid temperature (T_f = T_s = 350.0 K)
/// bnd_t - fluid and solid no flux (q_f = q_s = 0.0 W m-2)
///
/// Momentum Transfer
///
/// Properties:
/// dom - density (rho = 1.0 kg m-3)
/// dom - viscosity (mu = 2.3e-5 Pa s)
/// dom - body force (f = <0.0, 0.0> N m-3)
/// dom - porosity (epsl = 0.4)
/// dom - permeability (K = 2.96e-8 m2); Ergun with particle diameter diam = 5.0e-3 m
/// dom - Forchheimer coefficient (F = 0.565); Ergun
///
/// Boundary conditions:
/// bnd_l - superficial velocity (v = <0.5, 0.0> m s-1)
/// bnd_r - pressure (p = 0.0 Pa)
/// bnd_b, bnd_t - free-slip walls (v = <0.5, 0.0> m s-1); plug flow
///
fn main() -> Result<(), VFVChemError> {
    create_dir_all("examples/output_poreheat_packedbed").unwrap();
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.2, 0.05, 2.5e-3, 2.5e-3, 2.5e-3, 100, 42)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let templ = vars.add_scldom_unkn(dom, 400.0, "examples/output_poreheat_packedbed/templ.vtu")?;
    let temps = vars.add_scldom_unkn(dom, 400.0, "examples/output_poreheat_packedbed/temps.vtu")?;
    let vel = vars.add_vecdom_unkn(dom, (0.5, 0.0), "examples/output_poreheat_packedbed/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_poreheat_packedbed/pres.vtu")?;

    // flow properties
    let den = vars.add_scldom_cnst(dom, 1.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 2.3e-5, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    let epsl = vars.add_scldom_cnst(dom, 0.4, "")?;
    let perm = vars.add_scldom_cnst(dom, 2.96e-8, "")?;  // permeability K
    let forch = vars.add_scldom_cnst(dom, 0.565, "")?;  // Forchheimer coefficient F

    // heat properties
    let vlcpl = vars.add_scldom_cnst(dom, 1010.0, "")?;  // (rho cp)_f
    let condl = vars.add_scldom_cnst(dom, 0.03, "")?;  // k_f
    let vlcps = vars.add_scldom_cnst(dom, 2.3e6, "")?;  // (rho c)_s; not used in steady state
    let conds = vars.add_scldom_cnst(dom, 1.0, "")?;  // k_s
    let htc = vars.add_scldom_cnst(dom, 110.0, "")?;  // h
    let area = vars.add_scldom_cnst(dom, 720.0, "")?;  // a_v
    let displ = vars.add_scldom_cnst(dom, 2.5e-3, "")?;  // alpha_L
    let dispt = vars.add_scldom_cnst(dom, 5e-4, "")?;  // alpha_T
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let asrc = vars.add_scldom_func(dom, Box::new(|_t, (_x, _y), unkn| 40.0 * (-3000.0 * (1.0 / unkn[0] - 1.0 / 400.0)).exp()), vec![temps], "")?;

    // boundary values
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.5, 0.0), "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.5, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.5, 0.0), "")?;
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let templ_l = vars.add_sclbnd_cnst(bnd_l, 400.0, "")?;
    let templ_b = vars.add_sclbnd_cnst(bnd_b, 350.0, "")?;
    let temps_b = vars.add_sclbnd_cnst(bnd_b, 350.0, "")?;
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;
    let hflx_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;

    // steady heat transfer in porous media
    let mut phys = SteadyHeatPore::new()?;
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, forch)?;  // arguments: domain, v, p, rho, mu, f, epsl, K, F
    phys.add_heatl_dom(dom, templ, vlcpl, condl, zero, zero)?;  // arguments: domain, T_f, (rho cp)_f, k_f, s_f, q_f
    phys.add_heats_dom(dom, temps, vlcps, conds, zero, zero, zero, asrc)?;  // arguments: domain, T_s, (rho c)_s, k_s, s_s, q_s, s_a, q_a
    phys.add_matrix_dom(dom, htc, area, displ, dispt)?;  // arguments: domain, h, a_v, alpha_L, alpha_T
    phys.add_templ_bnd(bnd_l, templ_l)?;
    phys.add_houtl_bnd(bnd_r)?;
    phys.add_templ_bnd(bnd_b, templ_b)?;
    phys.add_hflxl_bnd(bnd_t, hflx_t)?;
    phys.add_hflxs_bnd(bnd_l, hflx_l)?;
    phys.add_hflxs_bnd(bnd_r, hflx_r)?;
    phys.add_temps_bnd(bnd_b, temps_b)?;
    phys.add_hflxs_bnd(bnd_t, hflx_t)?;
    phys.add_vel_bnd(bnd_l, vel_l)?;
    phys.add_pres_bnd(bnd_r, pres_r)?;
    phys.add_vel_bnd(bnd_b, vel_b)?;
    phys.add_vel_bnd(bnd_t, vel_t)?;

    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-8, 1.0)?;

    Ok(())
}
