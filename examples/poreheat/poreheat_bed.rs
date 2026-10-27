use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Steady-state local thermal non-equilibrium heat and momentum transfer in a channel with a catalytic packed bed under a clear fluid layer and an exothermic surface reaction.
/// Run with: `cargo run --release --example poreheat_bed`
///
/// Geometry:
/// Channel with a packed bed in the lower half; flow from the left inlets to the right outlets.
/// Most of the flow bypasses the bed through the clear fluid layer, which carries away the heat released in the bed.
/// Prints the peak solid temperature, the fraction of the flow that bypasses the bed, and the global heat balance
/// (heat released, net advective outflow, and heat conducted out through the inlets).
///
///                    bnd_t
///        +-----------------------+
/// bnd_lf |       dom_fluid       | bnd_rf
///        +------- itf_bed -------+
/// bnd_lb |        dom_bed        | bnd_rb
///        +-----------------------+
///                    bnd_b
///
/// Domains:
/// dom_bed - packed bed (0.1 m x 0.01 m)
/// dom_fluid - clear fluid (0.1 m x 0.01 m)
///
/// Boundaries:
/// bnd_lb - inlet of dom_bed (0.01 m)
/// bnd_lf - inlet of dom_fluid (0.01 m)
/// bnd_rb - outlet of dom_bed (0.01 m)
/// bnd_rf - outlet of dom_fluid (0.01 m)
/// bnd_b - bottom wall of dom_bed (0.1 m)
/// bnd_t - top wall of dom_fluid (0.1 m)
///
/// Interfaces:
/// itf_bed - boundary between dom_bed and dom_fluid (0.1 m)
///
/// Heat Transfer
/// T_f is the fluid temperature and T_s is the solid temperature.
/// The clear fluid is a porous medium with a_v = 0 and alpha_L = alpha_T = 0, and without a solid.
///
/// Properties:
/// [non-constant properties in SI units]
/// dom_bed - fluid volumetric heat capacity ((rho cp)_f = 1000.0 * 4180.0 J m-3 K-1); water
/// dom_bed - fluid conductivity (k_f = 0.6 W m-1 K-1)
/// dom_bed - solid volumetric heat capacity ((rho c)_s = 2.3e6 J m-3 K-1); not used in steady state
/// dom_bed - solid conductivity (k_s = 1.0 W m-1 K-1); alumina pellets
/// dom_bed - heat transfer coefficient (h = 1560.0 W m-2 K-1); Wakao-Kaguei with Re_p = 2 and Pr = 7
/// dom_bed - specific surface area (a_v = 6 * (1 - epsl) / diam = 1800.0 m-1)
/// dom_bed - longitudinal dispersivity (alpha_L = diam / 2 = 1.0e-3 m)
/// dom_bed - transverse dispersivity (alpha_T = diam / 10 = 2.0e-4 m)
/// dom_bed - volumetric production coefficients and heat sources are zero
/// dom_bed - surface heat release (q_a = 10.0 * exp(-3000 * (1 / T_s - 1 / 300)) W m-2); Arrhenius
/// dom_fluid - fluid volumetric heat capacity ((rho cp)_f = 1000.0 * 4180.0 J m-3 K-1)
/// dom_fluid - fluid conductivity (k_f = 0.6 W m-1 K-1)
/// dom_fluid - heat transfer coefficient (h = 0.0 W m-2 K-1)
/// dom_fluid - specific surface area (a_v = 0.0 m-1)
/// dom_fluid - longitudinal dispersivity (alpha_L = 0.0 m)
/// dom_fluid - transverse dispersivity (alpha_T = 0.0 m)
/// dom_fluid - production coefficient and heat source are zero
///
/// Boundary conditions:
/// bnd_lb - fluid temperature (T_f = 300.0 K); solid no flux (q_s = 0.0 W m-2)
/// bnd_lf - fluid temperature (T_f = 300.0 K)
/// bnd_rb - fluid outflow; solid no flux (q_s = 0.0 W m-2)
/// bnd_rf - fluid outflow
/// bnd_b - fluid and solid no flux (q_f = q_s = 0.0 W m-2)
/// bnd_t - fluid no flux (q_f = 0.0 W m-2)
///
/// Interface conditions:
/// itf_bed - fluid temperature and heat flux continuity; solid no flux (adiabatic solid surface)
///
/// Momentum Transfer
/// The clear fluid is a porous medium with epsl = 1, K = inf, and F = 0.
///
/// Properties:
/// dom_bed - density (rho = 1000.0 kg m-3)
/// dom_bed - viscosity (mu = 0.001 Pa s)
/// dom_bed - body force (f = <0.0, 0.0> N m-3)
/// dom_bed - porosity (epsl = 0.4)
/// dom_bed - permeability (K = 4.74e-9 m2); Ergun with particle diameter diam = 2.0e-3 m
/// dom_bed - Forchheimer coefficient (F = 0.565); Ergun
/// dom_fluid - density (rho = 1000.0 kg m-3)
/// dom_fluid - viscosity (mu = 0.001 Pa s)
/// dom_fluid - body force (f = <0.0, 0.0> N m-3)
/// dom_fluid - porosity (epsl = 1.0)
/// dom_fluid - permeability (K = inf)
/// dom_fluid - Forchheimer coefficient (F = 0.0)
///
/// Boundary conditions:
/// bnd_lb, bnd_lf - superficial velocity (v = <1.0e-3, 0.0> m s-1)
/// bnd_rb, bnd_rf - pressure (p = 0.0 Pa)
/// bnd_b, bnd_t - no-slip (v = <0.0, 0.0> m s-1)
///
/// Interface conditions:
/// itf_bed - velocity and pressure continuity
///
fn main() -> Result<(), VFVChemError> {
    create_dir_all("examples/output_poreheat_bed").unwrap();
    let mut vars = Variables::new_obj("examples/obj/obj_bed.obj", 5e-4, 5e-4, 5e-4, 100, 42)?;
    let dom_b = vars.add_dom("dom_bed")?;
    let dom_f = vars.add_dom("dom_fluid")?;
    let bnd_lb = vars.add_bnd(dom_b, "bnd_lb")?;
    let bnd_rb = vars.add_bnd(dom_b, "bnd_rb")?;
    let bnd_b = vars.add_bnd(dom_b, "bnd_b")?;
    let bnd_lf = vars.add_bnd(dom_f, "bnd_lf")?;
    let bnd_rf = vars.add_bnd(dom_f, "bnd_rf")?;
    let bnd_t = vars.add_bnd(dom_f, "bnd_t")?;
    let itf_bed = vars.add_itf(dom_b, dom_f, "itf_bed")?;

    // unknowns
    let templ_b = vars.add_scldom_unkn(dom_b, 300.0, "examples/output_poreheat_bed/templ_bed.vtu")?;
    let temps_b = vars.add_scldom_unkn(dom_b, 300.0, "examples/output_poreheat_bed/temps_bed.vtu")?;
    let vel_b = vars.add_vecdom_unkn(dom_b, (1e-3, 0.0), "examples/output_poreheat_bed/vel_bed.vtu")?;
    let pres_b = vars.add_scldom_unkn(dom_b, 0.0, "examples/output_poreheat_bed/pres_bed.vtu")?;
    let templ_f = vars.add_scldom_unkn(dom_f, 300.0, "examples/output_poreheat_bed/templ_fluid.vtu")?;
    let vel_f = vars.add_vecdom_unkn(dom_f, (1e-3, 0.0), "examples/output_poreheat_bed/vel_fluid.vtu")?;
    let pres_f = vars.add_scldom_unkn(dom_f, 0.0, "examples/output_poreheat_bed/pres_fluid.vtu")?;

    // bed properties
    let den_b = vars.add_scldom_cnst(dom_b, 1000.0, "")?;  // density
    let visc_b = vars.add_scldom_cnst(dom_b, 1e-3, "")?;  // viscosity
    let fce_b = vars.add_vecdom_cnst(dom_b, (0.0, 0.0), "")?;  // body force
    let epsl_b = vars.add_scldom_cnst(dom_b, 0.4, "")?;  // porosity
    let perm_b = vars.add_scldom_cnst(dom_b, 4.74e-9, "")?;  // permeability K
    let forch_b = vars.add_scldom_cnst(dom_b, 0.565, "")?;  // Forchheimer coefficient F
    let htc_b = vars.add_scldom_cnst(dom_b, 1560.0, "")?;  // heat transfer coefficient h
    let area_b = vars.add_scldom_cnst(dom_b, 1800.0, "")?;  // specific surface area a_v
    let displ_b = vars.add_scldom_cnst(dom_b, 1e-3, "")?;  // longitudinal dispersivity
    let dispt_b = vars.add_scldom_cnst(dom_b, 2e-4, "")?;  // transverse dispersivity

    // fluid properties; a clear fluid is a pore domain with epsl = 1, K = inf, F = 0
    let den_f = vars.add_scldom_cnst(dom_f, 1000.0, "")?;
    let visc_f = vars.add_scldom_cnst(dom_f, 1e-3, "")?;
    let fce_f = vars.add_vecdom_cnst(dom_f, (0.0, 0.0), "")?;
    let epsl_f = vars.add_scldom_cnst(dom_f, 1.0, "")?;
    let perm_f = vars.add_scldom_cnst(dom_f, f64::INFINITY, "")?;
    let forch_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;
    let htc_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;
    let area_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;
    let displ_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;
    let dispt_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;

    // heat transfer properties; the surface heat release is multiplied by a_v
    let vlcpl_b = vars.add_scldom_cnst(dom_b, 1000.0 * 4180.0, "")?;  // (rho cp)_f
    let condl_b = vars.add_scldom_cnst(dom_b, 0.6, "")?;  // k_f
    let vlcps_b = vars.add_scldom_cnst(dom_b, 2.3e6, "")?;  // (rho c)_s; not used in steady state
    let conds_b = vars.add_scldom_cnst(dom_b, 1.0, "")?;  // k_s
    let zero_b = vars.add_scldom_cnst(dom_b, 0.0, "")?;
    let asrcs_b = vars.add_scldom_func(dom_b, Box::new(|_t, (_x, _y), unkn| 10.0 * (-3000.0 * (1.0 / unkn[0] - 1.0 / 300.0)).exp()), vec![temps_b], "")?;  // q_a
    let vlcpl_f = vars.add_scldom_cnst(dom_f, 1000.0 * 4180.0, "")?;
    let condl_f = vars.add_scldom_cnst(dom_f, 0.6, "")?;
    let zero_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;

    // boundary values
    let vel_lb = vars.add_vecbnd_cnst(bnd_lb, (1e-3, 0.0), "")?;
    let vel_lf = vars.add_vecbnd_cnst(bnd_lf, (1e-3, 0.0), "")?;
    let vel_w_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_w_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;
    let pres_rb = vars.add_sclbnd_cnst(bnd_rb, 0.0, "")?;
    let pres_rf = vars.add_sclbnd_cnst(bnd_rf, 0.0, "")?;
    let templ_lb = vars.add_sclbnd_cnst(bnd_lb, 300.0, "")?;
    let templ_lf = vars.add_sclbnd_cnst(bnd_lf, 300.0, "")?;
    let hflx_lb = vars.add_sclbnd_cnst(bnd_lb, 0.0, "")?;
    let hflx_rb = vars.add_sclbnd_cnst(bnd_rb, 0.0, "")?;
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;

    // steady heat transfer in porous media; the clear fluid has no solid heat domain
    let mut phys = SteadyHeatPore::new()?;
    phys.add_heatl_dom(dom_b, templ_b, vlcpl_b, condl_b, zero_b, zero_b)?;  // arguments: domain, T_f, (rho cp)_f, k_f, s_f, q_f
    phys.add_heats_dom(dom_b, temps_b, vlcps_b, conds_b, zero_b, zero_b, zero_b, asrcs_b)?;  // arguments: domain, T_s, (rho c)_s, k_s, s_s, q_s, s_a, q_a
    phys.add_heatl_dom(dom_f, templ_f, vlcpl_f, condl_f, zero_f, zero_f)?;
    phys.add_matrix_dom(dom_b, htc_b, area_b, displ_b, dispt_b)?;  // arguments: domain, h, a_v, alpha_L, alpha_T
    phys.add_matrix_dom(dom_f, htc_f, area_f, displ_f, dispt_f)?;
    phys.add_pore_dom(dom_b, vel_b, pres_b, den_b, visc_b, fce_b, epsl_b, perm_b, forch_b)?;  // arguments: domain, v, p, rho, mu, f, epsl, K, F
    phys.add_pore_dom(dom_f, vel_f, pres_f, den_f, visc_f, fce_f, epsl_f, perm_f, forch_f)?;
    phys.add_templ_bnd(bnd_lb, templ_lb)?;
    phys.add_templ_bnd(bnd_lf, templ_lf)?;
    phys.add_houtl_bnd(bnd_rb)?;
    phys.add_houtl_bnd(bnd_rf)?;
    phys.add_hflxl_bnd(bnd_b, hflx_b)?;
    phys.add_hflxl_bnd(bnd_t, hflx_t)?;
    phys.add_hflxs_bnd(bnd_lb, hflx_lb)?;
    phys.add_hflxs_bnd(bnd_rb, hflx_rb)?;
    phys.add_hflxs_bnd(bnd_b, hflx_b)?;
    phys.add_hcntl_itf(itf_bed)?;  // fluid heat continuity across the interface; without add_hcnts_itf the solid surface is adiabatic
    phys.add_vel_bnd(bnd_lb, vel_lb)?;
    phys.add_vel_bnd(bnd_lf, vel_lf)?;
    phys.add_pres_bnd(bnd_rb, pres_rb)?;
    phys.add_pres_bnd(bnd_rf, pres_rf)?;
    phys.add_vel_bnd(bnd_b, vel_w_b)?;
    phys.add_vel_bnd(bnd_t, vel_w_t)?;
    phys.add_fcnt_itf(itf_bed)?;  // flow continuity across the interface

    let mut lin = FaerLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-8, 1.0)?;

    Ok(())
}
