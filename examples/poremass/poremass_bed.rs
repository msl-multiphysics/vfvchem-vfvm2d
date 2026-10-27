use vfvchem_vfvm2d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Steady-state mass and momentum transfer in a channel with a catalytic packed bed under a clear fluid layer and first-order surface reaction 0 -> 1.
/// Run with: `cargo run --release --example poremass_bed`
///
/// Geometry:
/// Channel with a packed bed in the lower half; flow from the left inlets to the right outlets.
/// Most of the flow bypasses the bed through the clear fluid layer.
/// Prints the fraction of the flow that bypasses the bed, the outlet conversion, and the global mass balance.
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
/// Mass Transfer
/// Concentrations are per unit fluid mass and diffusion coefficients include the density.
/// The clear fluid is a porous medium with tau = 1, a_v = 0, and alpha_L = alpha_T = 0.
///
/// Properties:
/// [non-constant properties in SI units]
/// dom_bed - diffusion coefficients (rho*D; kg m-1 s-1)
/// [D_00 D_01] = [1.0e-6 0.0]
/// [D_10 D_11]   [0.0    1.0e-6]
/// dom_bed - tortuosity (tau = 1.5)
/// dom_bed - specific surface area (a_v = 6 * (1 - epsl) / diam = 1800.0 m-1)
/// dom_bed - longitudinal dispersivity (alpha_L = diam / 2 = 1.0e-3 m)
/// dom_bed - transverse dispersivity (alpha_T = diam / 10 = 2.0e-4 m)
/// dom_bed - volumetric production coefficients and reaction rates are zero
/// dom_bed - surface production coefficients (kg m-2 s-1; positive if production)
/// [s_s,0] = [-rho * 5.0e-6]
/// [s_s,1]   [0.0]
/// dom_bed - surface reaction rates (mol m-2 s-1; positive if production)
/// [R_s,0] = [0.0]
/// [R_s,1]   [+rho * 5.0e-6 * c_0]
/// dom_fluid - diffusion coefficients (rho*D; kg m-1 s-1)
/// [D_00 D_01] = [1.0e-6 0.0]
/// [D_10 D_11]   [0.0    1.0e-6]
/// dom_fluid - tortuosity (tau = 1.0)
/// dom_fluid - specific surface area (a_v = 0.0 m-1)
/// dom_fluid - longitudinal dispersivity (alpha_L = 0.0 m)
/// dom_fluid - transverse dispersivity (alpha_T = 0.0 m)
/// dom_fluid - production coefficients and reaction rates are zero
///
/// Boundary conditions:
/// bnd_lb, bnd_lf - concentration (mol kg-1)
/// [c_0] = [1.0]
/// [c_1]   [0.0]
/// bnd_rb, bnd_rf - outflow
/// bnd_b, bnd_t - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
///
/// Interface conditions:
/// itf_bed - concentration and flux continuity
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
    create_dir_all("examples/output_poremass_bed").unwrap();
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
    let conc_0_b = vars.add_scldom_unkn(dom_b, 0.0, "examples/output_poremass_bed/conc_0_bed.vtu")?;
    let conc_1_b = vars.add_scldom_unkn(dom_b, 0.0, "examples/output_poremass_bed/conc_1_bed.vtu")?;
    let vel_b = vars.add_vecdom_unkn(dom_b, (1e-3, 0.0), "examples/output_poremass_bed/vel_bed.vtu")?;
    let pres_b = vars.add_scldom_unkn(dom_b, 0.0, "examples/output_poremass_bed/pres_bed.vtu")?;
    let conc_0_f = vars.add_scldom_unkn(dom_f, 0.0, "examples/output_poremass_bed/conc_0_fluid.vtu")?;
    let conc_1_f = vars.add_scldom_unkn(dom_f, 0.0, "examples/output_poremass_bed/conc_1_fluid.vtu")?;
    let vel_f = vars.add_vecdom_unkn(dom_f, (1e-3, 0.0), "examples/output_poremass_bed/vel_fluid.vtu")?;
    let pres_f = vars.add_scldom_unkn(dom_f, 0.0, "examples/output_poremass_bed/pres_fluid.vtu")?;

    // bed properties
    let den_b = vars.add_scldom_cnst(dom_b, 1000.0, "")?;  // density
    let visc_b = vars.add_scldom_cnst(dom_b, 1e-3, "")?;  // viscosity
    let fce_b = vars.add_vecdom_cnst(dom_b, (0.0, 0.0), "")?;  // body force
    let epsl_b = vars.add_scldom_cnst(dom_b, 0.4, "")?;  // porosity
    let perm_b = vars.add_scldom_cnst(dom_b, 4.74e-9, "")?;  // permeability K
    let forch_b = vars.add_scldom_cnst(dom_b, 0.565, "")?;  // Forchheimer coefficient F
    let tort_b = vars.add_scldom_cnst(dom_b, 1.5, "")?;  // tortuosity
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
    let tort_f = vars.add_scldom_cnst(dom_f, 1.0, "")?;
    let area_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;
    let displ_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;
    let dispt_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;

    // mass transfer properties; rho * D; volumetric terms are zero; surface terms are multiplied by a_v
    let diff_b = vars.add_scldom_cnst(dom_b, 1000.0 * 1e-9, "")?;
    let zero_b = vars.add_scldom_cnst(dom_b, 0.0, "")?;
    let sprd_0_b = vars.add_scldom_cnst(dom_b, -1000.0 * 5e-6, "")?;  // surface production coefficient (-rho k_s)
    let ssrc_1_b = vars.add_scldom_func(dom_b, Box::new(|_t, (_x, _y), unkn| 1000.0 * 5e-6 * unkn[0]), vec![conc_0_b], "")?;  // surface source (+rho k_s c_0)
    let diff_f = vars.add_scldom_cnst(dom_f, 1000.0 * 1e-9, "")?;
    let zero_f = vars.add_scldom_cnst(dom_f, 0.0, "")?;

    // boundary values
    let vel_lb = vars.add_vecbnd_cnst(bnd_lb, (1e-3, 0.0), "")?;
    let vel_lf = vars.add_vecbnd_cnst(bnd_lf, (1e-3, 0.0), "")?;
    let vel_w_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_w_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;
    let pres_rb = vars.add_sclbnd_cnst(bnd_rb, 0.0, "")?;
    let pres_rf = vars.add_sclbnd_cnst(bnd_rf, 0.0, "")?;
    let conc_lb = [vars.add_sclbnd_cnst(bnd_lb, 1.0, "")?, vars.add_sclbnd_cnst(bnd_lb, 0.0, "")?];
    let conc_lf = [vars.add_sclbnd_cnst(bnd_lf, 1.0, "")?, vars.add_sclbnd_cnst(bnd_lf, 0.0, "")?];
    let mflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let mflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;

    // steady mass transfer in porous media; arguments of add_mass_dom: component, domain, c, rho*D, s_v, R_v, s_s, R_s
    let mut phys = SteadyMassPore::new(2)?;
    phys.add_mass_dom(0, dom_b, conc_0_b, HashMap::from([(0, diff_b)]), zero_b, zero_b, sprd_0_b, zero_b)?;
    phys.add_mass_dom(1, dom_b, conc_1_b, HashMap::from([(1, diff_b)]), zero_b, zero_b, zero_b, ssrc_1_b)?;
    phys.add_mass_dom(0, dom_f, conc_0_f, HashMap::from([(0, diff_f)]), zero_f, zero_f, zero_f, zero_f)?;
    phys.add_mass_dom(1, dom_f, conc_1_f, HashMap::from([(1, diff_f)]), zero_f, zero_f, zero_f, zero_f)?;
    phys.add_matrix_dom(dom_b, tort_b, area_b, displ_b, dispt_b)?;  // arguments: domain, tau, a_v, alpha_L, alpha_T
    phys.add_matrix_dom(dom_f, tort_f, area_f, displ_f, dispt_f)?;
    phys.add_pore_dom(dom_b, vel_b, pres_b, den_b, visc_b, fce_b, epsl_b, perm_b, forch_b)?;  // arguments: domain, v, p, rho, mu, f, epsl, K, F
    phys.add_pore_dom(dom_f, vel_f, pres_f, den_f, visc_f, fce_f, epsl_f, perm_f, forch_f)?;
    for i in 0..2 {
        phys.add_conc_bnd(i, bnd_lb, conc_lb[i])?;
        phys.add_conc_bnd(i, bnd_lf, conc_lf[i])?;
        phys.add_mout_bnd(i, bnd_rb)?;
        phys.add_mout_bnd(i, bnd_rf)?;
        phys.add_mflx_bnd(i, bnd_b, mflx_b)?;
        phys.add_mflx_bnd(i, bnd_t, mflx_t)?;
        phys.add_mcnt_itf(i, itf_bed)?;  // mass continuity across the interface
    }
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
