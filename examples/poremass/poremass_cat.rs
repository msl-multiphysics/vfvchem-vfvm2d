use vfvchem_vfvm2d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Steady-state mass and momentum transfer in a catalytic packed bed with first-order surface reaction 0 -> 1.
/// Run with: `cargo run --release --example poremass_cat`
///
/// Geometry:
/// Packed bed; flow from the left inlet to the right outlet.
///
///                  bnd_t
///       +-----------------------+
/// bnd_l |          dom          | bnd_r
///       +-----------------------+
///                  bnd_b
///
/// Domains:
/// dom - packed bed (0.1 m x 0.02 m)
///
/// Boundaries:
/// bnd_l - inlet (0.02 m)
/// bnd_r - outlet (0.02 m)
/// bnd_b - bottom wall (0.1 m)
/// bnd_t - top wall (0.1 m)
///
/// Mass Transfer
/// Concentrations are per unit fluid mass and diffusion coefficients include the density.
///
/// Properties:
/// [non-constant properties in SI units]
/// dom - diffusion coefficients (rho*D; kg m-1 s-1)
/// [D_00 D_01] = [1.0e-6 0.0]
/// [D_10 D_11]   [0.0    1.0e-6]
/// dom - tortuosity (tau = 1.5)
/// dom - specific surface area (a_v = 6 * (1 - epsl) / diam = 1800.0 m-1)
/// dom - longitudinal dispersivity (alpha_L = diam / 2 = 1.0e-3 m)
/// dom - transverse dispersivity (alpha_T = diam / 10 = 2.0e-4 m)
/// dom - volumetric production coefficients and reaction rates are zero
/// dom - surface production coefficients (kg m-2 s-1; positive if production)
/// [s_s,0] = [-rho * 5.0e-6]
/// [s_s,1]   [0.0]
/// dom - surface reaction rates (mol m-2 s-1; positive if production)
/// [R_s,0] = [0.0]
/// [R_s,1]   [+rho * 5.0e-6 * c_0]
///
/// Boundary conditions:
/// bnd_l - concentration (mol kg-1)
/// [c_0] = [1.0]
/// [c_1]   [0.0]
/// bnd_r - outflow
/// bnd_b - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
/// bnd_t - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
///
/// Momentum Transfer
///
/// Properties:
/// dom - density (rho = 1000.0 kg m-3)
/// dom - viscosity (mu = 0.001 Pa s)
/// dom - body force (f = <0.0, 0.0> N m-3)
/// dom - porosity (epsl = 0.4)
/// dom - permeability (K = 3.95e-9 m2); Kozeny-Carman with particle diameter diam = 2.0e-3 m
/// dom - Forchheimer coefficient (F = 0.565); Ergun
///
/// Boundary conditions:
/// bnd_l - superficial velocity (v = <1.0e-3, 0.0> m s-1)
/// bnd_r - pressure (p = 0.0 Pa)
/// bnd_b - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_t - no-slip (v = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    create_dir_all("examples/output_poremass_cat").unwrap();
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.1, 0.02, 1e-3, 1e-3, 1e-3, 100, 42)?;
    let dom = vars.add_dom("dom")?;
    let bnd_i = vars.add_bnd(dom, "bnd_l")?;
    let bnd_o = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let conc_0 = vars.add_scldom_unkn(dom, 0.0, "examples/output_poremass_cat/conc_0.vtu")?;
    let conc_1 = vars.add_scldom_unkn(dom, 0.0, "examples/output_poremass_cat/conc_1.vtu")?;
    let vel = vars.add_vecdom_unkn(dom, (1e-3, 0.0), "examples/output_poremass_cat/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_poremass_cat/pres.vtu")?;

    // fluid and bed properties
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;  // density
    let visc = vars.add_scldom_cnst(dom, 1e-3, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;  // body force
    let epsl = vars.add_scldom_cnst(dom, 0.4, "")?;  // porosity
    let perm = vars.add_scldom_cnst(dom, 3.95e-9, "")?;  // permeability K
    let forch = vars.add_scldom_cnst(dom, 0.565, "")?;  // Forchheimer coefficient F
    let tort = vars.add_scldom_cnst(dom, 1.5, "")?;  // tortuosity
    let area = vars.add_scldom_cnst(dom, 6.0 * (1.0 - 0.4) / 2e-3, "")?;  // specific surface area a_v
    let displ = vars.add_scldom_cnst(dom, 1e-3, "")?;  // longitudinal dispersivity
    let dispt = vars.add_scldom_cnst(dom, 2e-4, "")?;  // transverse dispersivity

    // mass transfer properties; rho * D; volumetric terms are zero; surface terms are multiplied by a_v
    let diff_00 = vars.add_scldom_cnst(dom, 1000.0 * 1e-9, "")?;
    let diff_11 = vars.add_scldom_cnst(dom, 1000.0 * 1e-9, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let sprd_0 = vars.add_scldom_cnst(dom, -1000.0 * 5e-6, "")?;  // surface production coefficient (-rho k_s)
    let ssrc_1 = vars.add_scldom_func(dom, Box::new(|_t, (_x, _y), unkn| 1000.0 * 5e-6 * unkn[0]), vec![conc_0], "")?;  // surface source (+rho k_s c_0)

    // boundary values
    let vel_i = vars.add_vecbnd_cnst(bnd_i, (1e-3, 0.0), "")?;
    let vel_w = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;
    let conc_0_i = vars.add_sclbnd_cnst(bnd_i, 1.0, "")?;
    let conc_1_i = vars.add_sclbnd_cnst(bnd_i, 0.0, "")?;
    let mflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let mflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;

    // steady mass transfer in porous media; arguments of add_mass_dom: component, domain, c, rho*D, s_v, R_v, s_s, R_s
    let mut phys = SteadyMassPore::new(2)?;
    phys.add_mass_dom(0, dom, conc_0, HashMap::from([(0, diff_00)]), zero, zero, sprd_0, zero)?;
    phys.add_mass_dom(1, dom, conc_1, HashMap::from([(1, diff_11)]), zero, zero, zero, ssrc_1)?;
    phys.add_matrix_dom(dom, tort, area, displ, dispt)?;  // arguments: domain, tau, a_v, alpha_L, alpha_T
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, forch)?;  // arguments: domain, v, p, rho, mu, f, epsl, K, F
    for i in 0..2 {
        phys.add_conc_bnd(i, bnd_i, [conc_0_i, conc_1_i][i])?;
        phys.add_mout_bnd(i, bnd_o)?;
        phys.add_mflx_bnd(i, bnd_b, mflx_b)?;
        phys.add_mflx_bnd(i, bnd_t, mflx_t)?;
    }
    phys.add_vel_bnd(bnd_i, vel_i)?;
    phys.add_pres_bnd(bnd_o, pres_o)?;
    phys.add_vel_bnd(bnd_b, vel_w)?;
    phys.add_vel_bnd(bnd_t, vel_t)?;

    let mut lin = FaerLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-8, 1.0)?;
    
    Ok(())
}
