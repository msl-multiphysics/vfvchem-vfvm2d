use vfvchem_vfvm2d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Transient mass and momentum transfer with tracer injection into water-saturated sand and first-order reaction 0 -> 1 in the pore fluid.
/// Run with: `cargo run --release --example poremass_hom`
///
/// Geometry:
/// Sand layer; flow from the left inlet to the right outlet.
///
///                  bnd_t
///       +-----------------------+
/// bnd_l |          dom          | bnd_r
///       +-----------------------+
///                  bnd_b
///
/// Domains:
/// dom - sand (0.2 m x 0.05 m)
///
/// Boundaries:
/// bnd_l - inlet (0.05 m)
/// bnd_r - outlet (0.05 m)
/// bnd_b - bottom wall (0.2 m)
/// bnd_t - top wall (0.2 m)
///
/// Mass Transfer
/// Concentrations are per unit fluid mass and diffusion coefficients include the density.
/// Volumetric terms are per unit fluid volume and are multiplied by the porosity.
///
/// Properties:
/// [non-constant properties in SI units]
/// dom - diffusion coefficients (rho*D; kg m-1 s-1)
/// [D_00 D_01] = [1.0e-6 0.0]
/// [D_10 D_11]   [0.0    1.0e-6]
/// dom - tortuosity (tau = 1.4)
/// dom - specific surface area (a_v = 0.0 m-1; no surface reaction)
/// dom - longitudinal dispersivity (alpha_L = 2.5e-3 m)
/// dom - transverse dispersivity (alpha_T = 2.5e-4 m)
/// dom - volumetric production coefficients (kg m-3 s-1; positive if production)
/// [s_v,0] = [-rho * 0.01]
/// [s_v,1]   [0.0]
/// dom - volumetric reaction rates (mol m-3 s-1; positive if production)
/// [R_v,0] = [0.0]
/// [R_v,1]   [+rho * 0.01 * c_0]
///
/// Initial conditions:
/// dom - concentration (mol kg-1)
/// [c_0] = [0.0]
/// [c_1]   [0.0]
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
/// dom - porosity (epsl = 0.35)
/// dom - permeability (K = 1.41e-10 m2); Kozeny-Carman with particle diameter diam = 5.0e-4 m
/// dom - Forchheimer coefficient (F = 0.690); Ergun
///
/// Initial conditions:
/// dom - velocity (v = <5.0e-4, 0.0> m s-1)
/// dom - pressure (p = 0.0 Pa)
///
/// Boundary conditions:
/// bnd_l - superficial velocity (v = <5.0e-4, 0.0> m s-1)
/// bnd_r - pressure (p = 0.0 Pa)
/// bnd_b - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_t - no-slip (v = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    create_dir_all("examples/output_poremass_hom").unwrap();
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.2, 0.05, 2.5e-3, 2.5e-3, 2.5e-3, 100, 42)?;
    let dom = vars.add_dom("dom")?;
    let bnd_i = vars.add_bnd(dom, "bnd_l")?;
    let bnd_o = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let conc_0 = vars.add_scldom_unkn(dom, 0.0, "examples/output_poremass_hom/conc_0.vtu")?;
    let conc_1 = vars.add_scldom_unkn(dom, 0.0, "examples/output_poremass_hom/conc_1.vtu")?;
    let vel = vars.add_vecdom_unkn(dom, (5e-4, 0.0), "examples/output_poremass_hom/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_poremass_hom/pres.vtu")?;

    // fluid and matrix properties
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 1e-3, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    let epsl = vars.add_scldom_cnst(dom, 0.35, "")?;
    let perm = vars.add_scldom_cnst(dom, 1.41e-10, "")?;
    let forch = vars.add_scldom_cnst(dom, 0.690, "")?;
    let tort = vars.add_scldom_cnst(dom, 1.4, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;  // a_v and the surface terms
    let displ = vars.add_scldom_cnst(dom, 2.5e-3, "")?;
    let dispt = vars.add_scldom_cnst(dom, 2.5e-4, "")?;

    // mass transfer properties; rho * D; volumetric terms are multiplied by epsl
    let diff_00 = vars.add_scldom_cnst(dom, 1000.0 * 1e-9, "")?;
    let diff_11 = vars.add_scldom_cnst(dom, 1000.0 * 1e-9, "")?;
    let mprd_0 = vars.add_scldom_cnst(dom, -1000.0 * 0.01, "")?;  // production coefficient (-rho k)
    let msrc_1 = vars.add_scldom_func(dom, Box::new(|_t, (_x, _y), unkn| 1000.0 * 0.01 * unkn[0]), vec![conc_0], "")?;  // source (+rho k c_0)

    // boundary values
    let vel_i = vars.add_vecbnd_cnst(bnd_i, (5e-4, 0.0), "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;
    let conc_0_i = vars.add_sclbnd_cnst(bnd_i, 1.0, "")?;
    let conc_1_i = vars.add_sclbnd_cnst(bnd_i, 0.0, "")?;
    let mflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let mflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;

    // transient mass transfer in porous media; arguments of add_mass_dom: component, domain, c, rho*D, s_v, R_v, s_s, R_s
    let mut phys = TransientMassPore::new(2)?;
    phys.add_mass_dom(0, dom, conc_0, HashMap::from([(0, diff_00)]), mprd_0, zero, zero, zero)?;
    phys.add_mass_dom(1, dom, conc_1, HashMap::from([(1, diff_11)]), zero, msrc_1, zero, zero)?;
    phys.add_matrix_dom(dom, tort, zero, displ, dispt)?;  // arguments: domain, tau, a_v, alpha_L, alpha_T
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, forch)?;  // arguments: domain, v, p, rho, mu, f, epsl, K, F
    for i in 0..2 {
        phys.add_conc_bnd(i, bnd_i, [conc_0_i, conc_1_i][i])?;
        phys.add_mout_bnd(i, bnd_o)?;
        phys.add_mflx_bnd(i, bnd_b, mflx_b)?;
        phys.add_mflx_bnd(i, bnd_t, mflx_t)?;
    }
    phys.add_vel_bnd(bnd_i, vel_i)?;
    phys.add_pres_bnd(bnd_o, pres_o)?;
    phys.add_vel_bnd(bnd_b, vel_b)?;
    phys.add_vel_bnd(bnd_t, vel_t)?;

    // 50 s steps for 100 steps (pore velocity 1.4 mm s-1), write every 20; the deferred cross dispersion converges slowly in the first steps
    let mut lin = FaerLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 50.0, 100, 20, 300, 1e-6, 1.0)?;
    
    Ok(())
}
