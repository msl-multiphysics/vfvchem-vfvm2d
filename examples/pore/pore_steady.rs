use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Steady-state Darcy-Brinkman-Forchheimer equation with water flowing through a packed bed.
/// Run with: `cargo run --release --example pore_steady`
///
/// Geometry:
/// Packed column; upward flow from the bottom inlet to the top outlet.
///
///          bnd_t
///       +---------+
///       |         |
/// bnd_l |   dom   | bnd_r
///       |         |
///       +---------+
///          bnd_b
///
/// Domains:
/// dom - packed bed (0.02 m x 0.1 m)
///
/// Boundaries:
/// bnd_l - left wall (0.1 m)
/// bnd_r - right wall (0.1 m)
/// bnd_b - inlet (0.02 m)
/// bnd_t - outlet (0.02 m)
///
/// Properties:
/// dom - density (rho = 1000.0 kg m-3)
/// dom - viscosity (mu = 0.001 Pa s)
/// dom - body force (f = <0.0, -rho * 9.81> N m-3)
/// dom - porosity (epsl = 0.4)
/// dom - permeability (K = 4.7407e-9 m2); Ergun with particle diameter diam = 2.0e-3 m
/// dom - Forchheimer coefficient (F = 0.56481); Ergun
///
/// Boundary conditions:
/// bnd_l - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_r - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_b - superficial velocity (v = <0.0, 5.0e-3> m s-1); Re_p = rho * U * diam / mu = 10
/// bnd_t - pressure (p = 0.0 Pa)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_pore_steady").unwrap();

    // parameters
    let (den_c, visc_c, epsl_c, vel_c, grav) = (1000.0, 1.0e-3, 0.4, 5.0e-3, 9.81);

    // permeability and Forchheimer coefficient; computed by hand from the Ergun relation (diam = 2.0e-3 m, epsl = 0.4)
    // K = diam^2 * epsl^3 / (150 * (1 - epsl)^2) = 4.7407e-9 m2
    // F = 1.75 / sqrt(150 * epsl^3) = 0.56481
    let (perm_c, forch_c) = (4.7407e-9, 0.56481);
    let (wid, hgt) = (0.02, 0.1);

    // nodes, points, and variables
    // new_bounds - rectangle; arguments: x_min, y_min, x_max, y_max, h_min, h_max, ref_len, max_try, seed
    let mut vars = Variables::new_bounds(0.0, 0.0, wid, hgt, 1.0e-3, 1.0e-3, 1.0e-3, 100, 42)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknown domain values
    // vel is the superficial (Darcy) velocity
    let vel = vars.add_vecdom_unkn(dom, (0.0, vel_c), "examples/output_pore_steady/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_pore_steady/pres.vtu")?;

    // constant domain values
    let den = vars.add_scldom_cnst(dom, den_c, "")?;  // density
    let visc = vars.add_scldom_cnst(dom, visc_c, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom, (0.0, -den_c * grav), "")?;  // body force (den * g; not acceleration g)
    let epsl = vars.add_scldom_cnst(dom, epsl_c, "")?;  // porosity
    let perm = vars.add_scldom_cnst(dom, perm_c, "")?;  // permeability
    let forch = vars.add_scldom_cnst(dom, forch_c, "")?;  // Forchheimer coefficient

    // constant boundary values
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0), "")?;  // velocity
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0), "")?;  // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, vel_c), "")?;  // velocity
    let pres_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;  // pressure

    // steady-state porous medium flow solver
    // add_pore_dom - register domain with porous medium flow; permeability K and Forchheimer coefficient F are user inputs
    // mu is the liquid viscosity
    let mut phys = SteadyPore::new()?;
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, perm, forch)?;  // arguments: domain, v, p, rho, mu, f, epsl, K, F
    phys.add_vel_bnd(bnd_l, vel_l)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_r, vel_r)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_b, vel_b)?;  // arguments: boundary, v
    phys.add_pres_bnd(bnd_t, pres_t)?;  // arguments: boundary, p

    // physics solver
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-8, 1.0)?;

    Ok(())
}
