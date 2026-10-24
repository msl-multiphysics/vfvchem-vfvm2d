use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Transient Darcy-Brinkman-Forchheimer equation with wall channeling in a packed bed.
/// Run with: `cargo run --release --example pore_transient`
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
/// [non-constant properties in SI units]
/// dom - density (rho = 1000.0 kg m-3)
/// dom - liquid viscosity (mu = 0.001 Pa s)
/// dom - body force (f = <0.0, 0.0> N m-3)
/// dom - porosity rising near the walls (epsl = 0.4 * (1 + 1.36 * exp(-5 * d_w / diam)))
/// d_w = min(x, 0.02 - x) is the distance from the nearest wall
/// dom - particle diameter (diam = 2.0e-3 m)
///
/// Initial conditions:
/// dom - velocity (v = <0.0, 0.0> m s-1)
/// dom - pressure (p = 0.0 Pa)
///
/// Boundary conditions:
/// bnd_l - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_r - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_b - superficial velocity (v = <0.0, 5.0e-3> m s-1)
/// bnd_t - pressure (p = 0.0 Pa)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_pore_transient").unwrap();

    // nodes, points, and variables
    // porosity is evaluated at the vertices, so h_min resolves the near-wall porosity layer (decay length diam / 5)
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.02, 0.1, 2.5e-4, 1.0e-3, 2.0e-3, 100, 42)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknown domain values
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "examples/output_pore_transient/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_pore_transient/pres.vtu")?;

    // constant domain values
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;  // density
    let visc = vars.add_scldom_cnst(dom, 1.0e-3, "")?;  // liquid viscosity
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;  // body force (den * g; not acceleration g)
    let diam = vars.add_scldom_cnst(dom, 2.0e-3, "")?;  // particle diameter

    // non-constant domain values
    // porosity is evaluated at the vertices; K and F follow from the Ergun relation in the operators
    let epsl_func: SclDomFunc = Box::new(|_t, (x, _y), _unkn| {
        let dw = x.min(0.02 - x);
        0.4 * (1.0 + 1.36 * (-5.0 * dw / 2.0e-3).exp())
    });
    let epsl = vars.add_scldom_func(dom, epsl_func, vec![], "examples/output_pore_transient/epsl.vtu")?;  // porosity

    // constant boundary values
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0), "")?;  // velocity
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0), "")?;  // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 5.0e-3), "")?;  // velocity
    let pres_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;  // pressure

    // transient porous medium flow solver
    let mut phys = TransientPore::new()?;
    phys.add_pore_dom(dom, vel, pres, den, visc, fce, epsl, diam)?;  // arguments: domain, v, p, rho, mu, f, epsl, diam
    phys.add_vel_bnd(bnd_l, vel_l)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_r, vel_r)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_b, vel_b)?;  // arguments: boundary, v
    phys.add_pres_bnd(bnd_t, pres_t)?;  // arguments: boundary, p

    // physics solver
    // arguments: vars, linear_solver, dt, num_ts, num_ts_write, max_iter, tol, damping_factor
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 0.5, 20, 5, 30, 1e-6, 1.0)?;

    Ok(())
}
