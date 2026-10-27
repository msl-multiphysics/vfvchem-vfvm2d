use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Steady-state Navier-Stokes equation with plane Poiseuille flow of water.
/// Run with: `cargo run --release --example flow_analytic`
///
/// Geometry:
/// Rectangular channel defined by lengths along each axis.
/// Reynolds number Re = rho u_max H / mu = 50; analytic pressure gradient dp/dx = -8 mu u_max / H^2 = -3.2 Pa m-1.
///
///                 bnd_t
///       +-----------------------+
/// bnd_l |          dom          | bnd_r
///       +-----------------------+
///                 bnd_b
///
/// Domains:
/// dom - channel (0.02 m x 0.005 m)
///
/// Boundaries:
/// bnd_l - inlet (0.005 m)
/// bnd_r - outlet (0.005 m)
/// bnd_b - bottom wall (0.02 m)
/// bnd_t - top wall (0.02 m)
///
/// Properties:
/// dom - density (rho = 1000.0 kg m-3)
/// dom - viscosity (mu = 0.001 Pa s)
/// dom - body force (f = <0.0, 0.0> N m-3)
///
/// Boundary conditions:
/// [non-constant properties in SI units]
/// bnd_l - parabolic velocity (v = <u, 0.0> m s-1)
/// u = 4 * 0.01 * y * (0.005 - y) / 0.005^2
/// bnd_r - pressure (p = 0.0 Pa)
/// bnd_b - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_t - no-slip (v = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_analytic").unwrap();

    // nodes, points, and variables
    // new_bounds - create boundary nodes and internal points inside a bounding box
    // arguments: x_min, y_min, x_max, y_max, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.02, 0.005, 2.5e-4, 2.5e-4, 1e-3, 100, 42)?;

    // geometry
    // if using new_bounds, the domain is "dom"
    // the left, right, bottom, and top boundaries are "bnd_l", "bnd_r", "bnd_b", and "bnd_t"
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // vector values are given as (x, y) tuples
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., vel_0000.vtu)
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "examples/output_flow_analytic/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_flow_analytic/pres.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;  // density
    let visc = vars.add_scldom_cnst(dom, 0.001, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // non-constant boundary values
    // arguments: boundary, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position (x, y), and *unknown* values
    // - unknown values are given in the same order as in unknown_ids
    // parabolic inlet on y in [0.0, 0.005], maximum 0.01
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 0.01 * y * (0.005 - y) / (0.005 * 0.005);
        (u, 0.0)
    });
    let vel_l = vars.add_vecbnd_func(bnd_l, vel_func, vec![], "")?;  // velocity

    // constant boundary values
    // arguments: boundary, value, output_file
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;  // pressure
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;  // velocity
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0), "")?;  // velocity

    // steady-state momentum transfer solver
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l)?;  // arguments: boundary, v
    phys.add_pres_bnd(bnd_r, pres_r)?;  // arguments: boundary, p
    phys.add_vel_bnd(bnd_b, vel_b)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_t, vel_t)?;  // arguments: boundary, v

    // physics solver
    // RslabLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-6, 1.0)?;
    Ok(())
}
