use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Convergence scaling options for lid-driven cavity flow of water at absolute pressure.
/// Run with: `cargo run --release --example feat_scale`
///
/// The iteration stops when the change in every unknown field, measured against the field's own scale, is below tol.
/// default   - per-field scaling; the pressure is measured against its range and rho v^2, not its absolute value
/// floor     - set_scldom_scale and set_vecdom_scale set lower bounds on the field scales
/// scale_off - unscaled |x_new - x_old| / |x_old|, dominated by the absolute pressure, so the iteration stops early
///
/// Geometry:
/// Square domain with smaller point spacing near the walls.
/// Reynolds number Re = rho u_lid L / mu = 100.
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
/// dom - main domain (0.01 m x 0.01 m)
///
/// Boundaries:
/// bnd_l - left boundary (0.01 m)
/// bnd_r - right boundary (0.01 m)
/// bnd_b - bottom boundary (0.01 m)
/// bnd_t - top boundary (0.01 m)
///
/// Properties:
/// dom - density (rho = 1000.0 kg m-3)
/// dom - viscosity (mu = 0.001 Pa s)
/// dom - body force (f = <0.0, 0.0> N m-3)
///
/// Boundary conditions:
/// bnd_l - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_r - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_b - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_t - velocity (v = <0.01, 0.0> m s-1)
/// Reference pressure (p = 101325.0 Pa) at node 0 of dom
///
fn run(mode: &str) -> Result<(), VFVChemError> {
    // output directory
    let out = format!("examples/output_feat_scale/{mode}");
    create_dir_all(&out).unwrap();

    // nodes, points, and variables
    // new_bounds - create boundary nodes and internal points inside a bounding box
    // arguments: x_min, y_min, x_max, y_max, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    // points are refined near the walls, where the velocity gradients are largest
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.01, 0.01, 2.5e-4, 5e-4, 2e-3, 100, 42)?;

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
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), &format!("{out}/vel.vtu"))?;
    let pres = vars.add_scldom_unkn(dom, 101325.0, &format!("{out}/pres.vtu"))?;

    // constant domain values
    // arguments: domain, value, output_file
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;  // density
    let visc = vars.add_scldom_cnst(dom, 0.001, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // constant boundary values
    // arguments: boundary, value, output_file
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0), "")?;  // velocity
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0), "")?;  // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;  // velocity
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.01, 0.0), "")?;  // velocity

    // steady-state momentum transfer solver
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // set_pres_ref - set reference pressure at a node; needed if no boundary sets the pressure
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_r, vel_r)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_b, vel_b)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_t, vel_t)?;  // arguments: boundary, v
    phys.set_pres_ref(dom, 0, 101325.0)?;  // arguments: domain, node_index, p

    // convergence scaling
    // set_scldom_scale, set_vecdom_scale - set a lower bound on the scale of a field; arguments: field, scale
    // scale_off - if true, the unscaled change |x_new - x_old| / |x_old| is used instead
    if mode == "floor" {
        vars.set_scldom_scale(pres, 0.1)?;  // pressure difference (Pa)
        vars.set_vecdom_scale(vel, 0.01)?;  // velocity (m s-1)
    } else if mode == "scale_off" {
        vars.scale_off = true;
    }

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-6, 0.8)?;

    Ok(())
}

fn main() -> Result<(), VFVChemError> {
    // compare the number of iterations of each mode in the printed residuals
    run("default")?;
    run("floor")?;
    run("scale_off")?;

    Ok(())
}
