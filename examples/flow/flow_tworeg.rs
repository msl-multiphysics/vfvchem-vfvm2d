use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Steady-state Navier-Stokes equation with lid-driven cavity flow in two regions joined by an interface.
/// Run with: `cargo run --release --example flow_tworeg`
///
/// Geometry:
/// Square domain split along its diagonal into two triangular regions.
///
///            bnd_t
///       +-------------+
///       |  dom_up   / |
///       |         /   |
/// bnd_l |    itf_d    | bnd_r
///       |   /         |
///       | /   dom_lo  |
///       +-------------+
///            bnd_b
///
/// Domains:
/// dom_lo - lower right region (triangle with legs of 1.0 m)
/// dom_up - upper left region (triangle with legs of 1.0 m)
///
/// Boundaries:
/// bnd_l - left boundary of dom_up (1.0 m)
/// bnd_r - right boundary of dom_lo (1.0 m)
/// bnd_b - bottom boundary of dom_lo (1.0 m)
/// bnd_t - top boundary of dom_up (1.0 m)
///
/// Interfaces:
/// itf_d - diagonal between dom_lo and dom_up
///
/// Properties:
/// dom_lo - density (rho = 1000.0 kg m-3)
/// dom_lo - viscosity (mu = 0.001 Pa s)
/// dom_lo - body force (f = <0.0, 0.0> N m-3)
/// dom_up - density (rho = 1000.0 kg m-3)
/// dom_up - viscosity (mu = 0.001 Pa s)
/// dom_up - body force (f = <0.0, 0.0> N m-3)
///
/// Boundary conditions:
/// bnd_l - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_r - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_b - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_t - velocity (v = <4.0e-4, 0.0> m s-1)
/// Zero reference pressure at node 0 of dom_lo
///
/// Interface conditions:
/// itf_d - velocity and pressure continuity
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_tworeg").unwrap();

    // nodes, points, and variables
    // new_obj - import boundary, interface, and domain polylines from an obj file
    // arguments: input_file, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_obj("examples/obj/obj_tworeg.obj", 0.02, 0.02, 0.02, 100, 42)?;

    // geometry
    // domains, boundaries, and interfaces are selected by their group names in the obj file
    // add_bnd - arguments: domain, group_name
    // add_itf - arguments: domain_a, domain_b, group_name
    let dom_lo = vars.add_dom("dom_lo")?;  // lower right
    let dom_up = vars.add_dom("dom_up")?;  // upper left
    let bnd_b = vars.add_bnd(dom_lo, "bnd_b")?;  // bottom
    let bnd_r = vars.add_bnd(dom_lo, "bnd_r")?;  // right
    let bnd_t = vars.add_bnd(dom_up, "bnd_t")?;  // top
    let bnd_l = vars.add_bnd(dom_up, "bnd_l")?;  // left
    let itf_d = vars.add_itf(dom_lo, dom_up, "itf_d")?;  // diagonal

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // vector values are given as (x, y) tuples
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., vel_lo_0000.vtu)
    let vel_lo = vars.add_vecdom_unkn(dom_lo, (0.0, 0.0), "examples/output_flow_tworeg/vel_lo.vtu")?;
    let vel_up = vars.add_vecdom_unkn(dom_up, (0.0, 0.0), "examples/output_flow_tworeg/vel_up.vtu")?;
    let pres_lo = vars.add_scldom_unkn(dom_lo, 0.0, "examples/output_flow_tworeg/pres_lo.vtu")?;
    let pres_up = vars.add_scldom_unkn(dom_up, 0.0, "examples/output_flow_tworeg/pres_up.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    let den_lo = vars.add_scldom_cnst(dom_lo, 1000.0, "")?;  // density
    let den_up = vars.add_scldom_cnst(dom_up, 1000.0, "")?;  // density
    let visc_lo = vars.add_scldom_cnst(dom_lo, 0.001, "")?;  // viscosity
    let visc_up = vars.add_scldom_cnst(dom_up, 0.001, "")?;  // viscosity
    let fce_lo = vars.add_vecdom_cnst(dom_lo, (0.0, 0.0), "")?;  // body force (den * g; not acceleration g)
    let fce_up = vars.add_vecdom_cnst(dom_up, (0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // constant boundary values
    // arguments: boundary, value, output_file
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0), "")?;  // velocity
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0), "")?;  // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;  // velocity
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (4e-4, 0.0), "")?;  // velocity

    // steady-state momentum transfer solver
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_fcnt_itf - register interface with velocity and pressure continuity
    // set_pres_ref - set reference pressure at a node; needed if no boundary sets the pressure
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom_lo, vel_lo, pres_lo, den_lo, visc_lo, fce_lo)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_flow_dom(dom_up, vel_up, pres_up, den_up, visc_up, fce_up)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_r, vel_r)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_b, vel_b)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_t, vel_t)?;  // arguments: boundary, v
    phys.add_fcnt_itf(itf_d)?;  // arguments: interface
    phys.set_pres_ref(dom_lo, 0, 0.0)?;  // arguments: domain, node_index, p

    // physics solver
    // FaerLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = FaerLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;

    Ok(())
}
