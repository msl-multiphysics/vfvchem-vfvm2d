use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Steady-state heat and momentum transfer through a channel with four heated cylinders.
/// Run with: `cargo run --release --example heatflow_fourcyl`
///
/// Geometry:
/// Channel (0.2 m x 0.05 m) with four cylinders (radius = 0.005 m) in a 2 x 2 arrangement.
/// Cylinder centers at <0.1125, 0.0125>, <0.1375, 0.0125>, <0.1125, 0.0375>, and <0.1375, 0.0375> m.
///
/// Domains:
/// dom_c - channel
/// dom_h - cylinders
///
/// Boundaries:
/// bnd_i - inlet (x = 0.0 m; 0.0 m < y < 0.05 m)
/// bnd_o - outlet (x = 0.2 m; 0.0 m < y < 0.05 m)
/// bnd_w - bottom and top walls of channel
/// bnd_h - cylinder walls
///
/// Interfaces:
/// itf_h - cylinder walls (between dom_c and dom_h)
///
/// Heat Transfer
///
/// Properties:
/// dom_c - volumetric heat capacity (rho*cp = 1000.0 J m-3 K-1)
/// dom_c - thermal conductivity (k = 0.026 W m-1 K-1)
/// dom_c - heat source (Q = 0.0 W m-3)
/// dom_h - thermal conductivity (k = 50.0 W m-1 K-1)
/// dom_h - heat source (Q = 1.0e6 W m-3)
///
/// Boundary conditions:
/// bnd_i - temperature (T = 300 K)
/// bnd_o - outflow
/// bnd_w - no flux (q = 0.0 W m-2)
///
/// Interface conditions:
/// itf_h - temperature and flux continuity
///
/// Momentum Transfer
///
/// Properties:
/// dom_c - density (rho = 1.2 kg m-3)
/// dom_c - viscosity (mu = 1.8e-5 Pa s)
/// dom_c - body force (f = <0.0, 0.0> N m-3)
/// dom_h - unused domain; solid
///
/// Boundary conditions:
/// [non-constant properties in SI units]
/// bnd_i - parabolic velocity (v = <u, 0.0> m s-1)
/// u = 4 * 0.3 * y * (0.05 - y) / 0.05^2
/// bnd_o - pressure (p = 0.0 Pa)
/// bnd_w - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_h - no-slip (v = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatflow_fourcyl").unwrap();

    // nodes, points, and variables
    // new_obj - import boundary, interface, and domain polylines from an obj file
    // arguments: input_file, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_obj("examples/obj/obj_fourcyl.obj", 0.001, 0.001, 0.001, 100, 42)?;

    // geometry
    // domains, boundaries, and interfaces are selected by their group names in the obj file
    // add_bnd - arguments: domain, group_name
    // add_itf - arguments: domain_a, domain_b, group_name
    // the cylinders share one domain, one interface, and one no-slip boundary
    // the cylinder walls are both an interface (for heat transfer) and a boundary (for momentum transfer)
    let dom_c = vars.add_dom("dom_c")?;  // channel
    let dom_h = vars.add_dom("dom_h")?;  // cylinders
    let bnd_i = vars.add_bnd(dom_c, "bnd_i")?;  // inlet
    let bnd_o = vars.add_bnd(dom_c, "bnd_o")?;  // outlet
    let bnd_w = vars.add_bnd(dom_c, "bnd_w")?;  // channel walls
    let bnd_h = vars.add_bnd(dom_c, "bnd_h")?;  // cylinder walls (no-slip)
    let itf_h = vars.add_itf(dom_c, dom_h, "itf_h")?;  // cylinder walls

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // vector values are given as (x, y) tuples
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., temp_c_0000.vtu)
    let temp_c = vars.add_scldom_unkn(dom_c, 300.0, "examples/output_heatflow_fourcyl/temp_c.vtu")?;
    let temp_h = vars.add_scldom_unkn(dom_h, 300.0, "examples/output_heatflow_fourcyl/temp_h.vtu")?;
    let vel = vars.add_vecdom_unkn(dom_c, (1.0, 0.0), "examples/output_heatflow_fourcyl/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom_c, 0.0, "examples/output_heatflow_fourcyl/pres.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    // fluid in the channel, heated solid cylinders
    let vlcp_c = vars.add_scldom_cnst(dom_c, 1000.0, "")?;  // volumetric heat capacity (rho * cp)
    let cond_c = vars.add_scldom_cnst(dom_c, 0.026, "")?;  // thermal conductivity
    let hsrc_c = vars.add_scldom_cnst(dom_c, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let vlcp_h = vars.add_scldom_cnst(dom_h, 550.0, "")?;  // unused for steady conduction-only solid
    let cond_h = vars.add_scldom_cnst(dom_h, 50.0, "")?;  // thermal conductivity
    let hsrc_h = vars.add_scldom_cnst(dom_h, 1.0e6, "")?;  // heat source (Q; positive if source; negative if sink)
    let den = vars.add_scldom_cnst(dom_c, 1.2, "")?;  // density
    let visc = vars.add_scldom_cnst(dom_c, 1.8e-5, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom_c, (0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // non-constant boundary values
    // arguments: boundary, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position (x, y), and *unknown* values
    // - time is zero for steady-state problems
    // - unknown values are given in the same order as in unknown_ids
    // parabolic inlet on y in [0.0, 0.05], average 0.2, maximum 0.3
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 0.3 * y * (0.05 - y) / (0.05 * 0.05);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;  // velocity

    // constant boundary values
    // arguments: boundary, value, output_file
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0, "")?;  // temperature
    let hflx_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;  // heat flux (positive if outward; negative if inward)
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;  // pressure
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0), "")?;  // velocity
    let vel_h = vars.add_vecbnd_cnst(bnd_h, (0.0, 0.0), "")?;  // velocity

    // steady-state heat and momentum transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    // add_hout_bnd - register boundary with heat outflow (heat leaves by advection only)
    // add_hflx_bnd - register boundary with heat flux
    // add_hcnt_itf - register interface with temperature and flux continuity
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom_c, temp_c, vlcp_c, cond_c, hsrc_c)?;  // arguments: domain, T, rho*cp, k, Q
    phys.add_heat_dom(dom_h, temp_h, vlcp_h, cond_h, hsrc_h)?;  // arguments: domain, T, rho*cp, k, Q
    phys.add_temp_bnd(bnd_i, temp_i)?;  // arguments: boundary, T
    phys.add_hout_bnd(bnd_o)?;  // arguments: boundary
    phys.add_hflx_bnd(bnd_w, hflx_w)?;  // arguments: boundary, q
    phys.add_hcnt_itf(itf_h)?;  // arguments: interface
    phys.add_flow_dom(dom_c, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;  // arguments: boundary, v
    phys.add_pres_bnd(bnd_o, pres_o)?;  // arguments: boundary, p
    phys.add_vel_bnd(bnd_w, vel_w)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_h, vel_h)?;  // arguments: boundary, v

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for multiphysics problems, use a tolerance small enough to converge the results
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 0.5)?;

    Ok(())
}
