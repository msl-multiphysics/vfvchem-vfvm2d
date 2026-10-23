use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Transient heat and momentum transfer with backward-facing step flow.
/// Run with: `cargo run --release --example heatflow_step`
///
/// Geometry:
/// Channel (2.4 m x 0.4 m) with a step (0.4 m x 0.2 m) at the bottom left corner.
/// Inlet above the step, outlet at the right end, and walls elsewhere.
///
/// Domains:
/// dom_f - channel
///
/// Boundaries:
/// bnd_i - inlet (x = 0.0 m; 0.2 m < y < 0.4 m)
/// bnd_o - outlet (x = 2.4 m; 0.0 m < y < 0.4 m)
/// bnd_w - walls
///
/// Heat Transfer
///
/// Properties:
/// dom_f - volumetric heat capacity (rho*cp = 1000.0 J m-3 K-1)
/// dom_f - thermal conductivity (k = 0.026 W m-1 K-1)
/// dom_f - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_f - heat source (Q = 0.0 W m-3)
///
/// Initial conditions:
/// dom_f - temperature (T = 300 K)
///
/// Boundary conditions:
/// bnd_i - temperature (T = 300 K)
/// bnd_o - outflow
/// bnd_w - no flux (q = 0.0 W m-2)
///
/// Momentum Transfer
///
/// Properties:
/// dom_f - density (rho = 1.2 kg m-3)
/// dom_f - viscosity (mu = 1.8e-5 Pa s)
/// dom_f - body force (f = <0.0, 0.0> N m-3)
///
/// Initial conditions:
/// dom_f - velocity (v = <1.0, 0.0> m s-1)
/// dom_f - pressure (p = 0.0 Pa)
///
/// Boundary conditions:
/// [non-constant properties in SI units]
/// bnd_i - parabolic velocity (v = <u, 0.0> m s-1)
/// u = 4 * 0.0225 * (y - 0.2) * (0.4 - y) / 0.2^2
/// bnd_o - pressure (p = 0.0 Pa)
/// bnd_w - no-slip (v = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatflow_step").unwrap();

    // nodes, points, and variables
    // new_obj - import boundary, interface, and domain polylines from an obj file
    // arguments: input_file, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", 0.01, 0.02, 0.1, 100, 42)?;

    // geometry
    // domains, boundaries, and interfaces are selected by their group names in the obj file
    // add_bnd - arguments: domain, group_name
    let dom_f = vars.add_dom("dom_f")?;  // channel
    let bnd_i = vars.add_bnd(dom_f, "bnd_i")?;  // inlet
    let bnd_o = vars.add_bnd(dom_f, "bnd_o")?;  // outlet
    let bnd_w = vars.add_bnd(dom_f, "bnd_w")?;  // walls

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // vector values are given as (x, y) tuples
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., temp_0000.vtu)
    let temp = vars.add_scldom_unkn(dom_f, 300.0, "examples/output_heatflow_step/temp.vtu")?;
    let vel = vars.add_vecdom_unkn(dom_f, (1.0, 0.0), "examples/output_heatflow_step/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom_f, 0.0, "examples/output_heatflow_step/pres.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    // heat production is multiplied by the temperature; heat source adds the source directly
    // use heat source by default; use heat production (or both) if heat source does not converge
    let vlcp = vars.add_scldom_cnst(dom_f, 1000.0, "")?;  // volumetric heat capacity (rho * cp)
    let cond = vars.add_scldom_cnst(dom_f, 0.026, "")?;  // thermal conductivity
    let hprd = vars.add_scldom_cnst(dom_f, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc = vars.add_scldom_cnst(dom_f, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let den = vars.add_scldom_cnst(dom_f, 1.2, "")?;  // density
    let visc = vars.add_scldom_cnst(dom_f, 1.8e-5, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom_f, (0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // non-constant boundary values
    // arguments: boundary, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position (x, y), and *unknown* values
    // - unknown values are given in the same order as in unknown_ids
    // parabolic inlet on y in [0.2, 0.4], maximum 0.0225
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 0.0225 * (y - 0.2) * (0.4 - y) / (0.2 * 0.2);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;  // velocity

    // constant boundary values
    // arguments: boundary, value, output_file
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0, "")?;  // temperature
    let hflx_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;  // heat flux (positive if outward; negative if inward)
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;  // pressure
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0), "")?;  // velocity

    // transient heat and momentum transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    // add_hout_bnd - register boundary with heat outflow (heat leaves by advection only)
    // add_hflx_bnd - register boundary with heat flux
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    let mut phys = TransientHeatFlow::new()?;
    phys.add_heat_dom(dom_f, temp, vlcp, cond, hprd, hsrc)?;  // arguments: domain, T, rho*cp, k, s, Q
    phys.add_temp_bnd(bnd_i, temp_i)?;  // arguments: boundary, T
    phys.add_hout_bnd(bnd_o)?;  // arguments: boundary
    phys.add_hflx_bnd(bnd_w, hflx_w)?;  // arguments: boundary, q
    phys.add_flow_dom(dom_f, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;  // arguments: boundary, v
    phys.add_pres_bnd(bnd_o, pres_o)?;  // arguments: boundary, p
    phys.add_vel_bnd(bnd_w, vel_w)?;  // arguments: boundary, v

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, dt, num_ts, num_ts_write, max_iter, tol, damping_factor
    // num_ts_write - writes an output file every num_ts_write time steps; the initial condition is output 0
    // max_iter, tol - maximum number of iterations and tolerance in each time step
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for multiphysics problems, use a tolerance small enough to converge the results
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 0.1, 100, 1, 30, 1e-3, 0.5)?;

    Ok(())
}
