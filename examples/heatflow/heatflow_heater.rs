use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Steady-state heat and momentum transfer through a channel between heated solids.
/// Run with: `cargo run --release --example heatflow_heater`
///
/// Geometry:
/// Channel (3.5 m long; 0.5 m tall) that rises over a bump in the bottom wall (1.5 m < x < 2.0 m)
/// into a pocket in the top wall (1.0 m < x < 2.5 m), with solid regions below and above.
///
/// Domains:
/// dom_c - channel
/// dom_b - bottom solid (down to y = -0.5 m)
/// dom_t - top solid (up to y = 1.5 m)
///
/// Boundaries:
/// bnd_i - inlet (x = 0.0 m; 0.0 m < y < 0.5 m)
/// bnd_o - outlet (x = 3.5 m; 0.0 m < y < 0.5 m)
/// bnd_bc - bottom wall of channel
/// bnd_tc - top wall of channel
/// bnd_b - outer boundary of bottom solid
/// bnd_t - outer boundary of top solid
///
/// Interfaces:
/// itf_b - bottom wall of channel (between dom_c and dom_b)
/// itf_t - top wall of channel (between dom_c and dom_t)
///
/// Heat Transfer
///
/// Properties:
/// dom_c - volumetric heat capacity (rho*cp = 4.2e6 J m-3 K-1)
/// dom_c - thermal conductivity (k = 0.6 W m-1 K-1)
/// dom_c - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_c - heat source (Q = 0.0 W m-3)
/// dom_b - thermal conductivity (k = 50.0 W m-1 K-1)
/// dom_b - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_b - heat source (Q = 10.0 W m-3)
/// dom_t - thermal conductivity (k = 50.0 W m-1 K-1)
/// dom_t - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_t - heat source (Q = 10.0 W m-3)
///
/// Boundary conditions:
/// bnd_i - temperature (T = 300 K)
/// bnd_o - outflow
/// bnd_b - temperature (T = 400 K)
/// bnd_t - temperature (T = 400 K)
///
/// Interface conditions:
/// itf_b - temperature and flux continuity
/// itf_t - temperature and flux continuity
///
/// Momentum Transfer
///
/// Properties:
/// dom_c - density (rho = 1000.0 kg m-3)
/// dom_c - viscosity (mu = 0.001 Pa s)
/// dom_c - body force (f = <0.0, 0.0> N m-3)
/// dom_b - unused domain; solid
/// dom_t - unused domain; solid
///
/// Boundary conditions:
/// [non-constant properties in SI units]
/// bnd_i - parabolic velocity (v = <u, 0.0> m s-1)
/// u = 4 * 3.0e-4 * y * (0.5 - y) / 0.5^2
/// bnd_o - pressure (p = 0.0 Pa)
/// bnd_bc - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_tc - no-slip (v = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatflow_heater").unwrap();

    // nodes, points, and variables
    // new_obj - import boundary, interface, and domain polylines from an obj file
    // arguments: input_file, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_obj("examples/obj/obj_heater.obj", 0.01, 0.02, 0.1, 100, 42)?;

    // geometry
    // domains, boundaries, and interfaces are selected by their group names in the obj file
    // add_bnd - arguments: domain, group_name
    // add_itf - arguments: domain_a, domain_b, group_name
    // the channel walls are both interfaces (for heat transfer) and boundaries (for momentum transfer)
    let dom_c = vars.add_dom("dom_c")?;  // channel
    let dom_b = vars.add_dom("dom_b")?;  // bottom solid
    let dom_t = vars.add_dom("dom_t")?;  // top solid
    let bnd_i = vars.add_bnd(dom_c, "bnd_i")?;  // inlet
    let bnd_o = vars.add_bnd(dom_c, "bnd_o")?;  // outlet
    let bnd_bc = vars.add_bnd(dom_c, "bnd_bc")?;  // channel bottom wall (no-slip)
    let bnd_tc = vars.add_bnd(dom_c, "bnd_tc")?;  // channel top wall (no-slip)
    let bnd_b = vars.add_bnd(dom_b, "bnd_b")?;  // outer boundary of bottom solid
    let bnd_t = vars.add_bnd(dom_t, "bnd_t")?;  // outer boundary of top solid
    let itf_b = vars.add_itf(dom_c, dom_b, "itf_b")?;  // channel bottom wall
    let itf_t = vars.add_itf(dom_c, dom_t, "itf_t")?;  // channel top wall

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // vector values are given as (x, y) tuples
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., temp_c_0000.vtu)
    let temp_c = vars.add_scldom_unkn(dom_c, 300.0, "examples/output_heatflow_heater/temp_c.vtu")?;
    let temp_b = vars.add_scldom_unkn(dom_b, 300.0, "examples/output_heatflow_heater/temp_b.vtu")?;
    let temp_t = vars.add_scldom_unkn(dom_t, 300.0, "examples/output_heatflow_heater/temp_t.vtu")?;
    let vel = vars.add_vecdom_unkn(dom_c, (1.0, 0.0), "examples/output_heatflow_heater/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom_c, 0.0, "examples/output_heatflow_heater/pres.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    // heat production is multiplied by the temperature; heat source adds the source directly
    // use heat source by default; use heat production (or both) if heat source does not converge
    // fluid in the channel, heated solids above and below
    let vlcp_c = vars.add_scldom_cnst(dom_c, 4.2e6, "")?;  // volumetric heat capacity (rho * cp)
    let cond_c = vars.add_scldom_cnst(dom_c, 0.6, "")?;  // thermal conductivity
    let hprd_c = vars.add_scldom_cnst(dom_c, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc_c = vars.add_scldom_cnst(dom_c, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let vlcp_b = vars.add_scldom_cnst(dom_b, 550.0, "")?;  // unused for steady conduction-only solid
    let cond_b = vars.add_scldom_cnst(dom_b, 50.0, "")?;  // thermal conductivity
    let hprd_b = vars.add_scldom_cnst(dom_b, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc_b = vars.add_scldom_cnst(dom_b, 10.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let vlcp_t = vars.add_scldom_cnst(dom_t, 550.0, "")?;  // unused for steady conduction-only solid
    let cond_t = vars.add_scldom_cnst(dom_t, 50.0, "")?;  // thermal conductivity
    let hprd_t = vars.add_scldom_cnst(dom_t, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc_t = vars.add_scldom_cnst(dom_t, 10.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let den = vars.add_scldom_cnst(dom_c, 1000.0, "")?;  // density
    let visc = vars.add_scldom_cnst(dom_c, 0.001, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom_c, (0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // non-constant boundary values
    // arguments: boundary, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position (x, y), and *unknown* values
    // - time is zero for steady-state problems
    // - unknown values are given in the same order as in unknown_ids
    // parabolic inlet on y in [0.0, 0.5], maximum 3e-4
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 3e-4 * y * (0.5 - y) / (0.5 * 0.5);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;  // velocity

    // constant boundary values
    // arguments: boundary, value, output_file
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0, "")?;  // temperature
    let temp_bb = vars.add_sclbnd_cnst(bnd_b, 400.0, "")?;  // temperature
    let temp_tt = vars.add_sclbnd_cnst(bnd_t, 400.0, "")?;  // temperature
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;  // pressure
    let vel_bc = vars.add_vecbnd_cnst(bnd_bc, (0.0, 0.0), "")?;  // velocity
    let vel_tc = vars.add_vecbnd_cnst(bnd_tc, (0.0, 0.0), "")?;  // velocity

    // steady-state heat and momentum transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    // add_hout_bnd - register boundary with heat outflow (heat leaves by advection only)
    // add_hcnt_itf - register interface with temperature and flux continuity
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom_c, temp_c, vlcp_c, cond_c, hprd_c, hsrc_c)?;  // arguments: domain, T, rho*cp, k, s, Q
    phys.add_heat_dom(dom_b, temp_b, vlcp_b, cond_b, hprd_b, hsrc_b)?;  // arguments: domain, T, rho*cp, k, s, Q
    phys.add_heat_dom(dom_t, temp_t, vlcp_t, cond_t, hprd_t, hsrc_t)?;  // arguments: domain, T, rho*cp, k, s, Q
    phys.add_temp_bnd(bnd_i, temp_i)?;  // arguments: boundary, T
    phys.add_hout_bnd(bnd_o)?;  // arguments: boundary
    phys.add_temp_bnd(bnd_b, temp_bb)?;  // arguments: boundary, T
    phys.add_temp_bnd(bnd_t, temp_tt)?;  // arguments: boundary, T
    phys.add_hcnt_itf(itf_b)?;  // arguments: interface
    phys.add_hcnt_itf(itf_t)?;  // arguments: interface
    phys.add_flow_dom(dom_c, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;  // arguments: boundary, v
    phys.add_pres_bnd(bnd_o, pres_o)?;  // arguments: boundary, p
    phys.add_vel_bnd(bnd_bc, vel_bc)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_tc, vel_tc)?;  // arguments: boundary, v

    // physics solver
    // FaerLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for multiphysics problems, use a tolerance small enough to converge the results
    let mut lin = FaerLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 0.5)?;

    Ok(())
}
