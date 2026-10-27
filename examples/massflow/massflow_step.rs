use vfvchem_vfvm2d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Transient mass and momentum transfer with backward-facing step flow and first-order reaction 0 -> 1.
/// Run with: `cargo run --release --example massflow_step`
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
/// Mass Transfer
/// Concentrations are per unit mass and diffusion coefficients include the density on domains with flow.
///
/// Properties:
/// [non-constant properties in SI units]
/// dom_f - diffusion coefficients (rho*D; kg m-1 s-1)
/// [D_00 D_01] = [1.2e-4 0.0]
/// [D_10 D_11]   [0.0    2.4e-4]
/// dom_f - production coefficients (s-1; positive if production)
/// [s_0] = [0.0]
/// [s_1]   [0.0]
/// dom_f - reaction rates (mol m-3 s-1; positive if production)
/// [R_0] = [-1.2 * c_0]
/// [R_1]   [+1.2 * c_0]
///
/// Initial conditions:
/// dom_f - concentration (mol kg-1)
/// [c_0] = [0.0]
/// [c_1]   [0.0]
///
/// Boundary conditions:
/// bnd_i - concentration (mol kg-1)
/// [c_0] = [1.0 / 1.2]
/// [c_1]   [0.0      ]
/// bnd_o - outflow
/// bnd_w - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
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
    create_dir_all("examples/output_massflow_step").unwrap();

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
    // output index is automatically appended to end of output file name (e.g., conc_0_0000.vtu)
    let conc_0 = vars.add_scldom_unkn(dom_f, 0.0, "examples/output_massflow_step/conc_0.vtu")?;
    let conc_1 = vars.add_scldom_unkn(dom_f, 0.0, "examples/output_massflow_step/conc_1.vtu")?;
    let vel = vars.add_vecdom_unkn(dom_f, (1.0, 0.0), "examples/output_massflow_step/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom_f, 0.0, "examples/output_massflow_step/pres.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    // no need to declare zero diffusion coefficients
    // diffusion coefficients include the density on domains with flow
    // mass production is multiplied by the concentration; mass source adds the source directly
    // use mass source by default; use mass production (or both) if mass source does not converge
    let diff_00 = vars.add_scldom_cnst(dom_f, 1.2e-4, "")?;  // rho * D_00
    let diff_11 = vars.add_scldom_cnst(dom_f, 2.4e-4, "")?;  // rho * D_11
    let mprd_0 = vars.add_scldom_cnst(dom_f, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let mprd_1 = vars.add_scldom_cnst(dom_f, 0.0, "")?;  // mass production coefficient (s_1; positive if source; negative if sink)
    let den = vars.add_scldom_cnst(dom_f, 1.2, "")?;  // density
    let visc = vars.add_scldom_cnst(dom_f, 1.8e-5, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom_f, (0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // non-constant domain values
    // arguments: domain, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position (x, y), and *unknown* values
    // - unknown values are given in the same order as in unknown_ids
    // species 0 converts to species 1, first order in the reactant
    // the reaction rate is per unit volume, so the per unit mass concentration is multiplied by the density
    let msrc_0_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| { -1.2 * unkn[0] });  // unkn[0] is c_0
    let msrc_1_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {  1.2 * unkn[0] });  // unkn[0] is c_0
    let msrc_0 = vars.add_scldom_func(dom_f, msrc_0_func, vec![conc_0], "")?;  // mass source (R_0; positive if source; negative if sink)
    let msrc_1 = vars.add_scldom_func(dom_f, msrc_1_func, vec![conc_0], "")?;  // mass source (R_1; positive if source; negative if sink)

    // non-constant boundary values
    // arguments: boundary, value_func, unknown_ids, output_file
    // parabolic inlet on y in [0.2, 0.4], maximum 0.0225
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 0.0225 * (y - 0.2) * (0.4 - y) / (0.2 * 0.2);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;  // velocity

    // constant boundary values
    // arguments: boundary, value, output_file
    let conc_0_i = vars.add_sclbnd_cnst(bnd_i, 1.0 / 1.2, "")?;  // concentration (per unit mass)
    let conc_1_i = vars.add_sclbnd_cnst(bnd_i, 0.0, "")?;  // concentration
    let mflx_0_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;  // molar flux (positive if outward; negative if inward)
    let mflx_1_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;  // molar flux
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;  // pressure
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0), "")?;  // velocity

    // transient mass and momentum transfer solver
    // TransientMassFlow::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mout_bnd - register boundary with mass outflow (mass leaves by advection only)
    // add_mflx_bnd - register boundary with molar flux
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    // diffusion coefficients must be collected in a hashmap with driving concentration index as key
    // on domains with flow, c is per unit mass and D includes the density (rho * D)
    let mut phys = TransientMassFlow::new(2)?;
    phys.add_mass_dom(0, dom_f, conc_0, HashMap::from([(0, diff_00)]), mprd_0, msrc_0)?;  // arguments: component, domain, c, rho*D, s, R
    phys.add_mass_dom(1, dom_f, conc_1, HashMap::from([(1, diff_11)]), mprd_1, msrc_1)?;  // arguments: component, domain, c, rho*D, s, R
    phys.add_conc_bnd(0, bnd_i, conc_0_i)?;  // arguments: component, boundary, c
    phys.add_conc_bnd(1, bnd_i, conc_1_i)?;  // arguments: component, boundary, c
    phys.add_mout_bnd(0, bnd_o)?;  // arguments: component, boundary
    phys.add_mout_bnd(1, bnd_o)?;  // arguments: component, boundary
    phys.add_mflx_bnd(0, bnd_w, mflx_0_w)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(1, bnd_w, mflx_1_w)?;  // arguments: component, boundary, N
    phys.add_flow_dom(dom_f, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;  // arguments: boundary, v
    phys.add_pres_bnd(bnd_o, pres_o)?;  // arguments: boundary, p
    phys.add_vel_bnd(bnd_w, vel_w)?;  // arguments: boundary, v

    // physics solver
    // RslabLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, dt, num_ts, num_ts_write, max_iter, tol, damping_factor
    // num_ts_write - writes an output file every num_ts_write time steps; the initial condition is output 0
    // max_iter, tol - maximum number of iterations and tolerance in each time step
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for multiphysics problems, use a tolerance small enough to converge the results
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 1.0, 100, 1, 30, 1e-3, 0.5)?;

    Ok(())
}
