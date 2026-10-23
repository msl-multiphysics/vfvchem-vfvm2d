use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

/// Steady-state heat equation in two regions joined by an interface.
/// Run with: `cargo run --release --example heat_tworeg`
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
/// dom_lo - thermal conductivity (k = 1.0 W m-1 K-1)
/// dom_lo - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_lo - heat source (Q = -500.0 W m-3)
/// dom_up - thermal conductivity (k = 1.0 W m-1 K-1)
/// dom_up - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_up - heat source (Q = -500.0 W m-3)
///
/// Boundary conditions:
/// bnd_l - outward flux (q = 50.0 W m-2)
/// bnd_r - temperature (T = 300 K)
/// bnd_b - no flux (q = 0.0 W m-2)
/// bnd_t - temperature (T = 400 K)
///
/// Interface conditions:
/// itf_d - temperature and flux continuity
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heat_tworeg").unwrap();

    // nodes, points, and variables
    // new_obj - import boundary, interface, and domain polylines from an obj file
    // arguments: input_file, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_obj("examples/obj/obj_tworeg.obj", 0.05, 0.05, 0.05, 100, 42)?;

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
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., temp_lo_0000.vtu)
    let temp_lo = vars.add_scldom_unkn(dom_lo, 0.0, "examples/output_heat_tworeg/temp_lo.vtu")?;
    let temp_up = vars.add_scldom_unkn(dom_up, 0.0, "examples/output_heat_tworeg/temp_up.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    // heat production is multiplied by the temperature; heat source adds the source directly
    // use heat source by default; use heat production (or both) if heat source does not converge
    let cond_lo = vars.add_scldom_cnst(dom_lo, 1.0, "")?;  // thermal conductivity
    let cond_up = vars.add_scldom_cnst(dom_up, 1.0, "")?;  // thermal conductivity
    let hprd_lo = vars.add_scldom_cnst(dom_lo, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hprd_up = vars.add_scldom_cnst(dom_up, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc_lo = vars.add_scldom_cnst(dom_lo, -500.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let hsrc_up = vars.add_scldom_cnst(dom_up, -500.0, "")?;  // heat source (Q; positive if source; negative if sink)

    // constant boundary values
    // arguments: boundary, value, output_file
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 50.0, "")?;  // heat flux (positive if outward; negative if inward)
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0, "")?;  // temperature
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;  // heat flux
    let temp_t = vars.add_sclbnd_cnst(bnd_t, 400.0, "")?;  // temperature

    // steady-state heat transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_hflx_bnd - register boundary with heat flux
    // add_temp_bnd - register boundary with temperature
    // add_hcnt_itf - register interface with temperature and flux continuity
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom_lo, temp_lo, cond_lo, hprd_lo, hsrc_lo)?;  // arguments: domain, T, k, s, Q
    phys.add_heat_dom(dom_up, temp_up, cond_up, hprd_up, hsrc_up)?;  // arguments: domain, T, k, s, Q
    phys.add_hflx_bnd(bnd_l, hflx_l)?;  // arguments: boundary, q
    phys.add_temp_bnd(bnd_r, temp_r)?;  // arguments: boundary, T
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // arguments: boundary, q
    phys.add_temp_bnd(bnd_t, temp_t)?;  // arguments: boundary, T
    phys.add_hcnt_itf(itf_d)?;  // arguments: interface

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-3, 1.0)?;

    Ok(())
}
