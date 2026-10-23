use vfvchem_vfvm2d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Heat and mass transfer through two regions with different properties joined by an interface.
/// Run with: `cargo run --release --example feat_itf`
///
/// The interface conditions make the unknown and its flux continuous across the interface.
/// Heat: water below and glass above, joined with add_hcnt_itf.
/// Mass: two diffusion coefficients, joined with add_mcnt_itf.
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
/// Heat Transfer
///
/// Properties:
/// dom_lo - thermal conductivity (k = 0.6 W m-1 K-1)
/// dom_lo - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_lo - heat source (Q = 0.0 W m-3)
/// dom_up - thermal conductivity (k = 1.4 W m-1 K-1)
/// dom_up - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_up - heat source (Q = 0.0 W m-3)
///
/// Boundary conditions:
/// bnd_l - temperature (T = 350 K)
/// bnd_r - temperature (T = 300 K)
/// bnd_b - no flux (q = 0.0 W m-2)
/// bnd_t - no flux (q = 0.0 W m-2)
///
/// Interface conditions:
/// itf_d - temperature and flux continuity
///
/// Mass Transfer
///
/// Properties:
/// dom_lo - diffusion coefficient (D_00 = 1.0e-9 m2 s-1)
/// dom_lo - production coefficient (s_0 = 0.0 s-1)
/// dom_lo - reaction rate (R_0 = 0.0 mol m-3 s-1)
/// dom_up - diffusion coefficient (D_00 = 4.0e-9 m2 s-1)
/// dom_up - production coefficient (s_0 = 0.0 s-1)
/// dom_up - reaction rate (R_0 = 0.0 mol m-3 s-1)
///
/// Boundary conditions:
/// bnd_l - concentration (c_0 = 1.0 mol m-3)
/// bnd_r - concentration (c_0 = 0.0 mol m-3)
/// bnd_b - no flux (N_0 = 0.0 mol m-2 s-1)
/// bnd_t - no flux (N_0 = 0.0 mol m-2 s-1)
///
/// Interface conditions:
/// itf_d - concentration and flux continuity
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_feat_itf").unwrap();

    // heat transfer

    // nodes, points, and variables
    // new_obj - import boundary, interface, and domain polylines from an obj file
    // arguments: input_file, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_obj("examples/obj/obj_tworeg.obj", 0.02, 0.04, 0.2, 100, 42)?;

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
    let temp_lo = vars.add_scldom_unkn(dom_lo, 300.0, "examples/output_feat_itf/temp_lo.vtu")?;
    let temp_up = vars.add_scldom_unkn(dom_up, 300.0, "examples/output_feat_itf/temp_up.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    // heat production is multiplied by the temperature; heat source adds the source directly
    // use heat source by default; use heat production (or both) if heat source does not converge
    let cond_lo = vars.add_scldom_cnst(dom_lo, 0.6, "")?;  // thermal conductivity (water)
    let cond_up = vars.add_scldom_cnst(dom_up, 1.4, "")?;  // thermal conductivity (glass)
    let hprd_lo = vars.add_scldom_cnst(dom_lo, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hprd_up = vars.add_scldom_cnst(dom_up, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc_lo = vars.add_scldom_cnst(dom_lo, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let hsrc_up = vars.add_scldom_cnst(dom_up, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)

    // constant boundary values
    // arguments: boundary, value, output_file
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 350.0, "")?;  // temperature
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0, "")?;  // temperature
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;  // heat flux (positive if outward; negative if inward)
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;  // heat flux

    // steady-state heat transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    // add_hflx_bnd - register boundary with heat flux
    // add_hcnt_itf - register interface with temperature and flux continuity
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom_lo, temp_lo, cond_lo, hprd_lo, hsrc_lo)?;  // arguments: domain, T, k, s, Q
    phys.add_heat_dom(dom_up, temp_up, cond_up, hprd_up, hsrc_up)?;  // arguments: domain, T, k, s, Q
    phys.add_temp_bnd(bnd_l, temp_l)?;  // arguments: boundary, T
    phys.add_temp_bnd(bnd_r, temp_r)?;  // arguments: boundary, T
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // arguments: boundary, q
    phys.add_hflx_bnd(bnd_t, hflx_t)?;  // arguments: boundary, q
    phys.add_hcnt_itf(itf_d)?;  // arguments: interface

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-8, 1.0)?;

    // mass transfer
    // each problem uses its own variables

    // nodes, points, and variables
    let mut vars = Variables::new_obj("examples/obj/obj_tworeg.obj", 0.02, 0.04, 0.2, 100, 42)?;

    // geometry
    let dom_lo = vars.add_dom("dom_lo")?;  // lower right
    let dom_up = vars.add_dom("dom_up")?;  // upper left
    let bnd_b = vars.add_bnd(dom_lo, "bnd_b")?;  // bottom
    let bnd_r = vars.add_bnd(dom_lo, "bnd_r")?;  // right
    let bnd_t = vars.add_bnd(dom_up, "bnd_t")?;  // top
    let bnd_l = vars.add_bnd(dom_up, "bnd_l")?;  // left
    let itf_d = vars.add_itf(dom_lo, dom_up, "itf_d")?;  // diagonal

    // unknown domain values
    let conc_lo = vars.add_scldom_unkn(dom_lo, 0.0, "examples/output_feat_itf/conc_lo.vtu")?;
    let conc_up = vars.add_scldom_unkn(dom_up, 0.0, "examples/output_feat_itf/conc_up.vtu")?;

    // constant domain values
    // mass production is multiplied by the concentration; mass source adds the source directly
    // use mass source by default; use mass production (or both) if mass source does not converge
    let diff_lo = vars.add_scldom_cnst(dom_lo, 1.0e-9, "")?;  // D_00
    let diff_up = vars.add_scldom_cnst(dom_up, 4.0e-9, "")?;  // D_00
    let mprd_lo = vars.add_scldom_cnst(dom_lo, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let mprd_up = vars.add_scldom_cnst(dom_up, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let msrc_lo = vars.add_scldom_cnst(dom_lo, 0.0, "")?;  // mass source (R_0; positive if source; negative if sink)
    let msrc_up = vars.add_scldom_cnst(dom_up, 0.0, "")?;  // mass source (R_0; positive if source; negative if sink)

    // constant boundary values
    let conc_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;  // concentration
    let conc_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;  // concentration
    let mflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;  // molar flux (positive if outward; negative if inward)
    let mflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;  // molar flux

    // steady-state mass transfer solver
    // SteadyMass::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    // add_mcnt_itf - register interface with concentration and flux continuity
    // diffusion coefficients must be collected in a hashmap with driving concentration index as key
    let mut phys = SteadyMass::new(1)?;
    phys.add_mass_dom(0, dom_lo, conc_lo, HashMap::from([(0, diff_lo)]), mprd_lo, msrc_lo)?;  // arguments: component, domain, c, D, s, R
    phys.add_mass_dom(0, dom_up, conc_up, HashMap::from([(0, diff_up)]), mprd_up, msrc_up)?;  // arguments: component, domain, c, D, s, R
    phys.add_conc_bnd(0, bnd_l, conc_l)?;  // arguments: component, boundary, c
    phys.add_conc_bnd(0, bnd_r, conc_r)?;  // arguments: component, boundary, c
    phys.add_mflx_bnd(0, bnd_b, mflx_b)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_t, mflx_t)?;  // arguments: component, boundary, N
    phys.add_mcnt_itf(0, itf_d)?;  // arguments: component, interface

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-8, 1.0)?;

    Ok(())
}
