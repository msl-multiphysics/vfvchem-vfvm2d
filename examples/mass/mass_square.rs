use vfvchem_vfvm2d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Steady-state multicomponent diffusion equation with cross-diffusion of component 0.
/// Run with: `cargo run --release --example mass_square`
///
/// Geometry:
/// Square domain defined by lengths along each axis.
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
/// dom - main domain (1.0 m x 1.0 m)
///
/// Boundaries:
/// bnd_l - left boundary (1.0 m)
/// bnd_r - right boundary (1.0 m)
/// bnd_b - bottom boundary (1.0 m)
/// bnd_t - top boundary (1.0 m)
///
/// Properties:
/// dom - diffusion coefficients (m2 s-1)
/// [D_00 D_01] = [1.0e-4 0.5e-4]
/// [D_10 D_11]   [0.0    1.0e-4]
/// dom - reaction rates (mol m-3 s-1; positive if production)
/// [R_0] = [-0.2]
/// [R_1]   [+0.1]
///
/// Boundary conditions:
/// bnd_l - outward flux 0; concentration 1 (mol m-2 s-1; mol m-3)
/// [N_0] = [0.5]
/// [c_1]   [1.0e-4]
/// bnd_r - concentration (mol m-3)
/// [c_0] = [3.0e-4]
/// [c_1]   [0.0]
/// bnd_b - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
/// bnd_t - concentration 0; no flux 1 (mol m-3; mol m-2 s-1)
/// [c_0] = [4.0e-4]
/// [N_1]   [0.0]
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_mass_square").unwrap();

    // nodes, points, and variables
    // new_bounds - create boundary nodes and internal points inside a bounding box
    // arguments: x_min, y_min, x_max, y_max, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.05, 100, 42)?;

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
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., conc_0_0000.vtu)
    let conc_0 = vars.add_scldom_unkn(dom, 0.0, "examples/output_mass_square/conc_0.vtu")?;
    let conc_1 = vars.add_scldom_unkn(dom, 0.0, "examples/output_mass_square/conc_1.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    // no need to declare zero diffusion coefficients
    // component 0 diffuses down its own gradient and the gradient of component 1
    let diff_00 = vars.add_scldom_cnst(dom, 1.0e-4, "")?;  // D_00
    let diff_01 = vars.add_scldom_cnst(dom, 0.5e-4, "")?;  // D_01
    let diff_11 = vars.add_scldom_cnst(dom, 1.0e-4, "")?;  // D_11
    let msrc_0 = vars.add_scldom_cnst(dom, -0.2, "")?;  // mass source (R_0; positive if source; negative if sink)
    let msrc_1 = vars.add_scldom_cnst(dom,  0.1, "")?;  // mass source (R_1; positive if source; negative if sink)

    // constant boundary values
    // arguments: boundary, value, output_file
    let mflx_0_l = vars.add_sclbnd_cnst(bnd_l, 0.5, "")?;  // molar flux (positive if outward; negative if inward)
    let conc_0_r = vars.add_sclbnd_cnst(bnd_r, 3e-4, "")?;  // concentration
    let mflx_0_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;  // molar flux
    let conc_0_t = vars.add_sclbnd_cnst(bnd_t, 4e-4, "")?;  // concentration
    let conc_1_l = vars.add_sclbnd_cnst(bnd_l, 1e-4, "")?;  // concentration
    let conc_1_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;  // concentration
    let mflx_1_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;  // molar flux
    let mflx_1_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;  // molar flux

    // steady-state mass transfer solver
    // SteadyMass::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    // diffusion coefficients must be collected in a hashmap with driving concentration index as key
    let mut phys = SteadyMass::new(2)?;
    phys.add_mass_dom(0, dom, conc_0, HashMap::from([(0, diff_00), (1, diff_01)]), msrc_0)?;  // arguments: component, domain, c, D, R
    phys.add_mass_dom(1, dom, conc_1, HashMap::from([(1, diff_11)]), msrc_1)?;  // arguments: component, domain, c, D, R
    phys.add_mflx_bnd(0, bnd_l, mflx_0_l)?;  // arguments: component, boundary, N
    phys.add_conc_bnd(0, bnd_r, conc_0_r)?;  // arguments: component, boundary, c
    phys.add_mflx_bnd(0, bnd_b, mflx_0_b)?;  // arguments: component, boundary, N
    phys.add_conc_bnd(0, bnd_t, conc_0_t)?;  // arguments: component, boundary, c
    phys.add_conc_bnd(1, bnd_l, conc_1_l)?;  // arguments: component, boundary, c
    phys.add_conc_bnd(1, bnd_r, conc_1_r)?;  // arguments: component, boundary, c
    phys.add_mflx_bnd(1, bnd_b, mflx_1_b)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(1, bnd_t, mflx_1_t)?;  // arguments: component, boundary, N

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;

    Ok(())
}
