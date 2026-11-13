use vfvchem_vfvm2d::*;

/// Steady-state heat transfer in a square domain.
/// Run with: `cargo run --release --example quickstart`
///
/// Geometry
/// --------
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
/// dom - main domain
///
/// Boundaries:
/// bnd_l - left boundary
/// bnd_r - right boundary
/// bnd_b - bottom boundary
/// bnd_t - top boundary
///
/// Heat Transfer
/// -------------
/// 
/// Domain Properties:
/// dom - thermal conductivity (k = 0.1 W m-1 K-1)
/// dom - heat source (Q = -500.0 W m-3)
///
/// Boundary Conditions:
/// bnd_l - outward flux (q = 100.0 W m-2)
/// bnd_r - temperature (T = 300 K)
/// bnd_b - no flux (q = 0.0 W m-2)
/// bnd_t - temperature (T = 400 K)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    // create rectangular geometry with Variables::new_bounds
    // nds_space - spacing at the boundaries
    // pts_space - spacing far from the boundaries
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.1, 0.1, 0.002, 0.005, "examples/output_quickstart")?;

    // geometry
    // define domain and boundaries with vars.add_dom and vars.add_bnd
    // for new_bounds, the domain is "dom" and the boundaries are "bnd_l", "bnd_r", "bnd_b", and "bnd_t"
    let dom = vars.add_dom("dom")?;  // domain
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top boundary

    // unknowns
    // define unknowns vars.add_scldom_unkn
    let temp = vars.add_scldom_unkn(dom, 0.0, "temp")?;  // temperature (T); put initial guess for temperature

    // domain values
    // define domain properties with vars.add_scldom_cnst
    let cond = vars.add_scldom_cnst(dom, 0.1, "cond")?;     // thermal conductivity (k)
    let hlin = vars.add_scldom_cnst(dom, 0.0, "hlin")?;     // linear source (L); multiply by temperature to yield additional source
    let hsrc = vars.add_scldom_cnst(dom, -500.0, "hsrc")?;  // heat source (Q); positive if source, negative if sink

    // boundary values
    // define boundary properties with vars.add_sclbnd_cnst
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 100.0)?;  // heat flux (q); positive if outward, negative if inward
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0)?;  // temperature
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;    // heat flux
    let temp_t = vars.add_sclbnd_cnst(bnd_t, 400.0)?;  // temperature

    // physics
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom, temp, cond, hlin, hsrc)?;    // heat transfer domain: dom, T, k, L, Q
    phys.add_hflx_bnd(bnd_l, hflx_l)?;  // heat flux boundary: bnd, q
    phys.add_temp_bnd(bnd_r, temp_r)?;  // temperature boundary: bnd, T
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // heat flux boundary: bnd, q
    phys.add_temp_bnd(bnd_t, temp_t)?;  // temperature boundary: bnd, T

    // solver
    // use higher damping (0.6 ~ 0.9) for speed; lower (0.1 ~ 0.3) for stability
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 10, 1e-3, 1.0)?;  // physics solver

    Ok(())
}
