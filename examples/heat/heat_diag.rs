use vfvchem_vfvm2d::*;

/// Steady-state heat transfer in a divided domain.
/// Run with: `cargo run --release --example heat_diag`
///
/// Geometry
/// --------
/// 
///          bnd_t
///       +---------+
///       |dom_tl / |
/// bnd_l |    /    | bnd_r
///       | / dom_br|
///       +---------+
///          bnd_b
///
/// Domains:
/// dom_br - bottom-right triangle
/// dom_tl - top-left triangle
///
/// Boundaries:
/// bnd_l - left boundary
/// bnd_r - right boundary
/// bnd_b - bottom boundary
/// bnd_t - top boundary
///
/// Interfaces:
/// itf_d - diagonal interface
/// 
/// Heat Transfer
/// -------------
/// 
/// Domain Properties:
/// dom_br - thermal conductivity (k = 0.1 W m-1 K-1)
/// dom_br - heat source (Q = -500.0 W m-3)
/// dom_tl - thermal conductivity (k = 0.1 W m-1 K-1)
/// dom_tl - heat source (Q = -500.0 W m-3)
///
/// Boundary Conditions:
/// bnd_l - outward flux (q = 100.0 W m-2)
/// bnd_r - temperature (T = 300 K)
/// bnd_b - no flux (q = 0.0 W m-2)
/// bnd_t - temperature (T = 400 K)
/// 
/// Interface Conditions:
/// itf_d - temperature and flux continuity
///
fn main() -> Result<(), VFVChemError> {
    // variables
    // import custom geometry with Variables::new_obj
    let mut vars = Variables::new_obj("examples/obj/obj_diag.obj", 0.002, 0.005, "examples/output_heat_diag")?;

    // geometry
    // for new_obj, the domain, boundary, and interface names are defined in the obj file
    let dom_br = vars.add_dom("dom_br")?;  // bottom-right triangle
    let dom_tl = vars.add_dom("dom_tl")?;  // top-left triangle
    let bnd_l = vars.add_bnd(dom_tl, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom_br, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom_br, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom_tl, "bnd_t")?;  // top boundary
    let itf_d = vars.add_itf(dom_br, dom_tl, "itf_d")?;  // diagonal interface

    // unknowns
    let temp_br = vars.add_scldom_unkn(dom_br, 0.0, "temp")?;  // temperature (T)
    let temp_tl = vars.add_scldom_unkn(dom_tl, 0.0, "temp")?;  // temperature

    // domain values
    let cond_br = vars.add_scldom_cnst(dom_br, 0.1, "cond")?;     // thermal conductivity (k)
    let hlin_br = vars.add_scldom_cnst(dom_br, 0.0, "hlin")?;     // linear source (L); multiply by temperature to yield additional source
    let hsrc_br = vars.add_scldom_cnst(dom_br, -500.0, "hsrc")?;  // heat source (Q); positive if source, negative if sink
    let cond_tl = vars.add_scldom_cnst(dom_tl, 0.1, "cond")?;     // thermal conductivity
    let hlin_tl = vars.add_scldom_cnst(dom_tl, 0.0, "hlin")?;     // linear source
    let hsrc_tl = vars.add_scldom_cnst(dom_tl, -500.0, "hsrc")?;  // heat source

    // boundary values
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 100.0)?;  // heat flux (q); positive if outward, negative if inward
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0)?;  // temperature
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;    // heat flux
    let temp_t = vars.add_sclbnd_cnst(bnd_t, 400.0)?;  // temperature

    // physics
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom_br, temp_br, cond_br, hlin_br, hsrc_br)?;  // heat transfer domain: dom, T, k, L, Q
    phys.add_heat_dom(dom_tl, temp_tl, cond_tl, hlin_tl, hsrc_tl)?;  // heat transfer domain: dom, T, k, L, Q
    phys.add_hflx_bnd(bnd_l, hflx_l)?;  // heat flux boundary: bnd, q
    phys.add_temp_bnd(bnd_r, temp_r)?;  // temperature boundary: bnd, T
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // heat flux boundary: bnd, q
    phys.add_temp_bnd(bnd_t, temp_t)?;  // temperature boundary: bnd, T
    phys.add_hcnt_itf(itf_d)?;  // continuity interface: itf

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 10, 1e-3, 1.0)?;  // physics solver

    Ok(())
}
