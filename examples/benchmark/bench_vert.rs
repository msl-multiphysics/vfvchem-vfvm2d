use vfvchem_vfvm2d::*;

/// Composite-wall heat conduction with discontinuous conductivity.
/// Run with: `cargo run --release --example bench_vert`
///
/// Geometry
/// --------
/// 
///        bnd_tl   bnd_tr
///       +-------+-------+
///       |       |       |
/// bnd_l | dom_l | dom_r | bnd_r
///       |       |       |
///       +-------+-------+
///        bnd_bl   bnd_br
///
/// Domains:
/// dom_l - left block
/// dom_r - right block
///
/// Boundaries:
/// bnd_l - left boundary
/// bnd_r - right boundary
/// bnd_bl - bottom boundary of left block
/// bnd_tl - top boundary of left block
/// bnd_br - bottom boundary of right block
/// bnd_tr - top boundary of right block
///
/// Interfaces:
/// itf_v - vertical interface
///
/// Heat Transfer
/// -------------
/// 
/// Domain Properties:
/// dom_l - thermal conductivity (k = 1.0 W m-1 K-1)
/// dom_r - thermal conductivity (k = 1000.0 W m-1 K-1)
///
/// Boundary Conditions:
/// bnd_l - temperature (T = 400 K)
/// bnd_r - temperature (T = 300 K)
/// bnd_bl - no flux (q = 0.0 W m-2)
/// bnd_tl - no flux (q = 0.0 W m-2)
/// bnd_br - no flux (q = 0.0 W m-2)
/// bnd_tr - no flux (q = 0.0 W m-2)
///
/// Interface Conditions:
/// itf_v - continuous temperature and heat flux
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_obj("examples/obj/obj_vert.obj", 0.005, 0.005, "examples/output_bench_vert")?;

    // geometry
    let dom_l = vars.add_dom("dom_l")?;  // left block
    let dom_r = vars.add_dom("dom_r")?;  // right block
    let bnd_l = vars.add_bnd(dom_l, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom_r, "bnd_r")?;  // right boundary
    let bnd_bl = vars.add_bnd(dom_l, "bnd_bl")?;  // bottom boundary of left block
    let bnd_tl = vars.add_bnd(dom_l, "bnd_tl")?;  // top boundary of left block
    let bnd_br = vars.add_bnd(dom_r, "bnd_br")?;  // bottom boundary of right block
    let bnd_tr = vars.add_bnd(dom_r, "bnd_tr")?;  // top boundary of right block
    let itf_v = vars.add_itf(dom_l, dom_r, "itf_v")?;  // vertical interface

    // unknowns
    let temp_l = vars.add_scldom_unkn(dom_l, 350.0, "temp")?;  // temperature (T)
    let temp_r = vars.add_scldom_unkn(dom_r, 350.0, "temp")?;  // temperature (T)

    // constant domain values
    let cond_l = vars.add_scldom_cnst(dom_l, 1.0, "cond")?;  // thermal conductivity (k)
    let cond_r = vars.add_scldom_cnst(dom_r, 1000.0, "cond")?;  // thermal conductivity (k)
    let hlin_l = vars.add_scldom_cnst(dom_l, 0.0, "hlin")?;  // linear source (L)
    let hsrc_l = vars.add_scldom_cnst(dom_l, 0.0, "hsrc")?;  // heat source (Q)
    let hlin_r = vars.add_scldom_cnst(dom_r, 0.0, "hlin")?;  // linear source (L)
    let hsrc_r = vars.add_scldom_cnst(dom_r, 0.0, "hsrc")?;  // heat source (Q)

    // constant boundary values
    let tset_l = vars.add_sclbnd_cnst(bnd_l, 400.0)?;  // temperature
    let tset_r = vars.add_sclbnd_cnst(bnd_r, 300.0)?;  // temperature
    let hflx_bl = vars.add_sclbnd_cnst(bnd_bl, 0.0)?;  // heat flux (q)
    let hflx_tl = vars.add_sclbnd_cnst(bnd_tl, 0.0)?;  // heat flux
    let hflx_br = vars.add_sclbnd_cnst(bnd_br, 0.0)?;  // heat flux
    let hflx_tr = vars.add_sclbnd_cnst(bnd_tr, 0.0)?;  // heat flux

    // physics
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom_l, temp_l, cond_l, hlin_l, hsrc_l)?;  // heat transfer domain: dom, T, k, L, Q
    phys.add_heat_dom(dom_r, temp_r, cond_r, hlin_r, hsrc_r)?;  // heat transfer domain: dom, T, k, L, Q
    phys.add_temp_bnd(bnd_l, tset_l)?;  // temperature boundary: bnd, T
    phys.add_temp_bnd(bnd_r, tset_r)?;  // temperature boundary: bnd, T
    phys.add_hflx_bnd(bnd_bl, hflx_bl)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_tl, hflx_tl)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_br, hflx_br)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_tr, hflx_tr)?;  // heat flux boundary: bnd, q
    phys.add_hcnt_itf(itf_v)?;  // continuity interface: itf

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 5, 1e-4, 1.0)?;  // physics solver

    Ok(())
}
