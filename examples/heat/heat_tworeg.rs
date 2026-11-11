use vfvchem_vfvm2d::*;

// cargo run --release --example heat_tworeg
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_obj("examples/obj/obj_tworeg.obj", 0.02, 0.05, "examples/output_heat_tworeg")?;

    // geometry
    let dom_lo = vars.add_dom("dom_lo")?;
    let dom_up = vars.add_dom("dom_up")?;
    let bnd_b = vars.add_bnd(dom_lo, "bnd_b")?;
    let bnd_r = vars.add_bnd(dom_lo, "bnd_r")?;
    let bnd_t = vars.add_bnd(dom_up, "bnd_t")?;
    let bnd_l = vars.add_bnd(dom_up, "bnd_l")?;
    let itf_d = vars.add_itf(dom_lo, dom_up, "itf_d")?;

    // unknowns
    let temp_lo = vars.add_scldom_unkn(dom_lo, 0.0, "temp")?;
    let temp_up = vars.add_scldom_unkn(dom_up, 0.0, "temp")?;

    // domain values
    let cond_lo = vars.add_scldom_cnst(dom_lo, 1.0, "cond")?;
    let cond_up = vars.add_scldom_cnst(dom_up, 1.0, "cond")?;
    let hlin_lo = vars.add_scldom_cnst(dom_lo, 0.0, "hlin")?;
    let hlin_up = vars.add_scldom_cnst(dom_up, 0.0, "hlin")?;
    let hsrc_lo = vars.add_scldom_cnst(dom_lo, -500.0, "hsrc")?;  // positive if source; negative if sink
    let hsrc_up = vars.add_scldom_cnst(dom_up, -500.0, "hsrc")?;

    // boundary values
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // positive if outward; negative if inward
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0)?;
    let temp_t = vars.add_sclbnd_cnst(bnd_t, 400.0)?;
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 50.0)?;

    // physics
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom_lo, temp_lo, cond_lo, hlin_lo, hsrc_lo)?;  // arguments: domain, T, k, Q
    phys.add_heat_dom(dom_up, temp_up, cond_up, hlin_up, hsrc_up)?;  // arguments: domain, T, k, Q
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // arguments: boundary, q
    phys.add_temp_bnd(bnd_r, temp_r)?;  // arguments: boundary, T
    phys.add_temp_bnd(bnd_t, temp_t)?;  // arguments: boundary, T
    phys.add_hflx_bnd(bnd_l, hflx_l)?;  // arguments: boundary, q
    phys.add_hcnt_itf(itf_d)?;  // arguments: interface

    // solver
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-3, 1.0)?;

    Ok(())
}
