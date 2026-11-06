use vfvchem_vfvm2d::*;

// cargo run --release --example quickstart
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(
        0.0, 0.0, 1.0, 1.0, 0.02, 0.05, 0.10,
        100, 42, "examples/output_quickstart"
    )?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let temp = vars.add_scldom_unkn(dom, 0.0, "temp")?;

    // domain values
    let cond = vars.add_scldom_cnst(dom, 1.0, "cond")?;
    let hlin = vars.add_scldom_cnst(dom, 0.0, "hlin")?;
    let hsrc = vars.add_scldom_cnst(dom, -500.0, "hsrc")?;  // positive if source; negative if sink

    // boundary values
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 50.0)?;  // positive if outward; negative if inward
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0)?;
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;
    let temp_t = vars.add_sclbnd_cnst(bnd_t, 400.0)?;

    // physics
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom, temp, cond, hlin, hsrc)?;  // arguments: domain, T, k, Q
    phys.add_hflx_bnd(bnd_l, hflx_l)?;  // arguments: boundary, q
    phys.add_temp_bnd(bnd_r, temp_r)?;  // arguments: boundary, T
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // arguments: boundary, q
    phys.add_temp_bnd(bnd_t, temp_t)?;  // arguments: boundary, T

    // solver
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-3, 1.0)?;

    Ok(())
}
