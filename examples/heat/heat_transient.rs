use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

const SPACE: f64 = 0.05;

// cargo run --release --example heat_transient
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heat_transient").unwrap();

    // create variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, SPACE, SPACE, SPACE, 100, 42)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknown domain values
    let temp = vars.add_scldom_unkn(dom, 0.0, "examples/output_heat_transient/temp.vtu")?;

    // constant domain values
    let vlcp = vars.add_scldom_cnst(dom, 10.0, "")?;
    let cond = vars.add_scldom_cnst(dom, 1.0, "")?;
    let hsrc = vars.add_scldom_cnst(dom, -500.0, "")?;

    // constant boundary values
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 50.0, "")?;
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0, "")?;
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let temp_t = vars.add_sclbnd_cnst(bnd_t, 400.0, "")?;

    // transient heat transfer solver
    let mut phys = TransientHeat::new()?;
    phys.add_heat_dom(dom, temp, vlcp, cond, hsrc)?;
    phys.add_hflx_bnd(bnd_l, hflx_l)?;
    phys.add_temp_bnd(bnd_r, temp_r)?;
    phys.add_hflx_bnd(bnd_b, hflx_b)?;
    phys.add_temp_bnd(bnd_t, temp_t)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 0.1, 100, 1, 10, 1e-3, 1.0)?;

    Ok(())
}
