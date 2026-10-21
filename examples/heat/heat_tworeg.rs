use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

const SPACE: f64 = 0.05;

// cargo run --release --example heat_tworeg
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heat_tworeg").unwrap();

    // create variables
    let mut vars = Variables::new_obj("examples/obj/obj_tworeg.obj", SPACE, 100, 42)?;

    // geometry
    let dom_lo = vars.add_dom("dom_lo")?;
    let dom_up = vars.add_dom("dom_up")?;
    let bnd_b = vars.add_bnd(dom_lo, "bnd_b")?;
    let bnd_r = vars.add_bnd(dom_lo, "bnd_r")?;
    let bnd_t = vars.add_bnd(dom_up, "bnd_t")?;
    let bnd_l = vars.add_bnd(dom_up, "bnd_l")?;
    let itf = vars.add_itf(dom_lo, dom_up, "itf")?;

    // unknown domain values
    let temp_lo = vars.add_scldom_unkn(dom_lo, 0.0, "examples/output_heat_tworeg/temp_lo.vtu")?;
    let temp_up = vars.add_scldom_unkn(dom_up, 0.0, "examples/output_heat_tworeg/temp_up.vtu")?;

    // constant domain values
    let cond_lo = vars.add_scldom_cnst(dom_lo, 1.0, "")?;
    let cond_up = vars.add_scldom_cnst(dom_up, 1.0, "")?;
    let hsrc_lo = vars.add_scldom_cnst(dom_lo, -500.0, "")?;
    let hsrc_up = vars.add_scldom_cnst(dom_up, -500.0, "")?;

    // constant boundary values
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 50.0, "")?;
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0, "")?;
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let temp_t = vars.add_sclbnd_cnst(bnd_t, 400.0, "")?;

    // steady-state heat transfer solver
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom_lo, temp_lo, cond_lo, hsrc_lo)?;
    phys.add_heat_dom(dom_up, temp_up, cond_up, hsrc_up)?;
    phys.add_hflx_bnd(bnd_l, hflx_l)?;
    phys.add_temp_bnd(bnd_r, temp_r)?;
    phys.add_hflx_bnd(bnd_b, hflx_b)?;
    phys.add_temp_bnd(bnd_t, temp_t)?;
    phys.add_hcnt_itf(itf)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-3, 1.0)?;

    Ok(())
}
