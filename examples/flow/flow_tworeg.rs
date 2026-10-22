use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

const SPACE: f64 = 0.02;

// cargo run --release --example flow_tworeg
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_tworeg").unwrap();

    // create variables
    let mut vars = Variables::new_obj("examples/obj/obj_tworeg.obj", SPACE, SPACE, SPACE, 100, 42)?;

    // geometry
    let dom_lo = vars.add_dom("dom_lo")?;
    let dom_up = vars.add_dom("dom_up")?;
    let bnd_b = vars.add_bnd(dom_lo, "bnd_b")?;
    let bnd_r = vars.add_bnd(dom_lo, "bnd_r")?;
    let bnd_t = vars.add_bnd(dom_up, "bnd_t")?;
    let bnd_l = vars.add_bnd(dom_up, "bnd_l")?;
    let itf = vars.add_itf(dom_lo, dom_up, "itf")?;

    // unknown domain values
    let vel_lo = vars.add_vecdom_unkn(dom_lo, (0.0, 0.0), "examples/output_flow_tworeg/vel_lo.vtu")?;
    let vel_up = vars.add_vecdom_unkn(dom_up, (0.0, 0.0), "examples/output_flow_tworeg/vel_up.vtu")?;
    let pres_lo = vars.add_scldom_unkn(dom_lo, 0.0, "examples/output_flow_tworeg/pres_lo.vtu")?;
    let pres_up = vars.add_scldom_unkn(dom_up, 0.0, "examples/output_flow_tworeg/pres_up.vtu")?;

    // constant domain values
    let den_lo = vars.add_scldom_cnst(dom_lo, 1000.0, "")?;
    let den_up = vars.add_scldom_cnst(dom_up, 1000.0, "")?;
    let visc_lo = vars.add_scldom_cnst(dom_lo, 0.001, "")?;
    let visc_up = vars.add_scldom_cnst(dom_up, 0.001, "")?;
    let fce_lo = vars.add_vecdom_cnst(dom_lo, (0.0, 0.0), "")?;
    let fce_up = vars.add_vecdom_cnst(dom_up, (0.0, 0.0), "")?;

    // constant boundary values
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0), "")?;
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0), "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (4e-4, 0.0), "")?;

    // steady-state momentum transfer solver
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom_lo, vel_lo, pres_lo, den_lo, visc_lo, fce_lo)?;
    phys.add_flow_dom(dom_up, vel_up, pres_up, den_up, visc_up, fce_up)?;
    phys.add_vel_bnd(bnd_l, vel_l)?;
    phys.add_vel_bnd(bnd_r, vel_r)?;
    phys.add_vel_bnd(bnd_b, vel_b)?;
    phys.add_vel_bnd(bnd_t, vel_t)?;
    phys.add_fcnt_itf(itf)?;
    phys.set_pres_ref(dom_lo, 0, 0.0)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;

    Ok(())
}
