use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

const SPACE: f64 = 0.02;

// cargo run --release --example heatflow_step_uniform
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatflow_step_uniform").unwrap();

    // create variables
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", SPACE, 100, 42)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_i = vars.add_bnd(dom, "bnd_i")?;
    let bnd_o = vars.add_bnd(dom, "bnd_o")?;
    let bnd_w = vars.add_bnd(dom, "bnd_w")?;

    // unknown domain values
    let temp = vars.add_scldom_unkn(dom, 300.0, "examples/output_heatflow_step_uniform/temp.vtu")?;
    let vel = vars.add_vecdom_unkn(dom, (1.0, 0.0), "examples/output_heatflow_step_uniform/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_heatflow_step_uniform/pres.vtu")?;

    // constant domain values
    let vlcp = vars.add_scldom_cnst(dom, 1.0, "")?;
    let cond = vars.add_scldom_cnst(dom, 0.01, "")?;
    let hsrc = vars.add_scldom_cnst(dom, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom, 1.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 0.01, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;

    // parabolic inlet on y in [0.2, 0.4], maximum 1.5
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 1.5 * (y - 0.2) * (0.4 - y) / (0.2 * 0.2);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;

    // constant boundary values
    // no heat source and inlet at the initial temperature, so the temperature stays uniform
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0, "")?;
    let hflx_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0), "")?;

    // transient heat and momentum transfer solver
    let mut phys = TransientHeatFlow::new()?;
    phys.add_heat_dom(dom, temp, vlcp, cond, hsrc)?;
    phys.add_temp_bnd(bnd_i, temp_i)?;
    phys.add_hout_bnd(bnd_o)?;
    phys.add_hflx_bnd(bnd_w, hflx_w)?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
    phys.add_vel_bnd(bnd_i, vel_i)?;
    phys.add_pres_bnd(bnd_o, pres_o)?;
    phys.add_vel_bnd(bnd_w, vel_w)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 0.1, 100, 1, 30, 1e-3, 0.5)?;

    Ok(())
}
