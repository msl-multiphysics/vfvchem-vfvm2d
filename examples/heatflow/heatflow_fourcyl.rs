use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

const SPACE: f64 = 0.001;

// cargo run --release --example heatflow_fourcyl
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatflow_fourcyl").unwrap();

    // create variables
    let mut vars = Variables::new_obj("examples/obj/obj_fourcyl.obj", SPACE, 100, 42)?;

    // geometry
    // the heated cylinders share one domain, one interface, and one no-slip boundary
    let dom_c = vars.add_dom("dom_c")?;
    let dom_h = vars.add_dom("dom_h")?;
    let bnd_i = vars.add_bnd(dom_c, "bnd_i")?;
    let bnd_o = vars.add_bnd(dom_c, "bnd_o")?;
    let bnd_w = vars.add_bnd(dom_c, "bnd_w")?;
    let bnd_h = vars.add_bnd(dom_c, "bnd_h")?;
    let itf_h = vars.add_itf(dom_c, dom_h, "itf_h")?;

    // unknown domain values
    let temp_c = vars.add_scldom_unkn(dom_c, 300.0, "examples/output_heatflow_fourcyl/temp_c.vtu")?;
    let temp_h = vars.add_scldom_unkn(dom_h, 300.0, "examples/output_heatflow_fourcyl/temp_h.vtu")?;
    let vel = vars.add_vecdom_unkn(dom_c, (1.0, 0.0), "examples/output_heatflow_fourcyl/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom_c, 0.0, "examples/output_heatflow_fourcyl/pres.vtu")?;

    // constant domain values
    // fluid in the channel, heated cylinders
    let vlcp_c = vars.add_scldom_cnst(dom_c, 1000.0, "")?;
    let cond_c = vars.add_scldom_cnst(dom_c, 0.026, "")?;
    let hsrc_c = vars.add_scldom_cnst(dom_c, 0.0, "")?;
    let vlcp_h = vars.add_scldom_cnst(dom_h, 550.0, "")?;
    let cond_h = vars.add_scldom_cnst(dom_h, 50.0, "")?;
    let hsrc_h = vars.add_scldom_cnst(dom_h, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom_c, 1.2, "")?;
    let visc = vars.add_scldom_cnst(dom_c, 1.8e-5, "")?;
    let fce = vars.add_vecdom_cnst(dom_c, (0.0, 0.0), "")?;

    // parabolic inlet on y in [0.0, 0.05], average 0.2, maximum 0.3
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 0.3 * y * (0.05 - y) / (0.05 * 0.05);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;

    // constant boundary values
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0, "")?;
    let hflx_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0), "")?;
    let vel_h = vars.add_vecbnd_cnst(bnd_h, (0.0, 0.0), "")?;

    // steady-state heat and momentum transfer solver
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom_c, temp_c, vlcp_c, cond_c, hsrc_c)?;
    phys.add_heat_dom(dom_h, temp_h, vlcp_h, cond_h, hsrc_h)?;
    phys.add_temp_bnd(bnd_i, temp_i)?;
    phys.add_hout_bnd(bnd_o)?;
    phys.add_hflx_bnd(bnd_w, hflx_w)?;
    phys.add_hcnt_itf(itf_h)?;
    phys.add_flow_dom(dom_c, vel, pres, den, visc, fce)?;
    phys.add_vel_bnd(bnd_i, vel_i)?;
    phys.add_pres_bnd(bnd_o, pres_o)?;
    phys.add_vel_bnd(bnd_w, vel_w)?;
    phys.add_vel_bnd(bnd_h, vel_h)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 0.5)?;

    Ok(())
}
