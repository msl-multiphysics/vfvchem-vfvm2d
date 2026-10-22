use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

const SPACE: f64 = 0.05;

// cargo run --release --example heatflow_heater
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatflow_heater").unwrap();

    // create variables
    let mut vars = Variables::new_obj("examples/obj/obj_heater.obj", SPACE, SPACE, SPACE, 100, 42)?;

    // geometry
    let dom_c = vars.add_dom("dom_c")?;
    let dom_b = vars.add_dom("dom_b")?;
    let dom_t = vars.add_dom("dom_t")?;
    let bnd_i = vars.add_bnd(dom_c, "bnd_i")?;
    let bnd_o = vars.add_bnd(dom_c, "bnd_o")?;
    let bnd_bc = vars.add_bnd(dom_c, "bnd_bc")?;
    let bnd_tc = vars.add_bnd(dom_c, "bnd_tc")?;
    let bnd_b = vars.add_bnd(dom_b, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom_t, "bnd_t")?;
    let itf_b = vars.add_itf(dom_c, dom_b, "itf_b")?;
    let itf_t = vars.add_itf(dom_c, dom_t, "itf_t")?;

    // unknown domain values
    let temp_c = vars.add_scldom_unkn(dom_c, 300.0, "examples/output_heatflow_heater/temp_c.vtu")?;
    let temp_b = vars.add_scldom_unkn(dom_b, 300.0, "examples/output_heatflow_heater/temp_b.vtu")?;
    let temp_t = vars.add_scldom_unkn(dom_t, 300.0, "examples/output_heatflow_heater/temp_t.vtu")?;
    let vel = vars.add_vecdom_unkn(dom_c, (1.0, 0.0), "examples/output_heatflow_heater/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom_c, 0.0, "examples/output_heatflow_heater/pres.vtu")?;

    // constant domain values
    // fluid in the channel, heated solids above and below
    let vlcp_c = vars.add_scldom_cnst(dom_c, 4.2e6, "")?;
    let cond_c = vars.add_scldom_cnst(dom_c, 0.6, "")?;
    let hsrc_c = vars.add_scldom_cnst(dom_c, 0.0, "")?;
    let vlcp_b = vars.add_scldom_cnst(dom_b, 550.0, "")?;
    let cond_b = vars.add_scldom_cnst(dom_b, 50.0, "")?;
    let hsrc_b = vars.add_scldom_cnst(dom_b, 10.0, "")?;
    let vlcp_t = vars.add_scldom_cnst(dom_t, 550.0, "")?;
    let cond_t = vars.add_scldom_cnst(dom_t, 50.0, "")?;
    let hsrc_t = vars.add_scldom_cnst(dom_t, 10.0, "")?;
    let den = vars.add_scldom_cnst(dom_c, 1000.0, "")?;
    let visc = vars.add_scldom_cnst(dom_c, 0.001, "")?;
    let fce = vars.add_vecdom_cnst(dom_c, (0.0, 0.0), "")?;

    // parabolic inlet on y in [0.0, 0.5], maximum 3e-4
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 3e-4 * y * (0.5 - y) / (0.5 * 0.5);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;

    // constant boundary values
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0, "")?;
    let temp_bb = vars.add_sclbnd_cnst(bnd_b, 400.0, "")?;
    let temp_tt = vars.add_sclbnd_cnst(bnd_t, 400.0, "")?;
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;
    let vel_bc = vars.add_vecbnd_cnst(bnd_bc, (0.0, 0.0), "")?;
    let vel_tc = vars.add_vecbnd_cnst(bnd_tc, (0.0, 0.0), "")?;

    // steady-state heat and momentum transfer solver
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom_c, temp_c, vlcp_c, cond_c, hsrc_c)?;
    phys.add_heat_dom(dom_b, temp_b, vlcp_b, cond_b, hsrc_b)?;
    phys.add_heat_dom(dom_t, temp_t, vlcp_t, cond_t, hsrc_t)?;
    phys.add_temp_bnd(bnd_i, temp_i)?;
    phys.add_hout_bnd(bnd_o)?;
    phys.add_temp_bnd(bnd_b, temp_bb)?;
    phys.add_temp_bnd(bnd_t, temp_tt)?;
    phys.add_hcnt_itf(itf_b)?;
    phys.add_hcnt_itf(itf_t)?;
    phys.add_flow_dom(dom_c, vel, pres, den, visc, fce)?;
    phys.add_vel_bnd(bnd_i, vel_i)?;
    phys.add_pres_bnd(bnd_o, pres_o)?;
    phys.add_vel_bnd(bnd_bc, vel_bc)?;
    phys.add_vel_bnd(bnd_tc, vel_tc)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 0.5)?;

    Ok(())
}
