use vfvchem_vfvm2d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

const SPACE: f64 = 0.02;

// cargo run --release --example heatmassflow_step
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatmassflow_step").unwrap();

    // create variables
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", SPACE, SPACE, SPACE, 100, 42)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_i = vars.add_bnd(dom, "bnd_i")?;
    let bnd_o = vars.add_bnd(dom, "bnd_o")?;
    let bnd_w = vars.add_bnd(dom, "bnd_w")?;

    // unknown domain values
    let temp = vars.add_scldom_unkn(dom, 300.0, "examples/output_heatmassflow_step/temp.vtu")?;
    let conc_0 = vars.add_scldom_unkn(dom, 0.0, "examples/output_heatmassflow_step/conc_0.vtu")?;
    let conc_1 = vars.add_scldom_unkn(dom, 0.0, "examples/output_heatmassflow_step/conc_1.vtu")?;
    let vel = vars.add_vecdom_unkn(dom, (1.0, 0.0), "examples/output_heatmassflow_step/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_heatmassflow_step/pres.vtu")?;

    // constant domain values
    let vlcp = vars.add_scldom_cnst(dom, 1000.0, "")?;
    let cond = vars.add_scldom_cnst(dom, 0.026, "")?;
    let diff_00 = vars.add_scldom_cnst(dom, 1.0e-4, "")?;
    let diff_11 = vars.add_scldom_cnst(dom, 2.0e-4, "")?;
    let den = vars.add_scldom_cnst(dom, 1.2, "")?;
    let visc = vars.add_scldom_cnst(dom, 1.8e-5, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;

    // species 0 converts to species 1, first order in the reactant
    let msrc_0_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| { -1.0 * unkn[0] });
    let msrc_1_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {  1.0 * unkn[0] });
    let hsrc_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| { 1.0e4 * unkn[0] });
    let msrc_0 = vars.add_scldom_func(dom, msrc_0_func, vec![conc_0], "")?;
    let msrc_1 = vars.add_scldom_func(dom, msrc_1_func, vec![conc_0], "")?;
    let hsrc = vars.add_scldom_func(dom, hsrc_func, vec![conc_0], "")?;

    // parabolic inlet on y in [0.2, 0.4], maximum 0.0225
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 0.0225 * (y - 0.2) * (0.4 - y) / (0.2 * 0.2);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;

    // constant boundary values
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0, "")?;
    let hflx_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;
    let conc_0_i = vars.add_sclbnd_cnst(bnd_i, 1.0, "")?;
    let conc_1_i = vars.add_sclbnd_cnst(bnd_i, 0.0, "")?;
    let mflx_0_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;
    let mflx_1_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0), "")?;

    // transient heat, mass, and momentum transfer solver
    let mut phys = TransientHeatMassFlow::new(2)?;
    phys.add_heat_dom(dom, temp, vlcp, cond, hsrc)?;
    phys.add_temp_bnd(bnd_i, temp_i)?;
    phys.add_hout_bnd(bnd_o)?;
    phys.add_hflx_bnd(bnd_w, hflx_w)?;
    phys.add_mass_dom(0, dom, conc_0, HashMap::from([(0, diff_00)]), msrc_0)?;
    phys.add_mass_dom(1, dom, conc_1, HashMap::from([(1, diff_11)]), msrc_1)?;
    phys.add_conc_bnd(0, bnd_i, conc_0_i)?;
    phys.add_conc_bnd(1, bnd_i, conc_1_i)?;
    phys.add_mout_bnd(0, bnd_o)?;
    phys.add_mout_bnd(1, bnd_o)?;
    phys.add_mflx_bnd(0, bnd_w, mflx_0_w)?;
    phys.add_mflx_bnd(1, bnd_w, mflx_1_w)?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
    phys.add_vel_bnd(bnd_i, vel_i)?;
    phys.add_pres_bnd(bnd_o, pres_o)?;
    phys.add_vel_bnd(bnd_w, vel_w)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 1.0, 100, 1, 30, 1e-3, 0.5)?;

    Ok(())
}
