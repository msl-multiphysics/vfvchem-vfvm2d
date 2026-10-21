use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

const SPACE: f64 = 0.02;

// cargo run --release --example flow_step
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_step").unwrap();

    // create variables
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", SPACE, 100, 42)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_i = vars.add_bnd(dom, "bnd_i")?;
    let bnd_o = vars.add_bnd(dom, "bnd_o")?;
    let bnd_w = vars.add_bnd(dom, "bnd_w")?;

    // unknown domain values
    let vel = vars.add_vecdom_unkn(dom, (0.015, 0.0), "examples/output_flow_step/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_flow_step/pres.vtu")?;

    // constant domain values
    let den = vars.add_scldom_cnst(dom, 1.2, "")?;
    let visc = vars.add_scldom_cnst(dom, 1.8e-5, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;

    // parabolic inlet on y in [0.2, 0.4], maximum 0.0225
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 0.0225 * (y - 0.2) * (0.4 - y) / (0.2 * 0.2);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;

    // constant boundary values
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0), "")?;

    // transient momentum transfer solver
    let mut phys = TransientFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
    phys.add_vel_bnd(bnd_i, vel_i)?;
    phys.add_pres_bnd(bnd_o, pres_o)?;
    phys.add_vel_bnd(bnd_w, vel_w)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10.0, 100, 1, 30, 1e-3, 0.5)?;

    Ok(())
}
