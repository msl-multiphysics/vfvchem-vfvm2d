use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

const SPACE: f64 = 0.02;

// cargo run --release --example flow_lid
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_lid").unwrap();

    // create variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, SPACE, 50, 42)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknown domain values
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "examples/output_flow_lid/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_flow_lid/pres.vtu")?;

    // constant domain values
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;
    let visc = vars.add_scldom_cnst(dom, 0.001, "")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "")?;
    
    // constant boundary values
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0), "")?;
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0), "")?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0), "")?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (4e-4, 0.0), "")?;

    // steady-state momentum transfer solver
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
    phys.add_vel_bnd(bnd_l, vel_l)?;
    phys.add_vel_bnd(bnd_r, vel_r)?;
    phys.add_vel_bnd(bnd_b, vel_b)?;
    phys.add_vel_bnd(bnd_t, vel_t)?;
    phys.set_pres_ref(dom, 0, 0.0)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;

    Ok(())
}
