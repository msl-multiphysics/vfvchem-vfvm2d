use vfvchem_vfvm2d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

const SPACE: f64 = 0.05;

// cargo run --release --example mass_square
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_mass_square").unwrap();

    // create variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, SPACE, 100, 42)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknown domain values
    let conc_0 = vars.add_scldom_unkn(dom, 0.0, "examples/output_mass_square/conc_0.vtu")?;
    let conc_1 = vars.add_scldom_unkn(dom, 0.0, "examples/output_mass_square/conc_1.vtu")?;

    // constant domain values
    // component 0 diffuses down its own gradient and the gradient of component 1
    let diff_00 = vars.add_scldom_cnst(dom, 1.0e-4, "")?;
    let diff_01 = vars.add_scldom_cnst(dom, 0.5e-4, "")?;
    let diff_11 = vars.add_scldom_cnst(dom, 1.0e-4, "")?;
    let msrc_0 = vars.add_scldom_cnst(dom, -0.2, "")?;
    let msrc_1 = vars.add_scldom_cnst(dom,  0.1, "")?;

    // constant boundary values
    let mflx_0_l = vars.add_sclbnd_cnst(bnd_l, 0.5, "")?;
    let conc_0_r = vars.add_sclbnd_cnst(bnd_r, 3e-4, "")?;
    let mflx_0_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let conc_0_t = vars.add_sclbnd_cnst(bnd_t, 4e-4, "")?;
    let conc_1_l = vars.add_sclbnd_cnst(bnd_l, 1e-4, "")?;
    let conc_1_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let mflx_1_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let mflx_1_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;

    // steady-state mass transfer solver
    let mut phys = SteadyMass::new(2)?;
    phys.add_mass_dom(0, dom, conc_0, HashMap::from([(0, diff_00), (1, diff_01)]), msrc_0)?;
    phys.add_mass_dom(1, dom, conc_1, HashMap::from([(1, diff_11)]), msrc_1)?;
    phys.add_mflx_bnd(0, bnd_l, mflx_0_l)?;
    phys.add_conc_bnd(0, bnd_r, conc_0_r)?;
    phys.add_mflx_bnd(0, bnd_b, mflx_0_b)?;
    phys.add_conc_bnd(0, bnd_t, conc_0_t)?;
    phys.add_conc_bnd(1, bnd_l, conc_1_l)?;
    phys.add_conc_bnd(1, bnd_r, conc_1_r)?;
    phys.add_mflx_bnd(1, bnd_b, mflx_1_b)?;
    phys.add_mflx_bnd(1, bnd_t, mflx_1_t)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;

    Ok(())
}
