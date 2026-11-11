use vfvchem_vfvm2d::*;

// cargo run --release --example flow_tworeg
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_obj("examples/obj/obj_tworeg.obj", 0.01, 0.02, "examples/output_flow_tworeg")?;

    // geometry
    let dom_lo = vars.add_dom("dom_lo")?;
    let dom_up = vars.add_dom("dom_up")?;
    let bnd_b = vars.add_bnd(dom_lo, "bnd_b")?;
    let bnd_r = vars.add_bnd(dom_lo, "bnd_r")?;
    let bnd_t = vars.add_bnd(dom_up, "bnd_t")?;
    let bnd_l = vars.add_bnd(dom_up, "bnd_l")?;
    let itf_d = vars.add_itf(dom_lo, dom_up, "itf_d")?;

    // unknowns
    let vel_lo = vars.add_vecdom_unkn(dom_lo, (0.0, 0.0), "vel")?;
    let vel_up = vars.add_vecdom_unkn(dom_up, (0.0, 0.0), "vel")?;
    let pres_lo = vars.add_scldom_unkn(dom_lo, 0.0, "pres")?;
    let pres_up = vars.add_scldom_unkn(dom_up, 0.0, "pres")?;

    // domain values
    let den_lo = vars.add_scldom_cnst(dom_lo, 1000.0, "den")?;
    let den_up = vars.add_scldom_cnst(dom_up, 1000.0, "den")?;
    let visc_lo = vars.add_scldom_cnst(dom_lo, 0.001, "visc")?;
    let visc_up = vars.add_scldom_cnst(dom_up, 0.001, "visc")?;
    let fce_lo = vars.add_vecdom_cnst(dom_lo, (0.0, 0.0), "fce")?;
    let fce_up = vars.add_vecdom_cnst(dom_up, (0.0, 0.0), "fce")?;

    // boundary values
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0))?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (4e-4, 0.0))?;
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0))?;

    // physics
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom_lo, vel_lo, pres_lo, den_lo, visc_lo, fce_lo)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_flow_dom(dom_up, vel_up, pres_up, den_up, visc_up, fce_up)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_b, vel_b)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_r, vel_r)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_t, vel_t)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_l, vel_l)?;  // arguments: boundary, v
    phys.add_fcnt_itf(itf_d)?;  // arguments: interface
    phys.set_pres_ref(dom_lo, 0, 0.0)?;  // arguments: domain, node_index, p

    // solver
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;

    Ok(())
}
