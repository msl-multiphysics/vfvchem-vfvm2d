use vfvchem_vfvm2d::*;

// cargo run --release --example flow_lid
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, 0.01, 0.02, "examples/output_flow_lid")?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "vel")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "pres")?;

    // domain values
    let den = vars.add_scldom_cnst(dom, 1000.0, "den")?;
    let visc = vars.add_scldom_cnst(dom, 0.001, "visc")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "fce")?;

    // boundary values
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0))?;
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0))?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (4e-4, 0.0))?;

    // physics
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_r, vel_r)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_b, vel_b)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_t, vel_t)?;  // arguments: boundary, v
    phys.set_pres_ref(dom, 0, 0.0)?;  // arguments: domain, node_index, p

    // solver
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;

    Ok(())
}
