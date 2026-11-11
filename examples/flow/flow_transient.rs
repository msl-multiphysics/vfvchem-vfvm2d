use vfvchem_vfvm2d::*;

// cargo run --release --example flow_transient
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", 0.01, 0.02, "examples/output_flow_transient")?;

    // geometry
    let dom = vars.add_dom("dom_f")?;
    let bnd_i = vars.add_bnd(dom, "bnd_i")?;
    let bnd_w = vars.add_bnd(dom, "bnd_w")?;
    let bnd_o = vars.add_bnd(dom, "bnd_o")?;
    
    // unknowns
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "vel")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "pres")?;

    // domain values
    let den = vars.add_scldom_cnst(dom, 1.2, "den")?;
    let visc = vars.add_scldom_cnst(dom, 1.8e-5, "visc")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "fce")?;

    // boundary values
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 0.0225 * (y - 0.2) * (0.4 - y) / (0.2 * 0.2);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![])?;
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0))?;
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0)?;

    // physics
    let mut phys = TransientFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_w, vel_w)?;  // arguments: boundary, v
    phys.add_pres_bnd(bnd_o, pres_o)?;  // arguments: boundary, p

    // solver
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8, 0.1, 100, 1)?;

    Ok(())
}
