use vfvchem_vfvm2d::*;
use std::fs::create_dir_all;

const SPACE: f64 = 0.001;

// cargo run --release --example heatflow_cyl
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatflow_cyl").unwrap();

    // create variables
    let mut vars = Variables::new_obj("examples/obj/obj_cyl.obj", SPACE, 100, 42)?;

    // geometry
    // each heated cylinder is its own domain with its own interface and no-slip boundary
    let dom_c = vars.add_dom("dom_c")?;
    let bnd_i = vars.add_bnd(dom_c, "bnd_i")?;
    let bnd_o = vars.add_bnd(dom_c, "bnd_o")?;
    let bnd_w = vars.add_bnd(dom_c, "bnd_w")?;
    let mut dom_h = Vec::new();
    let mut bnd_h = Vec::new();
    let mut itf_h = Vec::new();
    for k in 0..4 {
        let dom = vars.add_dom(&format!("dom_h{k}"))?;
        bnd_h.push(vars.add_bnd(dom_c, &format!("bnd_h{k}"))?);
        itf_h.push(vars.add_itf(dom_c, dom, &format!("itf_h{k}"))?);
        dom_h.push(dom);
    }

    // unknown domain values
    let temp_c = vars.add_scldom_unkn(dom_c, 300.0, "examples/output_heatflow_cyl/temp_c.vtu")?;
    let vel = vars.add_vecdom_unkn(dom_c, (1.0, 0.0), "examples/output_heatflow_cyl/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom_c, 0.0, "examples/output_heatflow_cyl/pres.vtu")?;
    let mut temp_h = Vec::new();
    for (k, &dom) in dom_h.iter().enumerate() {
        temp_h.push(vars.add_scldom_unkn(dom, 300.0, &format!("examples/output_heatflow_cyl/temp_h{k}.vtu"))?);
    }

    // constant domain values
    let vlcp_c = vars.add_scldom_cnst(dom_c, 1.0, "")?;
    let cond_c = vars.add_scldom_cnst(dom_c, 0.01, "")?;
    let hsrc_c = vars.add_scldom_cnst(dom_c, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom_c, 1.0, "")?;
    let visc = vars.add_scldom_cnst(dom_c, 0.01, "")?;
    let fce = vars.add_vecdom_cnst(dom_c, (0.0, 0.0), "")?;

    // parabolic inlet on y in [0.0, 0.05], maximum 1.0
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let u = 4.0 * 1.0 * y * (0.05 - y) / (0.05 * 0.05);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;

    // constant boundary values
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0, "")?;
    let hflx_w = vars.add_sclbnd_cnst(bnd_w, 0.0, "")?;
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0), "")?;

    // steady-state heat and momentum transfer solver
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom_c, temp_c, vlcp_c, cond_c, hsrc_c)?;
    phys.add_temp_bnd(bnd_i, temp_i)?;
    phys.add_hout_bnd(bnd_o)?;
    phys.add_hflx_bnd(bnd_w, hflx_w)?;
    phys.add_flow_dom(dom_c, vel, pres, den, visc, fce)?;
    phys.add_vel_bnd(bnd_i, vel_i)?;
    phys.add_pres_bnd(bnd_o, pres_o)?;
    phys.add_vel_bnd(bnd_w, vel_w)?;

    // heated cylinders
    for k in 0..4 {
        let vlcp = vars.add_scldom_cnst(dom_h[k], 1.0, "")?;
        let cond = vars.add_scldom_cnst(dom_h[k], 1.0, "")?;
        let hsrc = vars.add_scldom_cnst(dom_h[k], 1e4, "")?;
        let vel_h = vars.add_vecbnd_cnst(bnd_h[k], (0.0, 0.0), "")?;
        phys.add_heat_dom(dom_h[k], temp_h[k], vlcp, cond, hsrc)?;
        phys.add_hcnt_itf(itf_h[k])?;
        phys.add_vel_bnd(bnd_h[k], vel_h)?;
    }

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 0.5)?;

    Ok(())
}
