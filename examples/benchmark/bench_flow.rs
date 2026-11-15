use vfvchem_vfvm2d::*;

/// Poiseuille flow with non-constant viscosity.
/// Run with: `cargo run --release --example bench_flow`
///
/// Geometry
/// --------
/// 
///          bnd_t
///       +---------+
///       |         |
/// bnd_l |   dom   | bnd_r
///       |         |
///       +---------+
///          bnd_b
///
/// Domains:
/// dom - main domain
///
/// Boundaries:
/// bnd_l - left boundary
/// bnd_r - right boundary
/// bnd_b - bottom boundary
/// bnd_t - top boundary
///
/// Momentum Transfer
/// -----------------
/// 
/// Domain Properties:
/// dom - density (rho = 1000.0 kg m-3)
/// dom - viscosity (mu = 0.1 * (1 + 3 y / 0.05))
///
/// Boundary Conditions:
/// bnd_l - velocity (fully developed profile)
/// bnd_r - pressure (p = 0 Pa)
/// bnd_b - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_t - no slip (vel = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.1, 0.05, 0.005, 0.005, "examples/output_bench_flow")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // domain
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top boundary

    // unknowns
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres = vars.add_scldom_unkn(dom, 0.0, "pres")?;  // pressure (p)

    // non-constant domain values
    let visc_func: SclDomFunc = Box::new(|_t, (_x, y), _unkn| {
        0.1 * (1.0 + 3.0 * y / 0.05)
    });
    let visc = vars.add_scldom_func(dom, visc_func, vec![], "visc")?;  // viscosity (mu)

    // constant domain values
    let den = vars.add_scldom_cnst(dom, 1000.0, "den")?;  // density (rho)
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "fce")?;  // body force (fce)

    // non-constant boundary values
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {
        let shear_0 = 8.0 * 0.05 / 4.0_f64.ln() - 8.0 * 0.05 / 3.0;
        let xi = 1.0 + 3.0 * y / 0.05;
        let gha = 8.0 * 0.05 / 3.0;
        let u = 0.05 / (3.0 * 0.1) * ((shear_0 + gha) * xi.ln() - gha * (xi - 1.0));
        (u, 0.0)
    });
    let vel_l = vars.add_vecbnd_func(bnd_l, vel_func, vec![])?;  // velocity

    // constant boundary values
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;  // velocity
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0))?;  // velocity
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;  // pressure

    // physics
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;  // flow domain: dom, vel, p, rho, mu, fce
    phys.add_vel_bnd(bnd_l, vel_l)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_b, vel_b)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_t, vel_t)?;  // velocity boundary: bnd, vel
    phys.add_pres_bnd(bnd_r, pres_r)?;  // pressure boundary: bnd, p

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 1.0)?;  // physics solver

    Ok(())
}
