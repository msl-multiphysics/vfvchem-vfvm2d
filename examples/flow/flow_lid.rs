use vfvchem_vfvm2d::*;

/// Steady-state lid-driven cavity flow in a square domain.
/// Run with: `cargo run --release --example flow_lid`
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
/// dom - viscosity (mu = 0.001 Pa s)
/// dom - body force (fce = <0.0, 0.0> N m-3)
///
/// Boundary Conditions:
/// bnd_l - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_r - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_b - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_t - velocity (vel = <0.01, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.1, 0.1, 0.001, 0.002, "examples/output_flow_lid")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // domain
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top boundary

    // unknowns
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres = vars.add_scldom_unkn(dom, 0.0, "pres")?;       // pressure (p)

    // domain values
    let den = vars.add_scldom_cnst(dom, 1000.0, "den")?;      // density (rho)
    let visc = vars.add_scldom_cnst(dom, 0.001, "visc")?;     // viscosity (mu)
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "fce")?;  // body force (fce)

    // boundary values
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0))?;   // velocity
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0))?;   // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;   // velocity
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.01, 0.0))?;  // velocity

    // physics
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;  // momentum transfer domain: dom, vel, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_r, vel_r)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_b, vel_b)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_t, vel_t)?;  // velocity boundary: bnd, vel
    phys.set_pres_ref(dom, 0, 0.0)?;  // sets pressure at domain node 0 (arbitrary node) to 0.0

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;  // physics solver

    Ok(())
}
