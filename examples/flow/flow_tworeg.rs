use vfvchem_vfvm2d::*;

/// Steady-state lid-driven cavity flow in a divided domain.
/// Run with: `cargo run --release --example flow_tworeg`
///
///          bnd_t
///       +---------+
///       |dom_tl / |
/// bnd_l |    /    | bnd_r
///       | / dom_br|
///       +---------+
///          bnd_b
///
/// Domains:
/// dom_br - bottom-right triangle
/// dom_tl - top-left triangle
///
/// Boundaries:
/// bnd_l - left boundary
/// bnd_r - right boundary
/// bnd_b - bottom boundary
/// bnd_t - top boundary
///
/// Interfaces:
/// itf_d - diagonal interface
///
/// Domain Properties:
/// dom_br - density (rho = 1000.0 kg m-3)
/// dom_br - viscosity (mu = 0.001 Pa s)
/// dom_br - body force (fce = <0.0, 0.0> N m-3)
/// dom_tl - density (rho = 1000.0 kg m-3)
/// dom_tl - viscosity (mu = 0.001 Pa s)
/// dom_tl - body force (fce = <0.0, 0.0> N m-3)
///
/// Boundary Conditions:
/// bnd_l - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_r - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_b - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_t - velocity (vel = <0.01, 0.0> m s-1)
///
/// Interface Conditions:
/// itf_d - velocity and traction continuity
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_obj("examples/obj/obj_tworeg.obj", 0.001, 0.002, "examples/output_flow_tworeg")?;

    // geometry
    let dom_br = vars.add_dom("dom_br")?;  // bottom-right triangle
    let dom_tl = vars.add_dom("dom_tl")?;  // top-left triangle
    let bnd_l = vars.add_bnd(dom_tl, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom_br, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom_br, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom_tl, "bnd_t")?;  // top boundary
    let itf_d = vars.add_itf(dom_br, dom_tl, "itf_d")?;  // diagonal interface

    // unknowns
    let vel_br = vars.add_vecdom_unkn(dom_br, (0.0, 0.0), "vel")?;  // velocity (vel)
    let vel_tl = vars.add_vecdom_unkn(dom_tl, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres_br = vars.add_scldom_unkn(dom_br, 0.0, "pres")?;  // pressure (p)
    let pres_tl = vars.add_scldom_unkn(dom_tl, 0.0, "pres")?;  // pressure (p)

    // domain values
    let den_br = vars.add_scldom_cnst(dom_br, 1000.0, "den")?;      // density (rho)
    let visc_br = vars.add_scldom_cnst(dom_br, 0.001, "visc")?;     // viscosity (mu)
    let fce_br = vars.add_vecdom_cnst(dom_br, (0.0, 0.0), "fce")?;  // body force (fce)
    let den_tl = vars.add_scldom_cnst(dom_tl, 1000.0, "den")?;      // density (rho)
    let visc_tl = vars.add_scldom_cnst(dom_tl, 0.001, "visc")?;     // viscosity (mu)
    let fce_tl = vars.add_vecdom_cnst(dom_tl, (0.0, 0.0), "fce")?;  // body force (fce)
    
    // boundary values
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0))?;   // velocity
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0))?;   // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;   // velocity
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.01, 0.0))?;  // velocity

    // physics
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom_br, vel_br, pres_br, den_br, visc_br, fce_br)?;  // momentum transfer domain: dom, vel, p, rho, mu, f
    phys.add_flow_dom(dom_tl, vel_tl, pres_tl, den_tl, visc_tl, fce_tl)?;  // momentum transfer domain: dom, vel, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_r, vel_r)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_b, vel_b)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_t, vel_t)?;  // velocity boundary: bnd, vel
    phys.add_fcnt_itf(itf_d)?;  // continuity interface: itf
    phys.set_pres_ref(dom_br, 0, 0.0)?;  // sets pressure at domain node 0 (arbitrary node) to 0.0

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;  // physics solver

    Ok(())
}
