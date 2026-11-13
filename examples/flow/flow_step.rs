use vfvchem_vfvm2d::*;

/// Steady-state flow past a backward-facing step.
/// Run with: `cargo run --release --example flow_step`
///
/// Geometry
/// --------
/// 
///            bnd_t
///       +-------------+
/// bnd_i |             |
///       +---+  dom_f  | bnd_o
///           |         |
///           +---------+
///              bnd_b
/// 
/// Domains:
/// dom_f - flow channel
///
/// Boundaries:
/// bnd_i - inlet
/// bnd_o - outlet
/// bnd_t - top boundary
/// bnd_b - bottom boundary
///
/// Momentum Transfer
/// -----------------
/// 
/// Domain Properties:
/// dom_f - density (rho = 1.2 kg m-3)
/// dom_f - viscosity (mu = 1.8e-5 Pa s)
/// dom_f - body force (fce = <0.0, 0.0> N m-3)
///
/// Boundary Conditions:
/// bnd_i - velocity (parabolic; average velocity = 0.05 m s-1)
/// bnd_o - pressure (p = 0.0 Pa)
/// bnd_t - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_b - no slip (vel = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", 0.002, 0.005, "examples/output_flow_step")?;

    // geometry
    let dom_f = vars.add_dom("dom_f")?;  // flow channel
    let bnd_i = vars.add_bnd(dom_f, "bnd_i")?;  // inlet
    let bnd_o = vars.add_bnd(dom_f, "bnd_o")?;  // outlet
    let bnd_t = vars.add_bnd(dom_f, "bnd_t")?;  // top boundary
    let bnd_b = vars.add_bnd(dom_f, "bnd_b")?;  // bottom boundary

    // unknowns
    let vel_f = vars.add_vecdom_unkn(dom_f, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres_f = vars.add_scldom_unkn(dom_f, 0.0, "pres")?;       // pressure (p)

    // domain values
    let den_f = vars.add_scldom_cnst(dom_f, 1.0, "den")?;         // density (rho)
    let visc_f = vars.add_scldom_cnst(dom_f, 1.0e-3, "visc")?;    // viscosity (mu)
    let fce_f = vars.add_vecdom_cnst(dom_f, (0.0, 0.0), "fce")?;  // body force (fce)

    // non-constant boundary values
    // define non-constant values with vars.add_vecbnd_func
    // function takes in time, position, and unknown scalars and returns the value at the boundary
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {  // t, (x, y), [unkn] -> (u, v)
        let u = 0.75 * 4.0 * (y - 0.05) * (0.10 - y) / (0.05 * 0.05);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![])?;  // unkn_id is empty since unknown scalars are not used

    // constant boundary values
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0)?;         // pressure
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0))?;   // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;   // velocity

    // physics
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom_f, vel_f, pres_f, den_f, visc_f, fce_f)?;  // momentum transfer domain: dom, vel, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;     // velocity boundary: bnd, vel
    phys.add_pres_bnd(bnd_o, pres_o)?;  // pressure boundary: bnd, p
    phys.add_vel_bnd(bnd_t, vel_t)?;     // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_b, vel_b)?;     // velocity boundary: bnd, vel

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 1.0)?;  // physics solver

    Ok(())
}
