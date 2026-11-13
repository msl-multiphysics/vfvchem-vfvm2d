use vfvchem_vfvm2d::*;

/// Steady-state heat and momentum transfer past a backward-facing step.
/// Run with: `cargo run --release --example heatflow_step`
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
/// Heat Transfer
/// -------------
/// 
/// Domain Properties:
/// dom_f - density (rho = 1.2 kg m-3)
/// dom_f - heat capacity (cp = 1005.0 J kg-1 K-1)
/// dom_f - thermal conductivity (k = 0.025 W m-1 K-1)
/// dom_f - heat source (Q = 0.0 W m-3)
///
/// Boundary Conditions:
/// bnd_i - temperature (T = 300.0 K)
/// bnd_o - outflow
/// bnd_t - temperature (T = 400.0 K)
/// bnd_b - no-flux (q = 0.0 W m-2)
/// 
/// Momentum Transfer
/// -----------------
/// 
/// Domain Properties:
/// dom - density (rho = 1.2 kg m-3)
/// dom - viscosity (mu = 1.8e-5 Pa s)
/// dom - body force (fce = <0.0, 0.0> N m-3)
///
/// Boundary Conditions:
/// bnd_i - velocity (parabolic; average velocity = 0.05 m s-1)
/// bnd_o - pressure (p = 0.0 Pa)
/// bnd_t - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_b - no slip (vel = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", 0.001, 0.002, "examples/output_heatflow_step")?;

    // geometry
    let dom_f = vars.add_dom("dom_f")?;  // flow channel
    let bnd_i = vars.add_bnd(dom_f, "bnd_i")?;  // inlet
    let bnd_o = vars.add_bnd(dom_f, "bnd_o")?;  // outlet
    let bnd_t = vars.add_bnd(dom_f, "bnd_t")?;  // top boundary
    let bnd_b = vars.add_bnd(dom_f, "bnd_b")?;  // bottom boundary

    // unknowns
    let temp_f = vars.add_scldom_unkn(dom_f, 0.0, "temp")?;       // temperature (T)
    let vel_f = vars.add_vecdom_unkn(dom_f, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres_f = vars.add_scldom_unkn(dom_f, 0.0, "pres")?;       // pressure (p)

    // domain values
    let den_f = vars.add_scldom_cnst(dom_f, 1.2, "den")?;         // density (rho)
    let hcap_f = vars.add_scldom_cnst(dom_f, 1005.0, "hcap")?;    // heat capacity (cp)
    let cond_f = vars.add_scldom_cnst(dom_f, 0.025, "cond")?;     // thermal conductivity (k)
    let hlin_f = vars.add_scldom_cnst(dom_f, 0.0, "hlin")?;       // linear source (L)
    let hsrc_f = vars.add_scldom_cnst(dom_f, 0.0, "hsrc")?;       // heat source (Q)
    let visc_f = vars.add_scldom_cnst(dom_f, 1.8e-5, "visc")?;    // viscosity (mu)
    let fce_f = vars.add_vecdom_cnst(dom_f, (0.0, 0.0), "fce")?;  // body force (fce)

    // non-constant boundary values
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {  // t, (x, y), [unkn] -> (u, v)
        let u = 0.075 * 4.0 * (y - 0.05) * (0.10 - y) / (0.05 * 0.05);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![])?;

    // constant boundary values
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0)?;       // temperature
    let temp_t = vars.add_sclbnd_cnst(bnd_t, 400.0)?;       // temperature
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;         // heat flux (q)
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0)?;         // pressure
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0))?;   // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;   // velocity

    // physics
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom_f, temp_f, den_f, hcap_f, vel_f,cond_f, hlin_f, hsrc_f)?;  // heat transfer domain: dom, T, rho, cp, k, L, Q
    phys.add_temp_bnd(bnd_i, temp_i)?;  // temperature boundary: bnd, T
    phys.add_temp_bnd(bnd_t, temp_t)?;  // temperature boundary: bnd, T
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // heat flux boundary: bnd, q
    phys.add_hout_bnd(bnd_o)?;                   // heat outflow boundary: bnd
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
