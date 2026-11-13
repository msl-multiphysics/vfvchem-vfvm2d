use vfvchem_vfvm2d::*;
use std::collections::HashMap;

/// Transient buoyancy-driven flow in a square cavity.
/// Run with: `cargo run --release --example massflow_bouy`
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
/// Mass Transfer
/// -------------
/// 
/// Domain Properties:
/// dom - density (rho = 1.0 * (1.0 + 0.05 * c))
/// dom - diffusion coefficient (D = 1.0e-4 m2 s-1)
/// dom - molar source (R = 0.0 mol kg-1 s-1)
///
/// Boundary Conditions:
/// bnd_l - concentration (c = 1.0 mol kg-1)
/// bnd_r - concentration (c = 0.0 mol kg-1)
/// bnd_b - no flux (N = 0.0 mol m-2 s-1)
/// bnd_t - no flux (N = 0.0 mol m-2 s-1)
///
/// Initial Conditions:
/// dom - concentration (c = 0.5 mol kg-1)
///
/// Momentum Transfer
/// -----------------
/// 
/// Domain Properties:
/// dom - density (rho = 1.0 * (1.0 + 0.05 * c) kg m-3)
/// dom - viscosity (mu = 1.0e-3 Pa s)
/// dom - body force (fce = <0.0, -9.81 * (rho - 1.0)>)
///
/// Boundary Conditions:
/// bnd_l - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_r - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_b - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_t - no slip (vel = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.1, 0.1, 0.002, 0.004, "examples/output_massflow_bouy")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // domain
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top boundary

    // unknowns
    let conc = vars.add_scldom_unkn(dom, 0.5, "conc")?;       // concentration (c)
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres = vars.add_scldom_unkn(dom, 0.0, "pres")?;       // pressure (p)

    // non-constant domain values
    let den_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {  // t, (x, y), [unkn] -> rho
        1.0 * (1.0 + 0.05 * unkn[0])
    });
    let den = vars.add_scldom_func(dom, den_func, vec![conc], "den")?;  // density (rho)
    let fce_func: VecDomFunc = Box::new(|_t, (_x, _y), unkn| {
        let rho = 1.0 * (1.0 + 0.05 * unkn[0]);
        (0.0, -9.81 * (rho - 1.0))
    });
    let fce = vars.add_vecdom_func(dom, fce_func, vec![conc], "fce")?;  // body force (fce)

    // constant domain values
    let diff = vars.add_scldom_cnst(dom, 1.0e-4, "diff")?;  // diffusion coefficient (D)
    let nlin = vars.add_scldom_cnst(dom, 0.0, "nlin")?;     // linear source (L)
    let nsrc = vars.add_scldom_cnst(dom, 0.0, "nsrc")?;     // molar source (S)
    let visc = vars.add_scldom_cnst(dom, 1.0e-3, "visc")?;  // viscosity (mu)

    // boundary values
    let conc_l = vars.add_sclbnd_cnst(bnd_l, 1.0)?;  // concentration
    let conc_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;  // concentration
    let nflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // molar flux (N)
    let nflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // molar flux
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0))?;  // velocity
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0))?;  // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;  // velocity
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0))?;  // velocity

    // physics
    let mut phys = TransientMoleFlow::new(1)?;  // number of species
    phys.add_mole_dom(0, dom, conc, den, vel, HashMap::from([(0, diff)]), nlin, nsrc)?;  // mass transfer domain: comp_id, dom, c, rho, vel, D, L, R
    phys.add_conc_bnd(0, bnd_l, conc_l)?;  // concentration boundary: comp_id, bnd, c
    phys.add_conc_bnd(0, bnd_r, conc_r)?;  // concentration boundary: comp_id, bnd, c
    phys.add_nflx_bnd(0, bnd_b, nflx_b)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(0, bnd_t, nflx_t)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;  // momentum transfer domain: dom, vel, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_r, vel_r)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_b, vel_b)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_t, vel_t)?;  // velocity boundary: bnd, vel
    phys.set_pres_ref(dom, 0, 0.0)?;  // sets pressure at domain node 0 (arbitrary node) to 0.0

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 1.0, 0.1, 200, 1)?;  // physics solver

    Ok(())
}
