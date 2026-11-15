use vfvchem_vfvm2d::*;
use std::collections::HashMap;

/// Steady-state mass and momentum transfer past a backward-facing step.
/// Run with: `cargo run --release --example moleflow_step`
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
/// Mass Transfer
/// -------------
/// 
/// Domain Properties:
/// dom_f - density (rho = 1.2 kg m-3)
/// dom_f - diffusion coefficient (D00 = 1.0e-5 m2 s-1; D01 = 1.0e-5 m2 s-1; D11 = 2.0e-5 m2 s-1)
/// dom_f - molar source (R0 = -0.5 * c0 * c1 ; R1 = -0.2 * c0 * c1)
///
/// Boundary Conditions:
/// bnd_i - concentration (c0 = 2.0 mol kg-1; c1 = 1.0 mol kg-1)
/// bnd_o - outflow
/// bnd_t - no-flux (N = 0.0 mol m-2 s-1)
/// bnd_b - no-flux (N = 0.0 mol m-2 s-1)
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
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", 0.001, 0.002, "examples/output_moleflow_step")?;

    // geometry
    let dom_f = vars.add_dom("dom_f")?;  // flow channel
    let bnd_i = vars.add_bnd(dom_f, "bnd_i")?;  // inlet
    let bnd_o = vars.add_bnd(dom_f, "bnd_o")?;  // outlet
    let bnd_t = vars.add_bnd(dom_f, "bnd_t")?;  // top boundary
    let bnd_b = vars.add_bnd(dom_f, "bnd_b")?;  // bottom boundary

    // unknowns
    let conc0_f = vars.add_scldom_unkn(dom_f, 0.0, "conc0")?;  // concentration (c0)
    let conc1_f = vars.add_scldom_unkn(dom_f, 0.0, "conc1")?;  // concentration (c1)
    let vel_f = vars.add_vecdom_unkn(dom_f, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres_f = vars.add_scldom_unkn(dom_f, 0.0, "pres")?;       // pressure (p)

    // non-constant domain values
    // split reaction rates into R_i = L_i * c_i + S_i to improve stability
    // R0 = (-0.5 * c1) * c0
    // R1 = (-0.2 * c0) * c1
    let nlin0_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {
        -0.5 * unkn[0]
    });
    let nlin0_f = vars.add_scldom_func(dom_f, nlin0_func, vec![conc1_f], "nlin0")?;  // linear source (L0)
    let nlin1_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {
        -0.2 * unkn[0]
    });
    let nlin1_f = vars.add_scldom_func(dom_f, nlin1_func, vec![conc0_f], "nlin1")?;  // linear source (L1)

    // constant domain values
    let den_f = vars.add_scldom_cnst(dom_f, 1.2, "den")?;         // density (rho)
    let diff00_f = vars.add_scldom_cnst(dom_f, 1.0e-5, "d00")?;   // diffusion coefficient (D00)
    let diff01_f = vars.add_scldom_cnst(dom_f, 1.0e-5, "d01")?;   // diffusion coefficient (D01)
    let nsrc0_f = vars.add_scldom_cnst(dom_f, 0.0, "nsrc0")?;     // molar source (S0)
    let diff11_f = vars.add_scldom_cnst(dom_f, 2.0e-5, "d11")?;   // diffusion coefficient (D11)
    let nsrc1_f = vars.add_scldom_cnst(dom_f, 0.0, "nsrc1")?;     // molar source (S1)
    let visc_f = vars.add_scldom_cnst(dom_f, 1.8e-5, "visc")?;    // viscosity (mu)
    let fce_f = vars.add_vecdom_cnst(dom_f, (0.0, 0.0), "fce")?;  // body force (fce)

    // non-constant boundary values
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {  // t, (x, y), [unkn] -> (u, v)
        let u = 0.075 * 4.0 * (y - 0.05) * (0.10 - y) / (0.05 * 0.05);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![])?;

    // constant boundary values
    let conc0_i = vars.add_sclbnd_cnst(bnd_i, 2.0)?;  // concentration
    let nflx0_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // molar flux (N0)
    let nflx0_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // molar flux
    let conc1_i = vars.add_sclbnd_cnst(bnd_i, 1.0)?;  // concentration
    let nflx1_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // molar flux (N1)
    let nflx1_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // molar flux
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0)?;        // pressure
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0))?;  // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;  // velocity

    // physics
    let mut phys = TransientMoleFlow::new(2)?;  // number of species
    phys.add_mole_dom(0, dom_f, conc0_f, den_f, vel_f, HashMap::from([(0, diff00_f), (1, diff01_f)]), nlin0_f, nsrc0_f)?;  // mass transfer domain: comp_id, dom, c, rho, D, L, R
    phys.add_conc_bnd(0, bnd_i, conc0_i)?;  // concentration boundary: comp_id, bnd, c
    phys.add_nflx_bnd(0, bnd_t, nflx0_t)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(0, bnd_b, nflx0_b)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nout_bnd(0, bnd_o)?;                    // molar outflow boundary: bnd, N
    phys.add_mole_dom(1, dom_f, conc1_f, den_f, vel_f, HashMap::from([(1, diff11_f)]), nlin1_f, nsrc1_f)?;  // mass transfer domain: comp_id, dom, c, rho, D, L, R
    phys.add_conc_bnd(1, bnd_i, conc1_i)?;  // concentration boundary: comp_id, bnd, c
    phys.add_nflx_bnd(1, bnd_t, nflx1_t)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(1, bnd_b, nflx1_b)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nout_bnd(1, bnd_o)?;                    // molar outflow boundary: bnd, N
    phys.add_flow_dom(dom_f, vel_f, pres_f, den_f, visc_f, fce_f)?;  // momentum transfer domain: dom, vel, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;     // velocity boundary: bnd, vel
    phys.add_pres_bnd(bnd_o, pres_o)?;  // pressure boundary: bnd, p
    phys.add_vel_bnd(bnd_t, vel_t)?;     // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_b, vel_b)?;     // velocity boundary: bnd, vel

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 1.0, 0.1, 100, 1)?;  // physics solver

    Ok(())
}
