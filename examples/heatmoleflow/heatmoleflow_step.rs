use vfvchem_vfvm2d::*;
use std::collections::HashMap;

/// Transient exothermic reacting flow past a backward-facing step.
/// Run with: `cargo run --release --example heatmoleflow_step`
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
/// dom_f - heat source (Q = 50000.0 * rho * k * c0)
///
/// Initial Conditions:
/// dom_f - temperature (T = 300.0 K)
///
/// Boundary Conditions:
/// bnd_i - temperature (T = 300.0 K)
/// bnd_o - outflow
/// bnd_t - no-flux (q = 0.0 W m-2)
/// bnd_b - no-flux (q = 0.0 W m-2)
///
/// Mole Transfer
/// -------------
///
/// Domain Properties:
/// dom_f - diffusion coefficient (D00 = D11 = 1.0e-4 m2 s-1)
/// dom_f - rate coefficient (k = 0.1 * exp[-10000.0 / 8.314 * (1 / T - 1 / 300.0)])
/// dom_f - molar sources (R0 = -rho * k * c0; R1 = rho * k * c0)
///
/// Initial Conditions:
/// dom_f - concentration (c0 = 0.5 mol kg-1; c1 = 0.0 mol kg-1)
///
/// Boundary Conditions:
/// bnd_i - concentration (c0 = 0.5 mol kg-1; c1 = 0.0 mol kg-1)
/// bnd_o - outflow
/// bnd_t - no-flux (N0 = 0.0 mol m-2 s-1; N1 = 0.0 mol m-2 s-1)
/// bnd_b - no-flux (N0 = 0.0 mol m-2 s-1; N1 = 0.0 mol m-2 s-1)
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
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", 0.001, 0.002, "examples/output_heatmoleflow_step")?;

    // geometry
    let dom_f = vars.add_dom("dom_f")?;  // flow channel
    let bnd_i = vars.add_bnd(dom_f, "bnd_i")?;  // inlet
    let bnd_o = vars.add_bnd(dom_f, "bnd_o")?;  // outlet
    let bnd_t = vars.add_bnd(dom_f, "bnd_t")?;  // top boundary
    let bnd_b = vars.add_bnd(dom_f, "bnd_b")?;  // bottom boundary

    // unknowns
    let temp_f = vars.add_scldom_unkn(dom_f, 300.0, "temp")?;     // temperature (T)
    let conc0_f = vars.add_scldom_unkn(dom_f, 0.5, "conc0")?;     // reactant concentration (c0)
    let conc1_f = vars.add_scldom_unkn(dom_f, 0.0, "conc1")?;     // product concentration (c1)
    let vel_f = vars.add_vecdom_unkn(dom_f, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres_f = vars.add_scldom_unkn(dom_f, 0.0, "pres")?;       // pressure (p)

    // non-constant domain values
    // split reaction rates into R_i = L_i * c_i + S_i to improve stability
    // R0 = (-rho * k(T)) * c0
    // R1 = rho * k(T) * c0
    let nlin0_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {
        let temp = unkn[0];
        let rate = 0.1 * (-10000.0 / 8.314 * (1.0 / temp - 1.0 / 300.0)).exp();
        -1.2 * rate
    });
    let nlin0_f = vars.add_scldom_func(dom_f, nlin0_func, vec![temp_f], "nlin0")?;  // linear molar source (L0)
    let nsrc1_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {
        let temp = unkn[0];
        let conc0 = unkn[1];
        let rate = 0.1 * (-10000.0 / 8.314 * (1.0 / temp - 1.0 / 300.0)).exp();
        1.2 * rate * conc0
    });
    let nsrc1_f = vars.add_scldom_func(dom_f, nsrc1_func, vec![temp_f, conc0_f], "nsrc1")?;  // molar source (S1)
    let hsrc_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {
        let temp = unkn[0];
        let conc0 = unkn[1];
        let rate = 0.1 * (-10000.0 / 8.314 * (1.0 / temp - 1.0 / 300.0)).exp();
        50000.0 * 1.2 * rate * conc0
    });
    let hsrc_f = vars.add_scldom_func(dom_f, hsrc_func, vec![temp_f, conc0_f], "hsrc")?;  // heat source (Q)

    // constant domain values
    let den_f = vars.add_scldom_cnst(dom_f, 1.2, "den")?;         // density (rho)
    let hcap_f = vars.add_scldom_cnst(dom_f, 1005.0, "hcap")?;    // heat capacity (cp)
    let cond_f = vars.add_scldom_cnst(dom_f, 0.025, "cond")?;     // thermal conductivity (k)
    let hlin_f = vars.add_scldom_cnst(dom_f, 0.0, "hlin")?;       // linear source (L)
    let diff00_f = vars.add_scldom_cnst(dom_f, 1.0e-4, "d00")?;   // diffusion coefficient (D00)
    let nsrc0_f = vars.add_scldom_cnst(dom_f, 0.0, "nsrc0")?;     // molar source (S0)
    let diff11_f = vars.add_scldom_cnst(dom_f, 1.0e-4, "d11")?;   // diffusion coefficient (D11)
    let nlin1_f = vars.add_scldom_cnst(dom_f, 0.0, "nlin1")?;     // linear molar source (L1)
    let visc_f = vars.add_scldom_cnst(dom_f, 1.8e-5, "visc")?;    // viscosity (mu)
    let fce_f = vars.add_vecdom_cnst(dom_f, (0.0, 0.0), "fce")?;  // body force (fce)

    // non-constant boundary values
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {  // t, (x, y), [unkn] -> (u, v)
        let u = 0.075 * 4.0 * (y - 0.05) * (0.10 - y) / (0.05 * 0.05);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![])?;

    // constant boundary values
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0)?;          // temperature
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;            // heat flux (q)
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;            // heat flux
    let conc0_i = vars.add_sclbnd_cnst(bnd_i, 0.5)?;           // reactant concentration
    let nflx0_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;           // reactant molar flux (N0)
    let nflx0_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;           // reactant molar flux
    let conc1_i = vars.add_sclbnd_cnst(bnd_i, 0.0)?;           // product concentration
    let nflx1_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;           // product molar flux (N1)
    let nflx1_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;           // product molar flux
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0)?;            // pressure
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0))?;      // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;      // velocity

    // physics
    let mut phys = TransientHeatMoleFlow::new(2)?;  // number of species
    phys.add_heat_dom(dom_f, temp_f, den_f, hcap_f, vel_f, cond_f, hlin_f, hsrc_f)?;  // heat transfer domain: dom, T, rho, cp, vel, k, L, Q
    phys.add_temp_bnd(bnd_i, temp_i)?;  // temperature boundary: bnd, T
    phys.add_hflx_bnd(bnd_t, hflx_t)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // heat flux boundary: bnd, q
    phys.add_hout_bnd(bnd_o)?;                  // heat outflow boundary: bnd
    phys.add_mole_dom(0, dom_f, conc0_f, den_f, vel_f, HashMap::from([(0, diff00_f)]), nlin0_f, nsrc0_f)?;  // mole transfer domain: comp_id, dom, c, rho, vel, D, L, S
    phys.add_conc_bnd(0, bnd_i, conc0_i)?;  // concentration boundary: comp_id, bnd, c
    phys.add_nflx_bnd(0, bnd_t, nflx0_t)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(0, bnd_b, nflx0_b)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nout_bnd(0, bnd_o)?;                    // molar outflow boundary: comp_id, bnd
    phys.add_mole_dom(1, dom_f, conc1_f, den_f, vel_f, HashMap::from([(1, diff11_f)]), nlin1_f, nsrc1_f)?;  // mole transfer domain: comp_id, dom, c, rho, vel, D, L, S
    phys.add_conc_bnd(1, bnd_i, conc1_i)?;  // concentration boundary: comp_id, bnd, c
    phys.add_nflx_bnd(1, bnd_t, nflx1_t)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(1, bnd_b, nflx1_b)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nout_bnd(1, bnd_o)?;                   // molar outflow boundary: comp_id, bnd
    phys.add_flow_dom(dom_f, vel_f, pres_f, den_f, visc_f, fce_f)?;  // momentum transfer domain: dom, vel, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;     // velocity boundary: bnd, vel
    phys.add_pres_bnd(bnd_o, pres_o)?;  // pressure boundary: bnd, p
    phys.add_vel_bnd(bnd_t, vel_t)?;     // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_b, vel_b)?;     // velocity boundary: bnd, vel

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 1.0, 0.1, 200, 2)?;  // physics solver

    Ok(())
}
