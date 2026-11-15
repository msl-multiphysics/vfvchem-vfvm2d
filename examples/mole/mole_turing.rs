use vfvchem_vfvm2d::*;
use std::collections::HashMap;

/// Turing pattern formation (Brusselator) in a square domain.
/// Run with: `cargo run --release --example mole_turing`
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
/// dom - density (rho = 1.0 kg m-3)
/// dom - diffusion coefficient (D00 = 1.0e-6 m2 s-1; D11 = 1.0e-5 m2 s-1)
/// dom - molar source (R0 = 1.0 - 2.8 * c0 + c0^2 * c1 ; R1 = 1.8 * c0 - c0^2 * c1)
///
/// Boundary Conditions:
/// bnd_l - no flux (N = 0.0 mol m-2 s-1)
/// bnd_r - no flux (N = 0.0 mol m-2 s-1)
/// bnd_b - no flux (N = 0.0 mol m-2 s-1)
/// bnd_t - no flux (N = 0.0 mol m-2 s-1)
///
/// Initial Conditions:
/// dom - concentration (c0 = 1.0 mol kg-1 with 1% spatial perturbation; c1 = 1.8 mol kg-1)
/// 
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.1, 0.1, 0.001, 0.001, "examples/output_mole_turing")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // domain
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top boundary

    // unknowns
    // perturb conc0 to trigger the instability; uniform IC is an exact steady state
    let conc0 = vars.add_scldom_unkn_nonuniform(dom, Box::new(|x, y| {
        let noise = ((x * 12.9898 + y * 78.233).sin() * 43758.5453).fract();  // pseudo-random in [0, 1)
        1.0 + 0.01 * noise
    }), "conc0")?;  // concentration (c0)
    let conc1 = vars.add_scldom_unkn(dom, 1.8, "conc1")?;  // concentration (c1)

    // non-constant domain values
    // split reaction rates into R_i = L_i * c_i + S_i to improve stability
    // R0 = (-2.8 + c0 * c1) * c0 + 1.0
    // R1 =          (-c0^2) * c1 + 1.8 * c0
    let nlin0_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {
        -2.8 + unkn[0] * unkn[1]
    });
    let nlin0 = vars.add_scldom_func(dom, nlin0_func, vec![conc0, conc1], "nlin0")?;  // linear source (L0)
    let nlin1_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {
        - unkn[0] * unkn[0]
    });
    let nlin1 = vars.add_scldom_func(dom, nlin1_func, vec![conc0, conc1], "nlin1")?;  // linear source (L1)
    let nsrc1_func: SclDomFunc = Box::new(|_t, (_x, _y), unkn| {
        1.8 * unkn[0]
    });
    let nsrc1 = vars.add_scldom_func(dom, nsrc1_func, vec![conc0, conc1], "nsrc1")?;  // molar source (S1)

    // constant domain values
    let den = vars.add_scldom_cnst(dom, 1.0, "den")?;        // density (rho)
    let diff00 = vars.add_scldom_cnst(dom, 1.0e-6, "d00")?;  // diffusion coefficient (D00)
    let nsrc0 = vars.add_scldom_cnst(dom, 1.0, "nsrc0")?;    // molar source (S0)
    let diff11 = vars.add_scldom_cnst(dom, 1.0e-5, "d11")?;  // diffusion coefficient (D11)
    
    // boundary values
    // define boundary properties with vars.add_sclbnd_cnst
    let nflx0_l = vars.add_sclbnd_cnst(bnd_l, 0.0)?;  // molar flux (N0)
    let nflx0_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;  // molar flux
    let nflx0_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // molar flux
    let nflx0_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // molar flux
    let nflx1_l = vars.add_sclbnd_cnst(bnd_l, 0.0)?;  // molar flux (N1)
    let nflx1_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;  // molar flux
    let nflx1_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // molar flux
    let nflx1_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // molar flux

    // physics
    let mut phys = TransientMole::new(2)?;  // number of species
    phys.add_mole_dom(0, dom, conc0, den, HashMap::from([(0, diff00)]), nlin0, nsrc0)?;  // mass transfer domain: comp_id, dom, c, rho, D, L, R
    phys.add_nflx_bnd(0, bnd_l, nflx0_l)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(0, bnd_r, nflx0_r)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(0, bnd_b, nflx0_b)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(0, bnd_t, nflx0_t)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_mole_dom(1, dom, conc1, den, HashMap::from([(1, diff11)]), nlin1, nsrc1)?;  // mass transfer domain: comp_id, dom, c, rho, D, L, R
    phys.add_nflx_bnd(1, bnd_l, nflx1_l)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(1, bnd_r, nflx1_r)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(1, bnd_b, nflx1_b)?;  // molar flux boundary: comp_id, bnd, N
    phys.add_nflx_bnd(1, bnd_t, nflx1_t)?;  // molar flux boundary: comp_id, bnd, N

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 1.0, 0.5, 200, 2)?;  // physics solver

    Ok(())
}
