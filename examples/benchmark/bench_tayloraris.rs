use std::collections::HashMap;
use vfvchem_vfvm2d::*;

/// Taylor-Aris dispersion of a solute pulse in plane Poiseuille flow.
/// Run with: `cargo run --release --example bench_tayloraris`
///
/// Geometry
/// --------
/// 
///           bnd_t
///       +-----------+
/// bnd_l |    dom    | bnd_r
///       +-----------+
///           bnd_b
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
/// dom - density (rho = 1000.0 kg m-3)
/// dom - velocity (vel = <1.5 * 2e-4 * (1 - (1 - y / 0.005)^2) - 2e-4, 0.0>)
/// dom - diffusion coefficient (D = 1e-7 m2 s-1)
///
/// Boundary Conditions:
/// bnd_l - no flux (N = 0.0 mol m-2 s-1)
/// bnd_r - no flux (N = 0.0 mol m-2 s-1)
/// bnd_b - no flux (N = 0.0 mol m-2 s-1)
/// bnd_t - no flux (N = 0.0 mol m-2 s-1)
///
/// Initial Condition:
/// dom - concentration (c = exp(-x^2 / (2 * 0.01^2)))
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(-0.14, 0.0, 0.14, 0.005, 0.0005, 0.0005, "examples/output_bench_tayloraris")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // half channel
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left end
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right end
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // wall
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // symmetry plane

    // non-uniform unknowns
    let conc_func: SclPosFunc = Box::new(|x, _y| {
        (-x * x / (2.0 * 0.01 * 0.01)).exp()
    });
    let conc = vars.add_scldom_unkn_nonuniform(dom, conc_func, "conc")?;  // concentration (c)

    // non-constant domain values
    let vel_func: VecDomFunc = Box::new(|_t, (_x, y), _unkn| {
        let eta = 1.0 - y / 0.005;
        (1.5 * 2e-4 * (1.0 - eta * eta) - 2e-4, 0.0)
    });
    let vel = vars.add_vecdom_func(dom, vel_func, vec![], "vel")?;  // velocity (vel); prescribed

    // constant domain values
    let den = vars.add_scldom_cnst(dom, 1000.0, "den")?;  // density (rho)
    let diff = vars.add_scldom_cnst(dom, 1e-7, "diff")?;  // diffusion coefficient (D)
    let zero = vars.add_scldom_cnst(dom, 0.0, "zero")?;  // no linear source or source

    // constant boundary values
    let nflx_l = vars.add_sclbnd_cnst(bnd_l, 0.0)?;  // mass flux
    let nflx_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;  // mass flux
    let nflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // mass flux
    let nflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // mass flux

    // physics
    let mut phys = TransientMoleFlow::new(1)?;  // one species
    phys.add_mole_dom(0, dom, conc, den, vel, HashMap::from([(0, diff)]), zero, zero)?;  // mass transfer domain: species, dom, c, rho, vel, {j: D_ij}, L, S
    phys.add_nflx_bnd(0, bnd_l, nflx_l)?;  // mass flux boundary: species, bnd, N
    phys.add_nflx_bnd(0, bnd_r, nflx_r)?;  // mass flux boundary: species, bnd, N
    phys.add_nflx_bnd(0, bnd_b, nflx_b)?;  // mass flux boundary: species, bnd, N
    phys.add_nflx_bnd(0, bnd_t, nflx_t)?;  // mass flux boundary: species, bnd, N

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 50, 1e-8, 1.0, 5.0, 200, 20)?;  // 1000 s in steps of 5 s, output every 20 steps

    Ok(())
}
