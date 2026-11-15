use vfvchem_vfvm2d::*;

/// Planar Poiseuille flow with viscosity mu(y).
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
/// dom - channel (0.1 m x 0.05 m)
///
/// Boundaries:
/// bnd_l - inlet
/// bnd_r - outlet
/// bnd_b - bottom wall
/// bnd_t - top wall
///
/// Momentum Transfer
/// -----------------
/// 
/// Domain Properties:
/// dom - density (rho = 1000.0 kg m-3)
/// dom - viscosity (mu = 0.1 * (1 + 3 y / H) Pa s)
///
/// Boundary Conditions:
/// bnd_l - velocity (fully developed u(y), v = 0)
/// bnd_r - pressure (p = 0 Pa)
/// bnd_b - no slip (u = 0, v = 0)
/// bnd_t - no slip (u = 0, v = 0)
///
fn main() -> Result<(), VFVChemError> {
    const LEN: f64 = 0.1;
    const HGT: f64 = 0.05;
    const GRAD: f64 = 8.0;      // G = -dp/dx (Pa m-1)
    const VISC_0: f64 = 0.1;    // mu_0 (Pa s)
    const VISC_A: f64 = 3.0;    // a; mu(H) / mu(0) = 4
    const DEN_0: f64 = 1000.0;  // rho_0 (kg m-3)

    // fully developed inlet profile for mu = mu_0 * (1 + a y / H)
    let vel_exact = |y: f64| {
        let shear_0 = GRAD * HGT / (1.0 + VISC_A).ln() - GRAD * HGT / VISC_A;
        let xi = 1.0 + VISC_A * y / HGT;
        let gha = GRAD * HGT / VISC_A;
        HGT / (VISC_A * VISC_0) * ((shear_0 + gha) * xi.ln() - gha * (xi - 1.0))
    };

    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, LEN, HGT, 0.005, 0.005, "examples/output_bench_flow")?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "vel")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "pres")?;

    // domain values
    let den = vars.add_scldom_cnst(dom, DEN_0, "den")?;
    let visc = vars.add_scldom_func(dom, Box::new(|_t, (_x, y), _u| VISC_0 * (1.0 + VISC_A * y / HGT)), vec![], "visc")?;
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0), "fce")?;

    // boundary values
    let vel_l = vars.add_vecbnd_func(bnd_l, Box::new(move |_t, (_x, y), _u| (vel_exact(y), 0.0)), vec![])?;
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0))?;
    let pres_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;

    // physics
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;
    phys.add_vel_bnd(bnd_l, vel_l)?;
    phys.add_vel_bnd(bnd_b, vel_b)?;
    phys.add_vel_bnd(bnd_t, vel_t)?;
    phys.add_pres_bnd(bnd_r, pres_r)?;

    // solver
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 1.0)?;

    Ok(())
}
