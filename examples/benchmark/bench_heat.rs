use vfvchem_vfvm2d::*;

/// Nonlinear heat conduction with k(T) = k0 * (1 + beta * T).
/// Run with: `cargo run --release --example bench_heat`
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
/// dom - main domain (0.1 m x 0.1 m)
///
/// Boundaries:
/// bnd_l - left boundary
/// bnd_r - right boundary
/// bnd_b - bottom boundary
/// bnd_t - top boundary
///
/// Heat Transfer
/// -------------
/// 
/// Domain Properties:
/// dom - thermal conductivity (k = 1.0 * (1 + 0.01 * T) W m-1 K-1)
///
/// Boundary Conditions:
/// bnd_l - temperature (T from Kirchhoff transform of 1D slab)
/// bnd_r - temperature (T = 0 K)
/// bnd_b - no flux (q = 0.0 W m-2)
/// bnd_t - no flux (q = 0.0 W m-2)
///
fn main() -> Result<(), VFVChemError> {
    const LEN: f64 = 0.1;
    const TEMP_HOT: f64 = 100.0;
    const COND_0: f64 = 1.0;
    const BETA: f64 = 0.01;
    let phi_hot = TEMP_HOT + 0.5 * BETA * TEMP_HOT * TEMP_HOT;

    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, LEN, LEN, 0.005, 0.005, "examples/output_bench_heat")?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let temp = vars.add_scldom_unkn(dom, 0.5 * TEMP_HOT, "temp")?;

    // domain values; conductivity depends on the temperature unknown
    let cond = vars.add_scldom_func(dom, Box::new(move |_t, _pos, u| COND_0 * (1.0 + BETA * u[0])), vec![temp], "cond")?;
    let hlin = vars.add_scldom_cnst(dom, 0.0, "hlin")?;
    let hsrc = vars.add_scldom_cnst(dom, 0.0, "hsrc")?;

    // boundary values; T = (sqrt(1 + 2 beta phi) - 1) / beta with phi = phi_hot * (1 - x / L)
    let tset_l = vars.add_sclbnd_func(bnd_l, Box::new(move |_t, (x, _y), _u| {
        let phi = phi_hot * (1.0 - x / LEN);
        ((1.0 + 2.0 * BETA * phi).sqrt() - 1.0) / BETA
    }), vec![])?;
    let tset_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;

    // physics
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom, temp, cond, hlin, hsrc)?;
    phys.add_temp_bnd(bnd_l, tset_l)?;
    phys.add_temp_bnd(bnd_r, tset_r)?;
    phys.add_hflx_bnd(bnd_b, hflx_b)?;
    phys.add_hflx_bnd(bnd_t, hflx_t)?;

    // solver
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 1.0)?;

    Ok(())
}
