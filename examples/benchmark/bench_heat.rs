use vfvchem_vfvm2d::*;

/// Heat conduction with non-constant thermal conductivity.
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
/// dom - main domain
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
/// dom - thermal conductivity (k = 1.0 * (1 + 0.01 * T))
///
/// Boundary Conditions:
/// bnd_l - temperature (T from Kirchhoff transform)
/// bnd_r - temperature (T = 0 K)
/// bnd_b - no flux (q = 0.0 W m-2)
/// bnd_t - no flux (q = 0.0 W m-2)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.1, 0.1, 0.005, 0.005, "examples/output_bench_heat")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // domain
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top boundary

    // unknowns
    let temp = vars.add_scldom_unkn(dom, 50.0, "temp")?;  // temperature (T)

    // non-constant domain values
    let cond_func: SclDomFunc = Box::new(|_t, _pos, unkn| {
        1.0 * (1.0 + 0.01 * unkn[0])
    });
    let cond = vars.add_scldom_func(dom, cond_func, vec![temp], "cond")?;  // thermal conductivity (k)

    // constant domain values
    let hlin = vars.add_scldom_cnst(dom, 0.0, "hlin")?;  // linear source (L)
    let hsrc = vars.add_scldom_cnst(dom, 0.0, "hsrc")?;  // heat source (Q)

    // non-constant boundary values
    // T = (sqrt(1 + 2 beta phi) - 1) / beta with phi = phi_hot * (1 - x / L)
    // phi_hot = T_hot + 0.5 * beta * T_hot^2 = 100 + 0.5 * 0.01 * 100^2 = 150
    let temp_func: SclBndFunc = Box::new(|_t, (x, _y), _unkn| {
        let phi = 150.0 * (1.0 - x / 0.1);
        ((1.0 + 2.0 * 0.01 * phi).sqrt() - 1.0) / 0.01
    });
    let temp_l = vars.add_sclbnd_func(bnd_l, temp_func, vec![])?;  // temperature

    // constant boundary values
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;  // temperature
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // heat flux (q)
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // heat flux

    // physics
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom, temp, cond, hlin, hsrc)?;  // heat transfer domain: dom, T, k, L, Q
    phys.add_temp_bnd(bnd_l, temp_l)?;  // temperature boundary: bnd, T
    phys.add_temp_bnd(bnd_r, temp_r)?;  // temperature boundary: bnd, T
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_t, hflx_t)?;  // heat flux boundary: bnd, q

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 1.0)?;  // physics solver

    Ok(())
}
