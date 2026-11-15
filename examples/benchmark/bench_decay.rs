use vfvchem_vfvm2d::*;

/// Adaptive time stepping with nonlinear decay.
/// Run with: `cargo run --release --example bench_decay`
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
/// dom - density (rho = 1.0 + 0.5 * x + 0.25 * y)
/// dom - heat capacity (cp = 2.0 J kg-1 K-1)
/// dom - thermal conductivity (k = 0.0 W m-1 K-1)
/// dom - heat source (Q = -rho * cp * T^2 / (1.0 + 0.5 * x))
///
/// Initial Conditions:
/// dom - temperature (T = 1.0 + 0.5 * x)
/// 
/// Boundary Conditions:
/// bnd_l - no flux (q = 0.0 W m-2)
/// bnd_r - no flux (q = 0.0 W m-2)
/// bnd_b - no flux (q = 0.0 W m-2)
/// bnd_t - no flux (q = 0.0 W m-2)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, 0.1, 0.1, "examples/output_bench_decay")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // domain
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top boundary

    // unknowns
    let temp = vars.add_scldom_unkn_nonuniform(dom, Box::new(|x, _y| 1.0 + 0.5 * x), "temp")?;  // temperature (T)

    // non-constant domain values
    let den_func: SclDomFunc = Box::new(|_t, (x, y), _unkn| {
        1.0 + 0.5 * x + 0.25 * y
    });
    let den = vars.add_scldom_func(dom, den_func, vec![], "den")?;  // density (rho)
    let hlin_func: SclDomFunc = Box::new(|_t, (x, y), unkn| {
        let temp = unkn[0];
        -(1.0 + 0.5 * x + 0.25 * y) * 2.0 * temp / (1.0 + 0.5 * x)
    });
    let hlin = vars.add_scldom_func(dom, hlin_func, vec![temp], "hlin")?;  // linear source (L)

    // constant domain values
    let hcap = vars.add_scldom_cnst(dom, 2.0, "hcap")?;  // heat capacity (cp)
    let cond = vars.add_scldom_cnst(dom, 0.0, "cond")?;  // thermal conductivity (k)
    let hsrc = vars.add_scldom_cnst(dom, 0.0, "hsrc")?;  // heat source (Q)

    // boundary values
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, 0.0)?;  // heat flux (q); positive if outward, negative if inward
    let hflx_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;  // heat flux
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // heat flux
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // heat flux

    // physics
    let mut phys = TransientHeat::new()?;
    phys.add_heat_dom(dom, temp, den, hcap, cond, hlin, hsrc)?;  // heat transfer domain: dom, T, rho, cp, k, L, Q
    phys.add_hflx_bnd(bnd_l, hflx_l)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_r, hflx_r)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_t, hflx_t)?;  // heat flux boundary: bnd, q

    // solver
    // max_iter, tol, damp, dt_init, dt_min, dt_max, t_end, num_step
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve_adaptive(&mut vars, &mut lin, 50, 1e-6, 1.0, 0.5, 1e-8, 2.0, 10.0, 10)?;  // physics solver

    Ok(())
}
