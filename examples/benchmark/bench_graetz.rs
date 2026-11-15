use vfvchem_vfvm2d::*;

/// Graetz problem: thermal entry length of plane Poiseuille flow.
/// Run with: `cargo run --release --example bench_graetz`
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
/// Heat Transfer
/// -------------
/// 
/// Domain Properties:
/// dom - density (rho = 1000.0 kg m-3)
/// dom - heat capacity (cp = 4000.0 J kg-1 K-1)
/// dom - velocity (vel = <1.5 * 0.005 * (1 - (1 - y / 0.005)^2), 0.0>)
/// dom - thermal conductivity (k = 0.4 W m-1 K-1)
///
/// Boundary Conditions:
/// bnd_l - temperature (T = 300.0 K)
/// bnd_b - temperature (T = 350.0 K)
/// bnd_t - heat flux (q = 0.0 W m-2)
/// bnd_r - outflow
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.3, 0.005, 0.0005, 0.0005, "examples/output_bench_graetz")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // half channel
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // inlet
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // outlet
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // wall
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // symmetry plane

    // unknowns
    let temp = vars.add_scldom_unkn(dom, 300.0, "temp")?;  // temperature (T)

    // non-constant domain values
    let vel_func: VecDomFunc = Box::new(|_t, (_x, y), _unkn| {
        let eta = 1.0 - y / 0.005;
        (1.5 * 0.005 * (1.0 - eta * eta), 0.0)
    });
    let vel = vars.add_vecdom_func(dom, vel_func, vec![], "vel")?;  // velocity (vel); prescribed

    // constant domain values
    let den = vars.add_scldom_cnst(dom, 1000.0, "den")?;  // density (rho)
    let hcap = vars.add_scldom_cnst(dom, 4000.0, "hcap")?;  // heat capacity (cp)
    let cond = vars.add_scldom_cnst(dom, 1000.0 * 4000.0 * 0.005 * 0.01 / 500.0, "cond")?;  // thermal conductivity (k)
    let zero = vars.add_scldom_cnst(dom, 0.0, "zero")?;  // no linear source or heat source

    // constant boundary values
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 300.0)?;  // temperature
    let temp_b = vars.add_sclbnd_cnst(bnd_b, 350.0)?;  // temperature
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // heat flux

    // physics
    // no flow domain is added; the velocity is prescribed
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom, temp, den, hcap, vel, cond, zero, zero)?;  // heat transfer domain: dom, T, rho, cp, vel, k, L, Q
    phys.add_temp_bnd(bnd_l, temp_l)?;  // temperature boundary: bnd, T
    phys.add_temp_bnd(bnd_b, temp_b)?;  // temperature boundary: bnd, T
    phys.add_hflx_bnd(bnd_t, hflx_t)?;  // heat flux boundary: bnd, q
    phys.add_hout_bnd(bnd_r)?;  // outflow boundary: bnd

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 100, 1e-8, 1.0)?;  // physics solver

    Ok(())
}
