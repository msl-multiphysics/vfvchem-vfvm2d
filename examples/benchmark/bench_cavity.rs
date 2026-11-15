use vfvchem_vfvm2d::*;

/// Natural convection in a differentially heated square cavity (de Vahl Davis 1983).
/// Run with: `cargo run --release --example bench_cavity`
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
/// dom - density (rho = 1.0 kg m-3)
/// dom - heat capacity (cp = 1.0 J kg-1 K-1)
/// dom - thermal conductivity (k = 1.0 W m-1 K-1)
///
/// Boundary Conditions:
/// bnd_l - temperature (T = 1.0 K)
/// bnd_r - temperature (T = 0.0 K)
/// bnd_b - heat flux (q = 0.0 W m-2)
/// bnd_t - heat flux (q = 0.0 W m-2)
///
/// Momentum Transfer
/// -----------------
/// 
/// Domain Properties:
/// dom - density (rho = 1.0 kg m-3)
/// dom - viscosity (mu = 0.71 Pa s)
/// dom - body force (fce = <0.0, Ra Pr (T - 0.5)>; Ra = 1e4)
///
/// Boundary Conditions:
/// bnd_l - velocity (vel = <0.0, 0.0> m s-1)
/// bnd_r - velocity (vel = <0.0, 0.0> m s-1)
/// bnd_b - velocity (vel = <0.0, 0.0> m s-1)
/// bnd_t - velocity (vel = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, 0.0125, 0.0125, "examples/output_bench_cavity")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // square cavity
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // hot wall
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // cold wall
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom wall
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top wall

    // unknowns
    let temp = vars.add_scldom_unkn(dom, 0.5, "temp")?;  // temperature (T)
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres = vars.add_scldom_unkn(dom, 0.0, "pres")?;  // pressure (p)

    // non-constant domain values
    let fce_func: VecDomFunc = Box::new(|_t, (_x, _y), unkn| {
        (0.0, 1e4 * 0.71 * (unkn[0] - 0.5))
    });
    let fce = vars.add_vecdom_func(dom, fce_func, vec![temp], "fce")?;  // body force (fce)

    // constant domain values
    let den = vars.add_scldom_cnst(dom, 1.0, "den")?;  // density (rho)
    let hcap = vars.add_scldom_cnst(dom, 1.0, "hcap")?;  // heat capacity (cp)
    let cond = vars.add_scldom_cnst(dom, 1.0, "cond")?;  // thermal conductivity (k)
    let zero = vars.add_scldom_cnst(dom, 0.0, "zero")?;  // no linear source or heat source
    let visc = vars.add_scldom_cnst(dom, 0.71, "visc")?;  // viscosity (mu)

    // constant boundary values
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 1.0)?;  // temperature
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;  // temperature
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;  // heat flux
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // heat flux
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0))?;  // velocity
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0))?;  // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0))?;  // velocity
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (0.0, 0.0))?;  // velocity

    // physics
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom, temp, den, hcap, vel, cond, zero, zero)?;  // heat transfer domain: dom, T, rho, cp, vel, k, L, Q
    phys.add_temp_bnd(bnd_l, temp_l)?;  // temperature boundary: bnd, T
    phys.add_temp_bnd(bnd_r, temp_r)?;  // temperature boundary: bnd, T
    phys.add_hflx_bnd(bnd_b, hflx_b)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_t, hflx_t)?;  // heat flux boundary: bnd, q
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;  // momentum transfer domain: dom, vel, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_r, vel_r)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_b, vel_b)?;  // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_t, vel_t)?;  // velocity boundary: bnd, vel
    phys.set_pres_ref(dom, 0, 0.0)?;  // pressure reference: dom, vertex, p

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 500, 1e-8, 0.3)?;  // physics solver

    Ok(())
}
