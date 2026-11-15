use std::f64::consts::PI;
use vfvchem_vfvm2d::*;

/// Manufactured diffusion-advection-reaction problem.
/// Run with: `cargo run --release --example bench_heatflow`
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
/// dom - density (rho = 1 * (1 + 0.3 s))
/// dom - heat capacity (cp = 1000 * (1 + 0.2 t))
/// dom - velocity (prescribed field)
/// dom - thermal conductivity (k = 5e-3 * (1 + 0.5 s + 0.3 t))
/// dom - linear source (L = -1.0 * (1 + s t))
/// dom - heat source (manufactured solution)
///
/// Boundary Conditions:
/// bnd_l - temperature (manufactured solution)
/// bnd_r - temperature (manufactured solution)
/// bnd_b - temperature (manufactured solution)
/// bnd_t - temperature (manufactured solution)
/// 
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.1, 0.1, 0.005, 0.005, "examples/output_bench_heatflow")?;

    // geometry
    let dom = vars.add_dom("dom")?;  // domain
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left boundary
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right boundary
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom boundary
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top boundary

    // unknowns
    let temp = vars.add_scldom_unkn(dom, 300.0, "temp")?;  // temperature (T)

    // non-constant domain values
    let den_func: SclDomFunc = Box::new(|_t, (x, y), _unkn| den_at(x, y));
    let den = vars.add_scldom_func(dom, den_func, vec![], "den")?;  // density (rho)
    let hcap_func: SclDomFunc = Box::new(|_t, (x, y), _unkn| hcap_at(x, y));
    let hcap = vars.add_scldom_func(dom, hcap_func, vec![], "hcap")?;  // heat capacity (cp)
    let vel_func: VecDomFunc = Box::new(|_t, (x, y), _unkn| vel_at(x, y));
    let vel = vars.add_vecdom_func(dom, vel_func, vec![], "vel")?;  // velocity (vel)
    let cond_func: SclDomFunc = Box::new(|_t, (x, y), _unkn| cond_at(x, y));
    let cond = vars.add_scldom_func(dom, cond_func, vec![], "cond")?;  // thermal conductivity (k)
    let hlin_func: SclDomFunc = Box::new(|_t, (x, y), _unkn| lin_at(x, y));
    let hlin = vars.add_scldom_func(dom, hlin_func, vec![], "hlin")?;  // linear source (L)
    let hsrc_func: SclDomFunc = Box::new(|_t, (x, y), _unkn| src_at(x, y));
    let hsrc = vars.add_scldom_func(dom, hsrc_func, vec![], "hsrc")?;  // heat source (Q)

    // non-constant boundary values
    let temp_l_func: SclBndFunc = Box::new(|_t, (x, y), _unkn| exact(x, y));
    let temp_l = vars.add_sclbnd_func(bnd_l, temp_l_func, vec![])?;  // temperature
    let temp_r_func: SclBndFunc = Box::new(|_t, (x, y), _unkn| exact(x, y));
    let temp_r = vars.add_sclbnd_func(bnd_r, temp_r_func, vec![])?;  // temperature
    let temp_b_func: SclBndFunc = Box::new(|_t, (x, y), _unkn| exact(x, y));
    let temp_b = vars.add_sclbnd_func(bnd_b, temp_b_func, vec![])?;  // temperature
    let temp_t_func: SclBndFunc = Box::new(|_t, (x, y), _unkn| exact(x, y));
    let temp_t = vars.add_sclbnd_func(bnd_t, temp_t_func, vec![])?;  // temperature

    // physics
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom, temp, den, hcap, vel, cond, hlin, hsrc)?;  // heat flow domain: dom, T, rho, cp, vel, k, L, Q
    phys.add_temp_bnd(bnd_l, temp_l)?;  // temperature boundary: bnd, T
    phys.add_temp_bnd(bnd_r, temp_r)?;  // temperature boundary: bnd, T
    phys.add_temp_bnd(bnd_b, temp_b)?;  // temperature boundary: bnd, T
    phys.add_temp_bnd(bnd_t, temp_t)?;  // temperature boundary: bnd, T

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 1.0)?;  // physics solver

    Ok(())
}

fn phi(x: f64, y: f64) -> (f64, (f64, f64), f64) {
    let (s, t) = (x / 0.1, y / 0.1);
    let (sa, ca) = (1.5 * s + 0.3).sin_cos();
    let (sb, cb) = (1.2 * t - 0.4).sin_cos();
    let val = sa * cb;
    let grad = (1.5 / 0.1 * ca * cb, -1.2 / 0.1 * sa * sb);
    let lap = -(1.5 * 1.5 + 1.2 * 1.2) / (0.1 * 0.1) * val;
    (val, grad, lap)
}

fn vel_at(x: f64, y: f64) -> (f64, f64) {
    let (s, t) = (x / 0.1, y / 0.1);
    (1e-3 * (1.0 + 0.5 * t) * (1.0 + 0.3 * s), 1e-3 * 0.4 * (PI * s).sin() * (1.0 - t))
}

fn den_at(x: f64, _y: f64) -> f64 { 1.0 * (1.0 + 0.3 * x / 0.1) }
fn hcap_at(_x: f64, y: f64) -> f64 { 1000.0 * (1.0 + 0.2 * y / 0.1) }
fn cond_at(x: f64, y: f64) -> f64 { 5e-3 * (1.0 + 0.5 * x / 0.1 + 0.3 * y / 0.1) }
fn lin_at(x: f64, y: f64) -> f64 { -1.0 * (1.0 + x * y / (0.1 * 0.1)) }
fn exact(x: f64, y: f64) -> f64 { 300.0 + 50.0 * phi(x, y).0 }

fn src_at(x: f64, y: f64) -> f64 {
    let (val, (px, py), lap) = phi(x, y);
    let (ux, uy, ulap) = (50.0 * px, 50.0 * py, 50.0 * lap);
    let u = 300.0 + 50.0 * val;
    let a = den_at(x, y) * hcap_at(x, y);
    let (vx, vy) = vel_at(x, y);
    let g = cond_at(x, y);
    let (gx, gy) = (0.5 * 5e-3 / 0.1, 0.3 * 5e-3 / 0.1);
    a * (vx * ux + vy * uy) - g * ulap - (gx * ux + gy * uy) - lin_at(x, y) * u
}
