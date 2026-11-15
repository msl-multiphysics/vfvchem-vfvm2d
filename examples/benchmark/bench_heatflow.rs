use std::f64::consts::PI;
use vfvchem_vfvm2d::*;

/// Manufactured advection-diffusion-reaction solution for SteadyHeatFlow.
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
/// The velocity is prescribed (no flow domain) and is not solenoidal.
/// All properties vary in space. Temperature and source follow a manufactured solution.
///
/// Domain Properties:
/// dom - density (rho = 1 * (1 + 0.3 s) kg m-3)
/// dom - heat capacity (cp = 1000 * (1 + 0.2 t) J kg-1 K-1)
/// dom - thermal conductivity (k = 5e-3 * (1 + 0.5 s + 0.3 t) W m-1 K-1)
/// dom - linear source (L = -1.0 * (1 + s t) W m-3 K-1)
/// dom - heat source (Q from manufactured solution)
/// dom - velocity (prescribed non-solenoidal field)
///
/// Boundary Conditions:
/// bnd_l - temperature (manufactured solution)
/// bnd_r - temperature (manufactured solution)
/// bnd_b - temperature (manufactured solution)
/// bnd_t - temperature (manufactured solution)
///
const LEN: f64 = 0.1;
const VEL: f64 = 1e-3;
// phi = sin(c0 s + c1) * cos(c2 t + c3); smooth maximum inside the domain
const CF: [f64; 4] = [1.5, 0.3, 1.2, -0.4];
const U0: f64 = 300.0;
const DU: f64 = 50.0;
const RHO_0: f64 = 1.0;
const CP_0: f64 = 1000.0;
const K_0: f64 = 5e-3;
const LIN_0: f64 = -1.0;

fn phi(x: f64, y: f64) -> (f64, (f64, f64), f64) {
    let (s, t) = (x / LEN, y / LEN);
    let (sa, ca) = (CF[0] * s + CF[1]).sin_cos();
    let (sb, cb) = (CF[2] * t + CF[3]).sin_cos();
    let val = sa * cb;
    let grad = (CF[0] / LEN * ca * cb, -CF[2] / LEN * sa * sb);
    let lap = -(CF[0] * CF[0] + CF[2] * CF[2]) / (LEN * LEN) * val;
    (val, grad, lap)
}

fn vel_at(x: f64, y: f64) -> (f64, f64) {
    let (s, t) = (x / LEN, y / LEN);
    (VEL * (1.0 + 0.5 * t) * (1.0 + 0.3 * s), VEL * 0.4 * (PI * s).sin() * (1.0 - t))
}

fn den_at(x: f64, _y: f64) -> f64 { RHO_0 * (1.0 + 0.3 * x / LEN) }
fn hcap_at(_x: f64, y: f64) -> f64 { CP_0 * (1.0 + 0.2 * y / LEN) }
fn cond_at(x: f64, y: f64) -> f64 { K_0 * (1.0 + 0.5 * x / LEN + 0.3 * y / LEN) }
fn lin_at(x: f64, y: f64) -> f64 { LIN_0 * (1.0 + x * y / (LEN * LEN)) }
fn exact(x: f64, y: f64) -> f64 { U0 + DU * phi(x, y).0 }

fn src_at(x: f64, y: f64) -> f64 {
    let (val, (px, py), lap) = phi(x, y);
    let (ux, uy, ulap) = (DU * px, DU * py, DU * lap);
    let u = U0 + DU * val;
    let a = den_at(x, y) * hcap_at(x, y);
    let (vx, vy) = vel_at(x, y);
    let g = cond_at(x, y);
    let (gx, gy) = (0.5 * K_0 / LEN, 0.3 * K_0 / LEN);
    a * (vx * ux + vy * uy) - g * ulap - (gx * ux + gy * uy) - lin_at(x, y) * u
}

fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, LEN, LEN, 0.005, 0.005, "examples/output_bench_heatflow")?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let temp = vars.add_scldom_unkn(dom, U0, "temp")?;

    // domain values
    let den = vars.add_scldom_func(dom, Box::new(|_t, (x, y), _u| den_at(x, y)), vec![], "den")?;
    let hcap = vars.add_scldom_func(dom, Box::new(|_t, (x, y), _u| hcap_at(x, y)), vec![], "hcap")?;
    let vel = vars.add_vecdom_func(dom, Box::new(|_t, (x, y), _u| vel_at(x, y)), vec![], "vel")?;
    let cond = vars.add_scldom_func(dom, Box::new(|_t, (x, y), _u| cond_at(x, y)), vec![], "cond")?;
    let hlin = vars.add_scldom_func(dom, Box::new(|_t, (x, y), _u| lin_at(x, y)), vec![], "hlin")?;
    let hsrc = vars.add_scldom_func(dom, Box::new(|_t, (x, y), _u| src_at(x, y)), vec![], "hsrc")?;

    // boundary values
    let temp_l = vars.add_sclbnd_func(bnd_l, Box::new(|_t, (x, y), _u| exact(x, y)), vec![])?;
    let temp_r = vars.add_sclbnd_func(bnd_r, Box::new(|_t, (x, y), _u| exact(x, y)), vec![])?;
    let temp_b = vars.add_sclbnd_func(bnd_b, Box::new(|_t, (x, y), _u| exact(x, y)), vec![])?;
    let temp_t = vars.add_sclbnd_func(bnd_t, Box::new(|_t, (x, y), _u| exact(x, y)), vec![])?;

    // physics
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom, temp, den, hcap, vel, cond, hlin, hsrc)?;
    phys.add_temp_bnd(bnd_l, temp_l)?;
    phys.add_temp_bnd(bnd_r, temp_r)?;
    phys.add_temp_bnd(bnd_b, temp_b)?;
    phys.add_temp_bnd(bnd_t, temp_t)?;

    // solver
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 1.0)?;

    Ok(())
}
