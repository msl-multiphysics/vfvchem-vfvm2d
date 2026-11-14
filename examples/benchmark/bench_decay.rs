use vfvchem_vfvm2d::*;

/// Adaptive time stepping benchmark with nonlinear decay.
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
/// dom - density (rho = 1.0 + 0.5 * x + 0.25 * y kg m-3)
/// dom - heat capacity (cp = 2.0 J kg-1 K-1)
/// dom - thermal conductivity (k = 0.0 W m-1 K-1)
/// dom - linear source (L = -rho * cp * T / (1.0 + 0.5 * x) W m-3 K-1)
///
/// Initial Conditions:
/// dom - temperature (T = 1.0 + 0.5 * x K)
/// 
/// Boundary Conditions:
/// bnd_l - no flux (q = 0.0 W m-2)
/// bnd_r - no flux (q = 0.0 W m-2)
/// bnd_b - no flux (q = 0.0 W m-2)
/// bnd_t - no flux (q = 0.0 W m-2)
///
/// Exact Solution
/// --------------
///
/// rho * cp * dT/dt = L * T = -rho * cp * T^2 / (1.0 + 0.5 * x)
/// T = (1.0 + 0.5 * x) / (1.0 + t)
/// the linear profile is reproduced exactly in space, so the error is purely temporal
/// L depends on T and is updated in the nonlinear iterations
/// the large initial step is rejected until the local error is within tolerance
/// with local error control, bdf2 gives a global error of about tol_time^(2/3)
///
fn main() -> Result<(), VFVChemError> {
    // refinement of the nonlinear tolerance; solve_adaptive uses tol_time = 100 * tol
    let tol_list = [1e-5, 1e-6, 1e-7, 1e-8];
    let t_end: f64 = 10.0;
    let mut err_list = Vec::new();  // (L1, L2, Linf) at t_end

    for (i, &tol) in tol_list.iter().enumerate() {
        // variables
        let folder = format!("examples/output_bench_decay_{i}");
        let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, 0.1, 0.1, &folder)?;

        // geometry
        let dom = vars.add_dom("dom")?;  // domain
        let bnd_l = vars.add_bnd(dom, "bnd_l")?;  // left boundary
        let bnd_r = vars.add_bnd(dom, "bnd_r")?;  // right boundary
        let bnd_b = vars.add_bnd(dom, "bnd_b")?;  // bottom boundary
        let bnd_t = vars.add_bnd(dom, "bnd_t")?;  // top boundary

        // unknowns
        let temp = vars.add_scldom_unkn(dom, 1.0, "temp")?;  // temperature (T)

        // set linear initial temperature
        // vertices are ordered as points, then nodes
        let pos: Vec<(f64, f64)> = vars.dom[dom].pts_gpt.iter().map(|&gpt| vars.pts.pts_pos[gpt])
            .chain(vars.dom[dom].nds_gnd.iter().map(|&gnd| vars.nds.nds_pos[gnd])).collect();
        for (vid, &(x, _y)) in pos.iter().enumerate() {
            vars.scldom[temp].value[vid] = 1.0 + 0.5 * x;
        }

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
        phys.solve_adaptive(&mut vars, &mut lin, 50, tol, 1.0, 0.5, 1e-8, 2.0, t_end, 10)?;  // physics solver

        // area-weighted errors at t_end
        let area = &vars.tri[dom].vert_area;
        let mut sum_area = 0.0;
        let mut sum_l1 = 0.0;
        let mut sum_l2 = 0.0;
        let mut err_linf: f64 = 0.0;
        for (vid, &(x, _y)) in pos.iter().enumerate() {
            let err = (vars.scldom[temp].value[vid] - (1.0 + 0.5 * x) / (1.0 + t_end)).abs();
            sum_area += area[vid];
            sum_l1 += area[vid] * err;
            sum_l2 += area[vid] * err * err;
            err_linf = err_linf.max(err);
        }
        err_list.push((sum_l1 / sum_area, (sum_l2 / sum_area).sqrt(), err_linf));
    }

    // print errors and observed order with respect to tol_time
    println!("Nonlinear decay: errors at t = {t_end}");
    println!("tol_time    L1          L2          Linf        order");
    for (i, &(err_l1, err_l2, err_linf)) in err_list.iter().enumerate() {
        let order = if i == 0 { f64::NAN } else { (err_linf / err_list[i - 1].2).ln() / (tol_list[i] / tol_list[i - 1]).ln() };
        println!("{:.2e}    {err_l1:.4e}  {err_l2:.4e}  {err_linf:.4e}  {order:.3}", 100.0 * tol_list[i]);
    }

    Ok(())
}
