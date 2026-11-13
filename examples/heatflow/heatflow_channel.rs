use vfvchem_vfvm2d::*;

/// Steady-state heat and momentum transfer through a channel.
/// Run with: `cargo run --release --example heatflow_channel`
///
/// Geometry
/// --------
/// 
///        bnd_tl           bnd_tr
///       +---------------+-------+
/// bnd_i |               | dom_r |
///       +-------+ dom_f +-------+
///       | dom_l |               | bnd_o
///       +-------+---------------+
///        bnd_bl           bnd_br
/// 
/// Domains:
/// dom_f - flow channel
/// dom_l - left block
/// dom_r - right block
///
/// Boundaries:
/// bnd_i - inlet
/// bnd_o - outlet
/// bnd_bl - bottom left boundary
/// bnd_br - bottom right boundary
/// bnd_tl - top left boundary
/// bnd_tr - top right boundary
/// bnd_lc - left block to channel boundary
/// bnd_rc - right block to channel boundary
///
/// Interfaces:
/// int_lc - left block to channel interface
/// int_rc - right block to channel interface
/// 
/// Heat Transfer
/// -------------
/// 
/// Domain Properties:
/// dom_f - density (rho = 1.2 kg m-3)
/// dom_f - heat capacity (cp = 1005.0 J kg-1 K-1)
/// dom_f - thermal conductivity (k = 0.025 W m-1 K-1)
/// dom_f - heat source (Q = 0.0 W m-3)
/// dom_l - density (rho = 2700.0 kg m-3)
/// dom_l - heat capacity (cp = 890.0 J kg-1 K-1)
/// dom_l - thermal conductivity (k = 237.0 W m-1 K-1)
/// dom_l - heat source (Q = 2.0e3 W m-3)
/// dom_r - density (rho = 2700.0 kg m-3)
/// dom_r - heat capacity (cp = 890.0 J kg-1 K-1)
/// dom_r - thermal conductivity (k = 237.0 W m-1 K-1)
/// dom_r - heat source (Q = -2.0e3 W m-3)
///
/// Boundary Conditions:
/// bnd_i - temperature (T = 300.0 K)
/// bnd_o - outflow
/// bnd_bl - no-flux (q = 0.0 W m-2)
/// bnd_br - no-flux (q = 0.0 W m-2)
/// bnd_tl - no-flux (q = 0.0 W m-2)
/// bnd_tr - no-flux (q = 0.0 W m-2)
/// 
/// Interfaces Conditions:
/// int_lc - temperature and flux continuity
/// int_rc - temperature and flux continuity
/// 
/// Momentum Transfer
/// -----------------
/// 
/// Domain Properties:
/// dom_f - density (rho = 1.2 kg m-3)
/// dom_f - viscosity (mu = 1.8e-5 Pa s)
/// dom_f - body force (fce = <0.0, 0.0> N m-3)
///
/// Boundary Conditions:
/// bnd_i - velocity (parabolic; average velocity = 0.05 m s-1)
/// bnd_o - pressure (p = 0.0 Pa)
/// bnd_lc - no slip (vel = <0.0, 0.0> m s-1)
/// bnd_rc - no slip (vel = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_obj("examples/obj/obj_channel.obj", 0.001, 0.002, "examples/output_heatflow_channel")?;

    // geometry
    let dom_f = vars.add_dom("dom_f")?;  // flow channel
    let dom_l = vars.add_dom("dom_l")?;  // left block
    let dom_r = vars.add_dom("dom_r")?;  // right block
    let bnd_i = vars.add_bnd(dom_f, "bnd_i")?;    // inlet
    let bnd_o = vars.add_bnd(dom_f, "bnd_o")?;    // outlet
    let bnd_bl = vars.add_bnd(dom_l, "bnd_bl")?;  // bottom left boundary
    let bnd_br = vars.add_bnd(dom_f, "bnd_br")?;  // bottom right boundary
    let bnd_tl = vars.add_bnd(dom_f, "bnd_tl")?;  // top left boundary
    let bnd_tr = vars.add_bnd(dom_r, "bnd_tr")?;  // top right boundary
    let bnd_lc = vars.add_bnd(dom_f, "bnd_lc")?;  // left block to channel boundary
    let bnd_rc = vars.add_bnd(dom_f, "bnd_rc")?;  // right block to channel boundary
    let itf_lc = vars.add_itf(dom_l, dom_f, "itf_lc")?;  // left block to channel interface
    let itf_rc = vars.add_itf(dom_r, dom_f, "itf_rc")?;  // right block to channel interface

    // unknowns
    let temp_f = vars.add_scldom_unkn(dom_f, 0.0, "temp")?;       // temperature (T)
    let temp_l = vars.add_scldom_unkn(dom_l, 0.0, "temp")?;       // temperature (T)
    let temp_r = vars.add_scldom_unkn(dom_r, 0.0, "temp")?;       // temperature (T)
    let vel_f = vars.add_vecdom_unkn(dom_f, (0.0, 0.0), "vel")?;  // velocity (vel)
    let pres_f = vars.add_scldom_unkn(dom_f, 0.0, "pres")?;       // pressure (p)

    // domain values
    let den_f = vars.add_scldom_cnst(dom_f, 1.2, "den")?;         // density (rho)
    let hcap_f = vars.add_scldom_cnst(dom_f, 1005.0, "hcap")?;    // heat capacity (cp)
    let cond_f = vars.add_scldom_cnst(dom_f, 0.025, "cond")?;     // thermal conductivity (k)
    let hlin_f = vars.add_scldom_cnst(dom_f, 0.0, "hlin")?;       // linear source (L)
    let hsrc_f = vars.add_scldom_cnst(dom_f, 0.0, "hsrc")?;       // heat source (Q)
    let den_l = vars.add_scldom_cnst(dom_l, 2700.0, "den")?;      // density
    let hcap_l = vars.add_scldom_cnst(dom_l, 890.0, "hcap")?;     // heat capacity
    let vel_l = vars.add_vecdom_cnst(dom_l, (0.0, 0.0), "vel")?;  // velocity; zero for solids
    let cond_l = vars.add_scldom_cnst(dom_l, 237.0, "cond")?;     // thermal conductivity
    let hlin_l = vars.add_scldom_cnst(dom_l, 0.0, "hlin")?;       // linear source
    let hsrc_l = vars.add_scldom_cnst(dom_l, 2.0e3, "hsrc")?;     // heat source
    let den_r = vars.add_scldom_cnst(dom_r, 2700.0, "den")?;      // density
    let hcap_r = vars.add_scldom_cnst(dom_r, 890.0, "hcap")?;     // heat capacity
    let vel_r = vars.add_vecdom_cnst(dom_r, (0.0, 0.0), "vel")?;  // velocity; zero for solids
    let cond_r = vars.add_scldom_cnst(dom_r, 237.0, "cond")?;     // thermal conductivity
    let hlin_r = vars.add_scldom_cnst(dom_r, 0.0, "hlin")?;       // linear source
    let hsrc_r = vars.add_scldom_cnst(dom_r, -2.0e3, "hsrc")?;    // heat source
    let visc_f = vars.add_scldom_cnst(dom_f, 1.8e-5, "visc")?;    // viscosity (mu)
    let fce_f = vars.add_vecdom_cnst(dom_f, (0.0, 0.0), "fce")?;  // body force (fce)

    // non-constant boundary values
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y), _unkn| {  // t, (x, y), [unkn] -> (u, v)
        let u = 0.075 * 4.0 * (y - 0.05) * (0.10 - y) / (0.05 * 0.05);
        (u, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![])?;

    // constant boundary values
    let temp_i = vars.add_sclbnd_cnst(bnd_i, 300.0)?;        // temperature
    let hflx_bl = vars.add_sclbnd_cnst(bnd_bl, 0.0)?;        // heat flux (q)
    let hflx_br = vars.add_sclbnd_cnst(bnd_br, 0.0)?;        // heat flux
    let hflx_tl = vars.add_sclbnd_cnst(bnd_tl, 0.0)?;        // heat flux
    let hflx_tr = vars.add_sclbnd_cnst(bnd_tr, 0.0)?;        // heat flux
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0)?;          // pressure
    let vel_tl = vars.add_vecbnd_cnst(bnd_tl, (0.0, 0.0))?;  // velocity
    let vel_br = vars.add_vecbnd_cnst(bnd_br, (0.0, 0.0))?;  // velocity
    let vel_lc = vars.add_vecbnd_cnst(bnd_lc, (0.0, 0.0))?;  // velocity
    let vel_rc = vars.add_vecbnd_cnst(bnd_rc, (0.0, 0.0))?;  // velocity

    // physics
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom_f, temp_f, den_f, hcap_f, vel_f, cond_f, hlin_f, hsrc_f)?;  // heat transfer domain: dom, T, rho, cp, k, L, Q
    phys.add_heat_dom(dom_l, temp_l, den_l, hcap_l, vel_l, cond_l, hlin_l, hsrc_l)?;  // heat transfer domain: dom, T, rho, cp, k, L, Q
    phys.add_heat_dom(dom_r, temp_r, den_r, hcap_r, vel_r, cond_r, hlin_r, hsrc_r)?;  // heat transfer domain: dom, T, rho, cp, k, L, Q
    phys.add_temp_bnd(bnd_i, temp_i)?;    // temperature boundary: bnd, T
    phys.add_hflx_bnd(bnd_bl, hflx_bl)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_br, hflx_br)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_tl, hflx_tl)?;  // heat flux boundary: bnd, q
    phys.add_hflx_bnd(bnd_tr, hflx_tr)?;  // heat flux boundary: bnd, q
    phys.add_hout_bnd(bnd_o)?;                     // heat outflow boundary: bnd
    phys.add_hcnt_itf(itf_lc)?;                    // continuity interface: itf
    phys.add_hcnt_itf(itf_rc)?;                    // continuity interface: itf
    phys.add_flow_dom(dom_f, vel_f, pres_f, den_f, visc_f, fce_f)?;  // momentum transfer domain: dom, vel, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;     // velocity boundary: bnd, vel
    phys.add_pres_bnd(bnd_o, pres_o)?;  // pressure boundary: bnd, p
    phys.add_vel_bnd(bnd_tl, vel_tl)?;   // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_br, vel_br)?;   // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_lc, vel_lc)?;   // velocity boundary: bnd, vel
    phys.add_vel_bnd(bnd_rc, vel_rc)?;   // velocity boundary: bnd, vel

    // solver
    let mut lin = RslabLu::new(4)?;  // matrix solver
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 1.0)?;  // physics solver

    Ok(())
}
