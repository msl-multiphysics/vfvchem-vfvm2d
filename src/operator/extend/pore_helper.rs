use crate::base::prelude::*;
use crate::operator::core::flow_helper;
use crate::operator::extend::pore_opvec_dom_drag::PoreOpVecDomDrag;

pub fn visc_eff(vars: &Variables, visc_id: usize, epsl_id: usize, vid: usize) -> (f64, (f64, f64)) {
    // effective (Brinkman) viscosity and its gradient at a vertex from the liquid viscosity mu
    // mu_eff = mu / epsl
    // grad(mu_eff) = grad(mu) / epsl - mu * grad(epsl) / epsl^2

    // get vertex properties
    let visc_dom = &vars.scldom[visc_id];
    let visc = visc_dom.value[vid];
    let (visc_x, visc_y) = visc_dom.grad[vid];
    let epsl_dom = &vars.scldom[epsl_id];
    let epsl = epsl_dom.value[vid];
    let (epsl_x, epsl_y) = epsl_dom.grad[vid];

    (visc / epsl, (visc_x / epsl - visc * epsl_x / (epsl * epsl), visc_y / epsl - visc * epsl_y / (epsl * epsl)))
}

pub fn diff_eff(vars: &Variables, diff_id: usize, epsl_id: usize, tort_id: usize, vid: usize) -> (f64, (f64, f64)) {
    // effective diffusion coefficient and its gradient at a vertex from the molecular diffusion coefficient G
    // G_eff = (epsl / tau) * G
    // grad(G_eff) = (epsl / tau) * grad(G) + G * (grad(epsl) / tau - epsl * grad(tau) / tau^2)

    // get vertex properties
    let diff_dom = &vars.scldom[diff_id];
    let diff = diff_dom.value[vid];
    let (diff_x, diff_y) = diff_dom.grad[vid];
    let epsl_dom = &vars.scldom[epsl_id];
    let epsl = epsl_dom.value[vid];
    let (epsl_x, epsl_y) = epsl_dom.grad[vid];
    let tort_dom = &vars.scldom[tort_id];
    let tort = tort_dom.value[vid];
    let (tort_x, tort_y) = tort_dom.grad[vid];

    let fac = epsl / tort;
    let fac_x = epsl_x / tort - epsl * tort_x / (tort * tort);
    let fac_y = epsl_y / tort - epsl * tort_y / (tort * tort);
    (fac * diff, (fac * diff_x + diff * fac_x, fac * diff_y + diff * fac_y))
}

pub fn weighted(vars: &Variables, val_id: usize, wgt_id: usize, vid: usize) -> (f64, (f64, f64)) {
    // weighted scalar and its gradient at a vertex
    // w * s; grad(w * s) = w * grad(s) + s * grad(w)

    // get vertex properties
    let val_dom = &vars.scldom[val_id];
    let val = val_dom.value[vid];
    let (val_x, val_y) = val_dom.grad[vid];
    let wgt_dom = &vars.scldom[wgt_id];
    let wgt = wgt_dom.value[vid];
    let (wgt_x, wgt_y) = wgt_dom.grad[vid];

    (wgt * val, (wgt * val_x + val * wgt_x, wgt * val_y + val * wgt_y))
}

pub fn disp_face(vars: &Variables, vel_id: usize, epsl_id: usize, displ_id: usize, dispt_id: usize, vid_v: usize, vid_n: usize, norm: (f64, f64)) -> (f64, f64) {
    // normal-normal and normal-tangential components of the dispersion tensor at a face
    // D_disp = alpha_T * |v| * I + (alpha_L - alpha_T) * v v^T / |v|, v = vel / epsl
    // D_nn = alpha_T * |v| + (alpha_L - alpha_T) * (v . n)^2 / |v| >= 0
    // D_nt = (alpha_L - alpha_T) * (v . n) * (v . t) / |v|, t = (-n_y, n_x)
    // v and the dispersivities are arithmetic means of the vertex values

    // get vertex properties
    let vel_dom = &vars.vecdom[vel_id];
    let epsl_dom = &vars.scldom[epsl_id];
    let displ_dom = &vars.scldom[displ_id];
    let dispt_dom = &vars.scldom[dispt_id];
    let (u_v, v_v) = vel_dom.value[vid_v];
    let (u_n, v_n) = vel_dom.value[vid_n];
    let epsl_v = epsl_dom.value[vid_v];
    let epsl_n = epsl_dom.value[vid_n];

    // get face properties
    let displ = 0.5 * (displ_dom.value[vid_v] + displ_dom.value[vid_n]);
    let dispt = 0.5 * (dispt_dom.value[vid_v] + dispt_dom.value[vid_n]);
    if displ == 0.0 && dispt == 0.0 {
        return (0.0, 0.0);
    }
    let vx = 0.5 * (u_v / epsl_v + u_n / epsl_n);
    let vy = 0.5 * (v_v / epsl_v + v_n / epsl_n);
    let speed = (vx * vx + vy * vy).sqrt();
    if speed == 0.0 {
        return (0.0, 0.0);
    }

    // get components
    let (nx, ny) = norm;
    let vn = vx * nx + vy * ny;
    let vt = -vx * ny + vy * nx;
    (dispt * speed + (displ - dispt) * vn * vn / speed, (displ - dispt) * vn * vt / speed)
}

pub fn update_flow(vars: &mut Variables, sclflw_id: usize, epsl_id: usize, perm_id: usize, forch_id: usize) {
    // updates the flow data of a porous medium; lagged from the current velocity and pressure
    //
    // Rhie-Chow coefficient
    // a_P = sum_ngh D_f * A(|Pe_f|) + max(F_f, 0) / epsl_v^2 + (mu_v / K_v + 2 * rho_v * F_v / sqrt(K_v) * |vel_v|) * A
    // D_f uses mu_eff and Pe_f = F_f / (epsl_v * epsl_f * D_f); see PoreOpVecDomAdvection
    // the drag makes d = A / a_P tend to K / mu in the Darcy limit; see PoreOpVecDomDrag
    // F_f is the central mass flux

    // get objects
    let sclflw = &vars.sclflw[sclflw_id];
    let dom = &vars.dom[sclflw.dom_id];
    let tri = &vars.tri[sclflw.dom_id];
    let vel_dom = &vars.vecdom[sclflw.vel_id];
    let den_dom = &vars.scldom[sclflw.den_id];
    let visc_dom = &vars.scldom[sclflw.visc_id];
    let epsl_dom = &vars.scldom[epsl_id];
    let perm_dom = &vars.scldom[perm_id];
    let forch_dom = &vars.scldom[forch_id];

    // get A / a_P at every vertex
    let mut rc_vec = Vec::with_capacity(dom.num_tot());  // vid -> A / a_P
    for vid_v in 0..dom.num_tot() {
        // get geometry
        let list_n = &tri.ngh_vid[vid_v];
        let dist_n = &tri.ngh_dist[vid_v];
        let face_n = &tri.face_len[vid_v];
        let norm_n = &tri.face_norm[vid_v];
        let skew_n = &tri.face_skew[vid_v];

        // get vertex properties
        let (u_v, v_v) = vel_dom.value[vid_v];
        let den_v = den_dom.value[vid_v];
        let (den_x_v, den_y_v) = den_dom.grad[vid_v];
        let (visc_v, (visc_x_v, visc_y_v)) = visc_eff(vars, sclflw.visc_id, epsl_id, vid_v);
        let epsl_v = epsl_dom.value[vid_v];

        // iterate through faces
        let mut coeff_p = 0.0;
        for (loc, &vid_n) in list_n.iter().enumerate() {
            // get geometry
            let face = face_n[loc];
            let dist = dist_n[loc];
            let (nx, ny) = norm_n[loc];
            let (skew_x, skew_y) = skew_n[loc];

            // get neighbor properties
            let (u_n, v_n) = vel_dom.value[vid_n];
            let den_n = den_dom.value[vid_n];
            let (den_x_n, den_y_n) = den_dom.grad[vid_n];
            let (visc_n, (visc_x_n, visc_y_n)) = visc_eff(vars, sclflw.visc_id, epsl_id, vid_n);
            let epsl_n = epsl_dom.value[vid_n];

            // get face density, viscosity, and inverse porosity
            let den_f = 0.5 * (den_v + den_n) + 0.5 * (skew_x * (den_x_v + den_x_n) + skew_y * (den_y_v + den_y_n));  // skew-corrected arithmetic mean
            let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
            let visc_f = visc_m + 0.5 * (skew_x * (visc_x_v + visc_x_n) + skew_y * (visc_y_v + visc_y_n));  // skew correction
            let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity
            let ieps_f = 0.5 * (1.0 / epsl_v + 1.0 / epsl_n);  // arithmetic mean of 1 / epsl

            // get conductance, mass flux, and Peclet number
            let cond = visc_f * face / dist;
            let flux = den_f * 0.5 * ((u_v + u_n) * nx + (v_v + v_n) * ny) * face;
            let peclet = if cond != 0.0 { ieps_f * flux / (epsl_v * cond) } else { 0.0 };

            // D_f * A(|Pe_f|) + max(F_f, 0) / epsl_v^2
            coeff_p += cond * (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5) + flux.max(0.0) / (epsl_v * epsl_v);
        }

        // implicit Darcy and linearized Forchheimer drag; same diagonal as PoreOpVecDomDrag
        let (kinv, fsk) = PoreOpVecDomDrag::coeff(perm_dom.value[vid_v], forch_dom.value[vid_v]);
        let speed = (u_v * u_v + v_v * v_v).sqrt();
        coeff_p += (visc_dom.value[vid_v] * kinv + 2.0 * den_v * fsk * speed) * tri.cvol_area[vid_v];
        rc_vec.push(if coeff_p > 0.0 { tri.cvol_area[vid_v] / coeff_p } else { 0.0 });
    }

    // get face mass flux
    flow_helper::update_flux(vars, sclflw_id, rc_vec);
}
