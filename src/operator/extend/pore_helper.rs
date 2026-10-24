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

pub fn update_flow(vars: &mut Variables, sclflw_id: usize, epsl_id: usize, diam_id: usize) {
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
    let diam_dom = &vars.scldom[diam_id];

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
        let (kinv, fsk) = PoreOpVecDomDrag::closure(epsl_v, diam_dom.value[vid_v]);
        let speed = (u_v * u_v + v_v * v_v).sqrt();
        coeff_p += (visc_dom.value[vid_v] * kinv + 2.0 * den_v * fsk * speed) * tri.cvol_area[vid_v];
        rc_vec.push(if coeff_p > 0.0 { tri.cvol_area[vid_v] / coeff_p } else { 0.0 });
    }

    // get face mass flux
    flow_helper::update_flux(vars, sclflw_id, rc_vec);
}
