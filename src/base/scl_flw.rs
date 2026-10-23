use crate::base::prelude::*;

#[derive(Default)]
pub struct ScalarFlow {
    // ids
    pub sclflw_id: usize,
    pub dom_id: usize,  // domain this is attached to

    // vectors
    pub vel_id: usize,  // velocity (vel)
    pub fce_id: usize,  // body force (fce)

    // scalars
    pub pres_id: usize,  // pressure (p)
    pub den_id: usize,   // density (rho)
    pub visc_id: usize,  // viscosity (mu)

    // flow data
    // shared by the momentum and continuity operators so that both use the same face mass flux
    pub rc_coeff: Vec<f64>,   // vid -> A / a_P for the Rhie-Chow correction
    pub flux: Vec<Vec<f64>>,  // vid, loc -> face mass flux
}

impl ScalarFlow {
    pub fn new(sclflw_id: usize, dom: &Domain, vel_id: usize, pres_id: usize, den_id: usize, visc_id: usize, fce_id: usize) -> Result<ScalarFlow, VFVChemError> {
        // default class
        let mut sclflw = ScalarFlow::default();
        sclflw.sclflw_id = sclflw_id;
        sclflw.dom_id = dom.dom_id;

        // field ids
        sclflw.vel_id = vel_id;
        sclflw.fce_id = fce_id;
        sclflw.pres_id = pres_id;
        sclflw.den_id = den_id;
        sclflw.visc_id = visc_id;

        // flow data is computed in update_flow after the triangulation
        sclflw.rc_coeff = Vec::new();
        sclflw.flux = Vec::new();

        Ok(sclflw)
    }

    pub fn update_flow(vars: &mut Variables, sclflw_id: usize) {
        // does not use &mut self to prevent borrowing problems
        // lagged from the current velocity and pressure
        //
        // Rhie-Chow coefficient
        // a_P = sum_ngh D_f * A(|Pe_f|) + max(F_f, 0); see OpVecDomAdvection
        // F_f is the central mass flux
        //
        // face mass flux of the continuity and momentum equations
        // F_f = rho_f * vel_f . n * dl
        // vel_f . n = ((vel_v + vel_n) / 2 + (grad(vel_v) + grad(vel_n)) / 2 . skew) . n - d_f * ((p_n - p_v) / dn - fce_f . n) + (d_v * (grad(p_v) - fce_v) + d_n * (grad(p_n) - fce_n)) / 2 . n
        // see OpSclDomDivergence, OpSclDomPressureStabilizer, and OpVecDomAdvection

        // get objects
        let sclflw = &vars.sclflw[sclflw_id];
        let dom = &vars.dom[sclflw.dom_id];
        let tri = &vars.tri[sclflw.dom_id];
        let vel_dom = &vars.vecdom[sclflw.vel_id];
        let fce_dom = &vars.vecdom[sclflw.fce_id];
        let pres_dom = &vars.scldom[sclflw.pres_id];
        let den_dom = &vars.scldom[sclflw.den_id];
        let visc_dom = &vars.scldom[sclflw.visc_id];

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
            let visc_v = visc_dom.value[vid_v];
            let (visc_x_v, visc_y_v) = visc_dom.grad[vid_v];

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
                let visc_n = visc_dom.value[vid_n];
                let (visc_x_n, visc_y_n) = visc_dom.grad[vid_n];

                // get face density and viscosity
                let den_f = 0.5 * (den_v + den_n) + 0.5 * (skew_x * (den_x_v + den_x_n) + skew_y * (den_y_v + den_y_n));  // skew-corrected arithmetic mean
                let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                let visc_f = visc_m + 0.5 * (skew_x * (visc_x_v + visc_x_n) + skew_y * (visc_y_v + visc_y_n));  // skew correction
                let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity

                // get conductance, mass flux, and Peclet number
                let cond = visc_f * face / dist;
                let flux = den_f * 0.5 * ((u_v + u_n) * nx + (v_v + v_n) * ny) * face;
                let peclet = if cond != 0.0 { flux / cond } else { 0.0 };

                // D_f * A(|Pe_f|) + max(F_f, 0)
                coeff_p += cond * (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5) + flux.max(0.0);
            }
            rc_vec.push(if coeff_p > 0.0 { tri.cvol_area[vid_v] / coeff_p } else { 0.0 });
        }

        // get d * (grad(p) - fce) of the momentum equation at every vertex
        // grad(p_v) = sum_m (c_x, c_y) * p_m / A_v; see Triangle
        // fce_v = fce at the control volume centroid; see OpVecDomSource
        let mut mom_vec = Vec::with_capacity(dom.num_tot());  // vid -> d * (grad(p) - fce)
        for vid_v in 0..dom.num_tot() {
            // get geometry
            let (cx, cy) = tri.cvol_skew[vid_v];
            let area = tri.cvol_area[vid_v];

            // get pressure gradient
            let mut pres_x = 0.0;
            let mut pres_y = 0.0;
            for &(vid_m, c_x, c_y) in &tri.grad_coeff[vid_v] {
                pres_x += c_x * pres_dom.value[vid_m] / area;
                pres_y += c_y * pres_dom.value[vid_m] / area;
            }

            // get body force at centroid
            let (fce_x, fce_y) = fce_dom.value[vid_v];
            let (fce_x_x, fce_x_y, fce_y_x, fce_y_y) = fce_dom.grad[vid_v];
            let fce_x_c = fce_x + fce_x_x * cx + fce_x_y * cy;
            let fce_y_c = fce_y + fce_y_x * cx + fce_y_y * cy;

            mom_vec.push((rc_vec[vid_v] * (pres_x - fce_x_c), rc_vec[vid_v] * (pres_y - fce_y_c)));
        }

        // get face mass flux at every face
        let mut flux_vec = Vec::with_capacity(dom.num_tot());  // vid, loc -> F_f
        for vid_v in 0..dom.num_tot() {
            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // get vertex properties
            let (u_v, v_v) = vel_dom.value[vid_v];
            let (u_x_v, u_y_v, v_x_v, v_y_v) = vel_dom.grad[vid_v];
            let pres_v = pres_dom.value[vid_v];
            let den_v = den_dom.value[vid_v];
            let (den_x_v, den_y_v) = den_dom.grad[vid_v];
            let (fce_x_v, fce_y_v) = fce_dom.value[vid_v];
            let rc_v = rc_vec[vid_v];
            let (mom_x_v, mom_y_v) = mom_vec[vid_v];

            // iterate through faces
            let mut flux_n = Vec::with_capacity(list_n.len());
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (nx, ny) = norm_n[loc];
                let (skew_x, skew_y) = skew_n[loc];

                // get neighbor properties
                let (u_n, v_n) = vel_dom.value[vid_n];
                let (u_x_n, u_y_n, v_x_n, v_y_n) = vel_dom.grad[vid_n];
                let pres_n = pres_dom.value[vid_n];
                let den_n = den_dom.value[vid_n];
                let (den_x_n, den_y_n) = den_dom.grad[vid_n];
                let (fce_x_n, fce_y_n) = fce_dom.value[vid_n];
                let rc_n = rc_vec[vid_n];
                let (mom_x_n, mom_y_n) = mom_vec[vid_n];

                // get face density and velocity
                let den_f = 0.5 * (den_v + den_n) + 0.5 * (skew_x * (den_x_v + den_x_n) + skew_y * (den_y_v + den_y_n));  // skew-corrected arithmetic mean
                let u_f = 0.5 * (u_v + u_n) + 0.5 * (skew_x * (u_x_v + u_x_n) + skew_y * (u_y_v + u_y_n));  // skew-corrected arithmetic mean
                let v_f = 0.5 * (v_v + v_n) + 0.5 * (skew_x * (v_x_v + v_x_n) + skew_y * (v_y_v + v_y_n));  // skew-corrected arithmetic mean

                // get Rhie-Chow correction
                // fce_f is the arithmetic mean at the edge midpoint, consistent with (p_n - p_v) / dn
                let fce_e = 0.5 * (nx * (fce_x_v + fce_x_n) + ny * (fce_y_v + fce_y_n));
                let mom_e = 0.5 * (nx * (mom_x_v + mom_x_n) + ny * (mom_y_v + mom_y_n));
                let vel_rc = -0.5 * (rc_v + rc_n) * ((pres_n - pres_v) / dist - fce_e) + mom_e;

                // rho_f * vel_f . n * dl
                flux_n.push(den_f * (u_f * nx + v_f * ny + vel_rc) * face);
            }
            flux_vec.push(flux_n);
        }

        // modify sclflw values
        let sclflw = &mut vars.sclflw[sclflw_id];
        sclflw.rc_coeff = rc_vec;
        sclflw.flux = flux_vec;
    }

}
