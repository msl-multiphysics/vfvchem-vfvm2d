use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct HeatDomAdvection {
    // domain
    pub dom_id: usize,

    // unknowns
    pub temp_id: usize,   // temperature (T)

    // others
    pub vel_id: usize,   // velocity (vel)
    pub cond_id: usize,  // thermal conductivity (k)
    pub den_id: usize,   // density (rho)
    pub hcap_id: usize,  // heat capacity (cp)
}

impl HeatDomAdvection {
    pub fn new(vars: &Variables, dom_id: usize, temp_id: usize, vel_id: usize, cond_id: usize, den_id: usize, hcap_id: usize) -> Result<HeatDomAdvection, VFVChemError> {
        // adds the advection term to the heat transfer equation
        // rho * cp * d(T)/dt + rho * cp * vel . grad(T) = -div( -k * grad(T) ) + L * T + R
        //                      ^^^^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, temp_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, cond_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, hcap_id)?;

        // create struct
        let mut oper = HeatDomAdvection::default();
        oper.dom_id = dom_id;
        oper.temp_id = temp_id;
        oper.vel_id = vel_id;
        oper.cond_id = cond_id;
        oper.den_id = den_id;
        oper.hcap_id = hcap_id;

        Ok(oper)
    }
}

impl OperatorBase for HeatDomAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the RHS of the PDE
        // -rho * cp * vel . grad(T)
        //
        // rewrite advection term
        // -rho * cp * vel . grad(T)
        // -div(rho * cp * vel * T) + T * div(rho * cp * vel)
        // -div(rho * cp * vel * T) + T * div(rho * cp * vel) + div(k * grad(T)) - div(k * grad(T))
        // -div(rho * cp * vel * T - k * grad(T)) + T * div(rho * cp * vel) - div(k * grad(T))
        //
        // get surface integral of 1st term
        // int_surf -div(rho * cp * vel * T - k * grad(T)) dA
        // int_lin (-rho * cp * vel * T + k * grad(T)) . n dl
        // sum_ngh (-rho * cp * vel * T + k * grad(T))_f . n dl
        //
        // apply power law scheme
        // let F = rho * cp * vel . n * dl; G = k * dl / dn; Pe = F / G; A = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (-F_f T_f + G_f * (T_n - T_v))
        // sum_ngh (-(G_f A_f + max(F_f, 0)) T_v + (G_f A_f + max(-F_f, 0)) T_n) -> 1st term
        //
        // get surface integral of 2nd term
        // int_surf T * div(rho * cp * vel) dA
        // T_c * div(rho * cp * vel)_c dA
        //
        // approximate divergence using divergence theorem
        // div(rho * cp * vel)_c = (1/dA) sum_ngh (rho * cp * vel)_f . n dl
        // T_c * sum_ngh (rho * cp * vel)_f . n dl
        // T_c * sum_ngh F_f
        // 
        // apply off-center correction
        // T_c = T_v + T_o 
        // (T_v + T_o) * sum_ngh F_f -> 2nd term
        //
        // get surface integral of 3rd term
        // int_surf -div(k * grad(T)) dA
        // int_lin -k * grad(T) . n dl
        // sum_ngh -k * (T_n - T_v) / dn * dl
        // sum_ngh -G_f * (T_n - T_v) -> 3rd term
        // 
        // add the terms
        // sum_ngh (-(G_f A_f + max(F_f, 0)) T_v + (G_f A_f + max(-F_f, 0)) T_n + (T_v + T_o) * F_f - G_f * (T_n - T_v))
        // sum_ngh ((G_f A_f + max(-F_f, 0) - G_f) * T_n + (-G_f A_f - max(F_f, 0) + F_f + G_f) * T_v + F_f * T_o)
        // sum_ngh ((G_f A_f + max(-F_f, 0) - G_f) * T_n - (G_f A_f + max(-F_f, 0) - G_f) * T_v + F_f * T_o)
        // sum_ngh ((G_f (A_f - 1) + max(-F_f, 0)) * T_n - (G_f (A_f - 1) + max(-F_f, 0)) * T_v + F_f * T_o)
        // sum_ngh ((G_f (A_f - 1) + max(-F_f, 0)) * (T_n - T_v) + F_f * T_o)

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let temp_dom = &vars.scldom[self.temp_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let cond_dom = &vars.scldom[self.cond_id];
        let den_dom = &vars.scldom[self.den_id];
        let hcap_dom = &vars.scldom[self.hcap_id];
        
        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let list_n = &tri.ngh_vid[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let offc = tri.vert_offc[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // get centroid properties
            let temp_o = temp_dom.correct_vert(offc, vid_v);

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let dist = dist_n[loc];
                let face = face_n[loc];
                let (nx, ny) = norm_n[loc];
                let skew = skew_n[loc];

                // get face properties
                let (velx_f, vely_f) = vel_dom.interp_face(skew, vid_v, vid_n);
                let den_f = den_dom.interp_face(skew, vid_v, vid_n);
                let cond_f = cond_dom.interp_face_harmonic(skew, vid_v, vid_n);
                let hcap_f = hcap_dom.interp_face(skew, vid_v, vid_n);
                
                // get power law term
                let diff = cond_f * face / dist;
                let flux = den_f * hcap_f * (velx_f * nx + vely_f * ny) * face;
                let peclet = flux / diff;
                let powlaw = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff_v = -diff * (powlaw - 1.0) - (-flux).max(0.0);
                self.add_sclscl_a(vars, a_triplet, self.temp_id, vid_v, self.temp_id, vid_n, -coeff_v)?;
                self.add_sclscl_a(vars, a_triplet, self.temp_id, vid_v, self.temp_id, vid_v,  coeff_v)?;
                self.add_scl_b(vars, b_vec, self.temp_id, vid_v, -flux * temp_o)?;
            }
        }

        Ok(())
    }
}
