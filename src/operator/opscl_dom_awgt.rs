use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpSclDomAdvectionWeight {
    // domain
    pub dom_id: usize,

    // unknowns
    pub unk_id: usize,   // unknown scalar (c_i)

    // others
    pub vel_id: usize,  // velocity (vel)
    pub den_id: usize,  // density (rho)
    pub diff_id: usize,  // diffusion coefficient (D_ij)
    pub wgt_id: usize,  // weight (m)
}

impl OpSclDomAdvectionWeight {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, vel_id: usize, den_id: usize, diff_id: usize, wgt_id: usize) -> Result<OpSclDomAdvectionWeight, VFVChemError> {
        // adds the weighted advection term to the scalar transport equations
        // rho * m * d(c_i)/dt + rho * m * vel . grad(c_i) = -div( -D_ij * grad(c_j) ) + L_i c_i + R_i
        //                       ^^^^^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, diff_id)?;
        check_dom_scldom(vars, dom_id, wgt_id)?;

        // create struct
        let mut oper = OpSclDomAdvectionWeight::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.diff_id = diff_id;
        oper.wgt_id = wgt_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclDomAdvectionWeight {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the RHS of the PDE
        // -rho * m * vel . grad(c)
        //
        // rewrite advection term
        // -rho * m * vel . grad(c)
        // -div(rho * m * vel * c) + c * div(rho * m * vel)
        // -div(rho * m * vel * c) + c * div(rho * m * vel) + div(D * grad(c)) - div(D * grad(c))
        // -div(rho * m * vel * c - D * grad(c)) + c * div(rho * m * vel) - div(D * grad(c))
        //
        // get surface integral of 1st term
        // int_surf div(rho * m * vel * c - D * grad(c)) dA
        // int_lin (rho * m * vel * c - D * grad(c)) . n dl
        // sum_ngh (rho * m * vel * c - D * grad(c))_f . n dl
        //
        // apply power law scheme
        // let F = rho * m * vel . n * dl; G = D * dl / dn; Pe = F / G; A = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (F_f c_f - G_f * (c_n - c_v))
        // sum_ngh ((G_f A_f + max(F_f, 0)) c_v - (G_f A_f + max(-F_f, 0)) c_n) -> 1st term
        //
        // get surface integral of 2nd term
        // int_surf c * div(rho * m * vel) dA
        // c_c * div(rho * m * vel)_c dA
        //
        // approximate divergence using divergence theorem
        // div(rho * m * vel)_c = (1/dA) sum_ngh (rho * m * vel)_f . n dl
        // c_c * sum_ngh (rho * m * vel)_f . n dl
        // c_c * sum_ngh F_f
        // 
        // apply off-center correction
        // c_c = c_v + c_o 
        // (c_v + c_o) * sum_ngh F_f -> 2nd term
        //
        // get surface integral of 3rd term
        // int_surf div(D * grad(c)) dA
        // int_lin D * grad(c) . n dl
        // sum_ngh D * (c_n - c_v) / dn * dl
        // sum_ngh G_f * (c_n - c_v) -> 3rd term
        // 
        // add the terms
        // sum_ngh (-(G_f A_f + max(F_f, 0)) c_v + (G_f A_f + max(-F_f, 0)) c_n + (c_v + c_o) * F_f - G_f * (c_n - c_v))
        // sum_ngh ((G_f A_f + max(-F_f, 0) - G_f) * c_n + (-G_f A_f - max(F_f, 0) + F_f + G_f) * c_v + F_f * c_o)
        // sum_ngh ((G_f A_f + max(-F_f, 0) - G_f) * c_n - (G_f A_f + max(-F_f, 0) - G_f) * c_v + F_f * c_o)
        // sum_ngh ((G_f (A_f - 1) + max(-F_f, 0)) * c_n - (G_f (A_f - 1) + max(-F_f, 0)) * c_v + F_f * c_o)
        // sum_ngh ((G_f (A_f - 1) + max(-F_f, 0)) * (c_n - c_v) + F_f * c_o)

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        let diff_dom = &vars.scldom[self.diff_id];
        let wgt_dom = &vars.scldom[self.wgt_id];
        
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
            let unk_o = unk_dom.correct_vert(offc, vid_v);

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
                let diff_f = diff_dom.interp_face_harmonic(skew, vid_v, vid_n);
                let wgt_f = wgt_dom.interp_face(skew, vid_v, vid_n);
                
                // get power law term
                let diff = diff_f * face / dist;
                let flux = den_f * wgt_f * (velx_f * nx + vely_f * ny) * face;
                let peclet = flux / diff;
                let powlaw = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff_v = -diff * (powlaw - 1.0) - (-flux).max(0.0);
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_n, -coeff_v)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v,  coeff_v)?;
                self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -flux * unk_o)?;
            }
        }

        Ok(())
    }
}
