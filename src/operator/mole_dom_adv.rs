use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct MoleDomAdvection {
    // domain
    pub dom_id: usize,

    // unknowns
    pub conc_id: usize,   // concentration (c_i); mole solute per mass solvent

    // others
    pub vel_id: usize,  // velocity (vel)
    pub den_id: usize,  // density (rho)
    pub diff_id: usize,  // diffusion coefficient (D_ij)
}

impl MoleDomAdvection {
    pub fn new(vars: &Variables, dom_id: usize, conc_id: usize, vel_id: usize, den_id: usize, diff_id: usize) -> Result<MoleDomAdvection, VFVChemError> {
        // adds the advection term to the mass transfer equation
        // rho * d(c_i)/dt + rho * vel . grad(c_i) = -div( -rho * D_ij * grad(c_j) ) + L_i * c_i + R_i
        //                   ^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, conc_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, diff_id)?;

        // create struct
        let mut oper = MoleDomAdvection::default();
        oper.dom_id = dom_id;
        oper.conc_id = conc_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.diff_id = diff_id;

        Ok(oper)
    }
}

impl OperatorBase for MoleDomAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the RHS of the PDE
        // -rho * vel . grad(c)
        //
        // rewrite advection term
        // -rho * vel . grad(c)
        // -div(rho * vel * c) + c * div(rho * vel)
        // -div(rho * vel * c) + c * div(rho * vel) + div(rho * D * grad(c)) - div(rho * D * grad(c))
        // -div(rho * vel * c - rho * D * grad(c)) + c * div(rho * vel) - div(rho * D * grad(c))
        //
        // get surface integral of 1st term
        // int_surf div(rho * vel * c - rho * D * grad(c)) dA
        // int_lin (rho * vel * c - rho * D * grad(c)) . n dl
        // sum_ngh (rho * vel * c - rho * D * grad(c))_f . n dl
        //
        // apply power law scheme
        // let F = rho * vel . n * dl; G = rho * D * dl / dn; Pe = F / G; A = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (F_f c_f - G_f * (c_n - c_v))
        // sum_ngh ((G_f A_f + max(F_f, 0)) c_v - (G_f A_f + max(-F_f, 0)) c_n) -> 1st term
        //
        // get surface integral of 2nd term
        // int_surf c * div(rho * vel) dA
        // c_c * div(rho * vel)_c dA
        //
        // approximate divergence using divergence theorem
        // div(rho * vel)_c = (1/dA) sum_ngh (rho * vel)_f . n dl
        // c_c * sum_ngh (rho * vel)_f . n dl
        // c_c * sum_ngh F_f
        // 
        // apply off-center correction
        // c_c = c_v + c_o 
        // (c_v + c_o) * sum_ngh F_f -> 2nd term
        //
        // get surface integral of 3rd term
        // int_surf div(rho * D * grad(c)) dA
        // int_lin rho * D * grad(c) . n dl
        // sum_ngh rho * D * (c_n - c_v) / dn * dl
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
        let conc_dom = &vars.scldom[self.conc_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        let diff_dom = &vars.scldom[self.diff_id];
        
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
            let conc_o = conc_dom.correct_vert(offc, vid_v);

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
                
                // get power law term
                let diff = den_f * diff_f * face / dist;
                let flux = den_f * (velx_f * nx + vely_f * ny) * face;
                let peclet = flux / diff;
                let powlaw = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff_v = -diff * (powlaw - 1.0) - (-flux).max(0.0);
                self.add_sclscl_a(vars, a_triplet, self.conc_id, vid_v, self.conc_id, vid_n, -coeff_v)?;
                self.add_sclscl_a(vars, a_triplet, self.conc_id, vid_v, self.conc_id, vid_v,  coeff_v)?;
                self.add_scl_b(vars, b_vec, self.conc_id, vid_v, -flux * conc_o)?;
            }
        }

        Ok(())
    }
}
