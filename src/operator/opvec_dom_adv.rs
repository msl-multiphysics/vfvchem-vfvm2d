use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpVecDomAdvection {
    // domain
    pub dom_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub den_id: usize,   // density (rho)
    pub visc_id: usize,  // viscosity (mu)
}

impl OpVecDomAdvection {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, den_id: usize, visc_id: usize) -> Result<OpVecDomAdvection, VFVChemError> {
        // adds the advection term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                   ^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;

        // create struct
        let mut oper = OpVecDomAdvection::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecDomAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the RHS of the PDE
        // -rho * vel . grad(vel)
        //
        // rewrite advection term
        // -rho * vel . grad(vel)
        // -div(rho * vel * vel) + vel * div(rho * vel)
        // -div(rho * vel * vel) + vel * div(rho * vel) + div(mu * grad(vel)) - div(mu * grad(vel))
        // -div(rho * vel * vel - mu * grad(vel)) + vel * div(rho * vel) - div(mu * grad(vel))
        //
        // get surface integral of 1st term
        // int_surf div(rho * vel * vel - mu * grad(vel)) dA
        // int_lin (rho * vel * vel - mu * grad(vel)) . n dl
        // sum_ngh (rho * vel * vel - mu * grad(vel))_f . n dl
        //
        // apply power law scheme
        // let F = rho * vel . n * dl; G = mu * dl / dn; Pe = F / G; A = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (F_f vel_f - G_f * (vel_n - vel_v))
        // sum_ngh ((G_f A_f + max(F_f, 0)) vel_v - (G_f A_f + max(-F_f, 0)) vel_n) -> 1st term
        //
        // get surface integral of 2nd term
        // int_surf vel * div(rho * vel) dA
        // vel_c * div(rho * vel)_c dA
        //
        // approximate divergence using divergence theorem
        // div(rho * vel)_c = (1/dA) sum_ngh (rho * vel)_f . n dl
        // vel_c * sum_ngh (rho * vel)_f . n dl
        // vel_c * sum_ngh F_f
        // 
        // apply off-center correction
        // vel_c = vel_v + vel_o 
        // (vel_v + vel_o) * sum_ngh F_f -> 2nd term
        //
        // get surface integral of 3rd term
        // int_surf div(mu * grad(vel)) dA
        // int_lin mu * grad(vel) . n dl
        // sum_ngh mu * (vel_n - vel_v) / dn * dl
        // sum_ngh G_f * (vel_n - vel_v) -> 3rd term
        // 
        // add the terms
        // sum_ngh (-(G_f A_f + max(F_f, 0)) vel_v + (G_f A_f + max(-F_f, 0)) vel_n + (vel_v + vel_o) * F_f - G_f * (vel_n - vel_v))
        // sum_ngh ((G_f A_f + max(-F_f, 0) - G_f) * vel_n + (-G_f A_f - max(F_f, 0) + F_f + G_f) * vel_v + F_f * vel_o)
        // sum_ngh ((G_f A_f + max(-F_f, 0) - G_f) * vel_n - (G_f A_f + max(-F_f, 0) - G_f) * vel_v + F_f * vel_o)
        // sum_ngh ((G_f (A_f - 1) + max(-F_f, 0)) * vel_n - (G_f (A_f - 1) + max(-F_f, 0)) * vel_v + F_f * vel_o)
        // sum_ngh ((G_f (A_f - 1) + max(-F_f, 0)) * (vel_n - vel_v) + F_f * vel_o)

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        let visc_dom = &vars.scldom[self.visc_id];
        
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
            let (velx_o, vely_o) = vel_dom.correct_vert(offc, vid_v);

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
                let visc_f = visc_dom.interp_face_harmonic(skew, vid_v, vid_n);
                
                // get power law term
                let diff = visc_f * face / dist;
                let flux = den_f * (velx_f * nx + vely_f * ny) * face;
                let peclet = flux / diff;
                let powlaw = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff_v = -diff * (powlaw - 1.0) - (-flux).max(0.0);
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_n, -coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v,  coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_n, -coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v,  coeff_v)?;
                self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -flux * velx_o)?;
                self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -flux * vely_o)?;
            }
        }

        Ok(())
    }
}
