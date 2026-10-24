use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclDomDispersion {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,    // unknown scalar (c_i)
    pub den_id: usize,    // density (rho)
    pub epsl_id: usize,   // porosity (epsl)
    pub displ_id: usize,  // longitudinal dispersivity (alpha_L)
    pub dispt_id: usize,  // transverse dispersivity (alpha_T)

    // vectors
    pub vel_id: usize,  // superficial velocity (vel)
}

impl PoreOpSclDomDispersion {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, den_id: usize, vel_id: usize, epsl_id: usize, displ_id: usize, dispt_id: usize) -> Result<PoreOpSclDomDispersion, VFVChemError> {
        // adds the mechanical dispersion term of a porous medium
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( -G_eff,ij * grad(c_j) - rho * D_disp * grad(c_i) ) + epsl * R_v,i + a_v * R_s,i
        //                                                                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // D_disp = alpha_T * |v| * I + (alpha_L - alpha_T) * v v^T / |v| with the interstitial velocity v = vel / epsl
        // D_disp is the same for all components; alpha_L = alpha_T = 0 turns it off

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;
        check_dom_scldom(vars, dom_id, displ_id)?;
        check_dom_scldom(vars, dom_id, dispt_id)?;

        // create struct
        let mut oper = PoreOpSclDomDispersion::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.den_id = den_id;
        oper.vel_id = vel_id;
        oper.epsl_id = epsl_id;
        oper.displ_id = displ_id;
        oper.dispt_id = dispt_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclDomDispersion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -div( -rho * D_disp * grad(c) ) dA
        // sum_ngh rho_f * (n . D_disp,f . grad(c)) * dl
        // sum_ngh rho_f * (D_nn,f * (c_n - c_v) / dn + D_nt,f * grad(c)_f . t) * dl
        //
        // t = (-n_y, n_x) is the face tangent
        // D_nn >= 0 is implicit and keeps the M-matrix; the cross term D_nt is deferred (lagged gradient)
        // rho_f * |D_nt| * (c_n - c_v) / dn is added implicitly and subtracted with the lagged values; this damps the deferred correction and cancels at convergence
        // D_disp is lagged (evaluated with the previous velocity)
        // the boundary half-edges carry no dispersive flux

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let den_dom = &vars.scldom[self.den_id];
        let unk_dom = &vars.scldom[self.unk_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let norm_n = &tri.face_norm[vid_v];

            // get vertex properties
            let den_v = den_dom.value[vid_v];
            let (unk_x_v, unk_y_v) = unk_dom.grad[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (nx, ny) = norm_n[loc];

                // get neighbor properties
                let den_n = den_dom.value[vid_n];
                let (unk_x_n, unk_y_n) = unk_dom.grad[vid_n];

                // get face dispersion coefficients and tangential derivative
                let (disp_nn, disp_nt) = pore_helper::disp_face(vars, self.vel_id, self.epsl_id, self.displ_id, self.dispt_id, vid_v, vid_n, (nx, ny));
                let den_f = 0.5 * (den_v + den_n);
                let dt_unk = 0.5 * (-ny * (unk_x_v + unk_x_n) + nx * (unk_y_v + unk_y_n));

                // load into matrix and vector
                let coeff = den_f * disp_nn * face / dist;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -coeff)?;
                self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -den_f * disp_nt * dt_unk * face)?;

                // stabilize the deferred cross term; implicit rho_f * |D_nt| * dl / dn, cancelled with the lagged values
                let coeff_s = den_f * disp_nt.abs() * face / dist;
                let dc_unk = unk_dom.value[vid_n] - unk_dom.value[vid_v];
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_n,  coeff_s)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -coeff_s)?;
                self.add_scl_b(vars, b_vec, self.unk_id, vid_v, coeff_s * dc_unk)?;
            }
        }

        Ok(())
    }
}
