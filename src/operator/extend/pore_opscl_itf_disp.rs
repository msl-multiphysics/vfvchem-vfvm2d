use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclItfDispersion {
    // interface
    pub itf_id: usize,

    // scalars - rows
    pub unk_b_id: usize,  // unknown scalar (c_i^b)

    // scalars - this side
    pub unk_id: usize,    // unknown scalar (c_i); unknown
    pub den_id: usize,    // density (rho)
    pub epsl_id: usize,   // porosity (epsl)
    pub displ_id: usize,  // longitudinal dispersivity (alpha_L)
    pub dispt_id: usize,  // transverse dispersivity (alpha_T)

    // vectors - this side
    pub vel_id: usize,  // superficial velocity (vel)
}

impl PoreOpSclItfDispersion {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, unk_b_id: usize, den_id: usize, vel_id: usize, epsl_id: usize, displ_id: usize, dispt_id: usize) -> Result<PoreOpSclItfDispersion, VFVChemError> {
        // adds the mechanical dispersion term of a porous medium to the combined control volumes of interface nodes
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( -G_eff,ij * grad(c_j) - rho * D_disp * grad(c_i) ) + epsl * R_v,i + a_v * R_s,i
        //                                                                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // D_disp = alpha_T * |v| * I + (alpha_L - alpha_T) * v v^T / |v| with the interstitial velocity v = vel / epsl
        // covers the half of the control volume on the side of the velocity; loads into the rows of domain b
        // use once for each side; the interface flux leaves one half and enters the other, so it cancels

        // input checking
        // the unknown on this side is either unk_a_id or unk_b_id
        let dom_id = vars.vecdom[vel_id].dom_id;
        let unk_id = if vars.scldom[unk_a_id].dom_id == dom_id { unk_a_id } else { unk_b_id };
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;
        check_dom_scldom(vars, dom_id, displ_id)?;
        check_dom_scldom(vars, dom_id, dispt_id)?;

        // create struct
        let mut oper = PoreOpSclItfDispersion::default();
        oper.itf_id = itf_id;
        oper.unk_b_id = unk_b_id;
        oper.unk_id = unk_id;
        oper.den_id = den_id;
        oper.vel_id = vel_id;
        oper.epsl_id = epsl_id;
        oper.displ_id = displ_id;
        oper.dispt_id = dispt_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclItfDispersion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -div( -rho * D_disp * grad(c) ) dA
        // sum_ngh rho_f * (D_nn,f * (c_n - c_v) / dn + D_nt,f * grad(c)_f . t) * dl
        //
        // D_nn is implicit; the cross term D_nt is deferred and stabilized as in PoreOpSclDomDispersion
        // the boundary half-edges carry no dispersive flux

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_id = vars.scldom[self.unk_id].dom_id;
        let dom = &vars.dom[dom_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];
        let tri = &vars.tri[dom_id];
        let den_dom = &vars.scldom[self.den_id];
        let unk_dom = &vars.scldom[self.unk_id];

        // iterate through nodes
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
            let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

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
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_id, vid_v, -coeff)?;
                self.add_scl_b(vars, b_vec, self.unk_b_id, vid_b_v, -den_f * disp_nt * dt_unk * face)?;

                // stabilize the deferred cross term; implicit rho_f * |D_nt| * dl / dn, cancelled with the lagged values
                let coeff_s = den_f * disp_nt.abs() * face / dist;
                let dc_unk = unk_dom.value[vid_n] - unk_dom.value[vid_v];
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_id, vid_n,  coeff_s)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_id, vid_v, -coeff_s)?;
                self.add_scl_b(vars, b_vec, self.unk_b_id, vid_b_v, coeff_s * dc_unk)?;
            }
        }

        Ok(())
    }
}
