use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclDomHeatDiffusion {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (T_f)
    pub cond_id: usize,  // fluid thermal conductivity (k_f)
    pub epsl_id: usize,  // porosity (epsl)
}

impl PoreOpSclDomHeatDiffusion {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, cond_id: usize, epsl_id: usize) -> Result<PoreOpSclDomHeatDiffusion, VFVChemError> {
        // adds the fluid conduction term of a porous medium
        // epsl * (rho cp)_f * dT_f/dt + ... = -div( -epsl * k_f * grad(T_f) ) + ...
        //                                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // same as OpSclDomDiffusion with the weighted conductivity epsl * k_f; a face with zero harmonic mean has no conduction

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, cond_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclDomHeatDiffusion::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.cond_id = cond_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclDomHeatDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -div( -epsl * k_f * grad(T_f) ) dA
        // sum_ngh (epsl * k_f)_f * (T_n - T_v) / dn * dl

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // get vertex properties
            let (diff_v, (diff_x_v, diff_y_v)) = pore_helper::weighted(vars, self.cond_id, self.epsl_id, vid_v);

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (sx, sy) = skew_n[loc];

                // get neighbor properties
                let (diff_n, (diff_x_n, diff_y_n)) = pore_helper::weighted(vars, self.cond_id, self.epsl_id, vid_n);

                // get face conductivity
                let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n));  // skew correction
                let diff_f = if diff_f > 0.0 { diff_f } else { diff_m };  // guard against negative diffusion
                let diff_f = if diff_m != 0.0 { diff_f } else { 0.0 };  // guard against conduction into a region without fluid

                // load into matrix
                let coeff = diff_f * face / dist;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -coeff)?;
            }
        }

        Ok(())
    }
}
