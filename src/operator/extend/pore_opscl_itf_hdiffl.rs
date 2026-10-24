use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclItfHeatDiffusion {
    // interface
    pub itf_id: usize,

    // scalars - domain a
    pub unk_a_id: usize,   // unknown scalar (T_f^a)
    pub cond_a_id: usize,  // fluid thermal conductivity (k_f^a)
    pub epsl_a_id: usize,  // porosity (epsl^a)

    // scalars - domain b
    pub unk_b_id: usize,   // unknown scalar (T_f^b)
    pub cond_b_id: usize,  // fluid thermal conductivity (k_f^b)
    pub epsl_b_id: usize,  // porosity (epsl^b)
}

impl PoreOpSclItfHeatDiffusion {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, cond_a_id: usize, epsl_a_id: usize, unk_b_id: usize, cond_b_id: usize, epsl_b_id: usize) -> Result<PoreOpSclItfHeatDiffusion, VFVChemError> {
        // adds the fluid conduction term of a porous medium to the combined control volumes of interface nodes
        // epsl * (rho cp)_f * dT_f/dt + ... = -div( -epsl * k_f * grad(T_f) ) + ...
        //                                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // same as PoreOpSclDomHeatDiffusion; each domain contributes the faces of its half of the control volume; loads into the rows of domain b
        // the interface flux leaves one half and enters the other, so it cancels

        // input checking
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_itf_scldom(vars, itf_id, cond_a_id, cond_b_id)?;
        check_itf_scldom(vars, itf_id, epsl_a_id, epsl_b_id)?;

        // create struct
        let mut oper = PoreOpSclItfHeatDiffusion::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.cond_a_id = cond_a_id;
        oper.epsl_a_id = epsl_a_id;
        oper.unk_b_id = unk_b_id;
        oper.cond_b_id = cond_b_id;
        oper.epsl_b_id = epsl_b_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclItfHeatDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -div( -epsl * k_f * grad(T_f) ) dA
        // sum_ngh (epsl * k_f)_f * (T_n - T_v) / dn * dl

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];

        // iterate through both sides of the interface
        for (unk_id, cond_id, epsl_id) in [(self.unk_a_id, self.cond_a_id, self.epsl_a_id), (self.unk_b_id, self.cond_b_id, self.epsl_b_id)] {
            // get objects on this side
            let dom_id = vars.scldom[unk_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];

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
                let skew_n = &tri.face_skew[vid_v];

                // get vertex properties
                let (diff_v, (diff_x_v, diff_y_v)) = pore_helper::weighted(vars, cond_id, epsl_id, vid_v);

                // iterate through faces
                for (loc, &vid_n) in list_n.iter().enumerate() {
                    // get geometry
                    let face = face_n[loc];
                    let dist = dist_n[loc];
                    let (sx, sy) = skew_n[loc];

                    // get neighbor properties
                    let (diff_n, (diff_x_n, diff_y_n)) = pore_helper::weighted(vars, cond_id, epsl_id, vid_n);

                    // get face conductivity
                    let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                    let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n));  // skew correction
                    let diff_f = if diff_f > 0.0 { diff_f } else { diff_m };  // guard against negative diffusion
                    let diff_f = if diff_m != 0.0 { diff_f } else { 0.0 };  // guard against conduction into a region without fluid

                    // load into matrix
                    let coeff = diff_f * face / dist;
                    self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, unk_id, vid_n,  coeff)?;
                    self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, unk_id, vid_v, -coeff)?;
                }

            }
        }

        Ok(())
    }
}
