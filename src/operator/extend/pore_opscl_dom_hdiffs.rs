use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclDomSolidDiffusion {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (T_s)
    pub cond_id: usize,  // solid thermal conductivity (k_s)
    pub epsl_id: usize,  // porosity (epsl); weight is 1 - epsl
}

impl PoreOpSclDomSolidDiffusion {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, cond_id: usize, epsl_id: usize) -> Result<PoreOpSclDomSolidDiffusion, VFVChemError> {
        // adds the solid conduction term of a porous medium
        // (1 - epsl) * (rho c)_s * dT_s/dt = -div( -(1 - epsl) * k_s * grad(T_s) ) + ...
        //                                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // same as OpSclDomDiffusion with the weighted conductivity (1 - epsl) * k_s; a face with zero harmonic mean has no conduction
        // the guard keeps the skew correction from leaking heat into a region without solid (epsl = 1)

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, cond_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclDomSolidDiffusion::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.cond_id = cond_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclDomSolidDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf -div( -(1 - epsl) * k_s * grad(T_s) ) dA
        // sum_ngh ((1 - epsl) * k_s)_f * (T_n - T_v) / dn * dl

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
            let (diff_v, (diff_x_v, diff_y_v)) = pore_helper::solid_weighted(vars, self.cond_id, self.epsl_id, vid_v);

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (sx, sy) = skew_n[loc];

                // get neighbor properties
                let (diff_n, (diff_x_n, diff_y_n)) = pore_helper::solid_weighted(vars, self.cond_id, self.epsl_id, vid_n);

                // get face conductivity
                let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n));  // skew correction
                let diff_f = if diff_f > 0.0 { diff_f } else { diff_m };  // guard against negative diffusion
                let diff_f = if diff_m != 0.0 { diff_f } else { 0.0 };  // guard against conduction into a region without solid

                // load into matrix
                let coeff = diff_f * face / dist;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -coeff)?;
            }
        }

        Ok(())
    }
}
