use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclItfSolidDiffusion {
    // interface
    pub itf_id: usize,

    // scalars - rows
    pub row_id: usize,  // row scalar (T_s^b with continuity; T_s on this side without)

    // scalars - this side
    pub unk_id: usize,   // unknown scalar (T_s)
    pub cond_id: usize,  // solid thermal conductivity (k_s)
    pub epsl_id: usize,  // porosity (epsl); weight is 1 - epsl
}

impl PoreOpSclItfSolidDiffusion {
    pub fn new(vars: &Variables, itf_id: usize, unk_id: usize, row_id: usize, cond_id: usize, epsl_id: usize) -> Result<PoreOpSclItfSolidDiffusion, VFVChemError> {
        // adds the solid conduction term of a porous medium to the interface nodes
        // (1 - epsl) * (rho c)_s * dT_s/dt = -div( -(1 - epsl) * k_s * grad(T_s) ) + ...
        //                                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // same as PoreOpSclDomSolidDiffusion; covers the faces of the half of the control volume on the side of the unknown; loads into the rows of row_id
        // row_id is T_s^b for a solid temperature with continuity (use once for each side), or unk_id for an adiabatic solid surface
        // the interface flux leaves one half and enters the other with continuity, so it cancels; without continuity, it is zero

        // input checking
        let dom_id = vars.scldom[unk_id].dom_id;
        let itf = &vars.itf[itf_id];
        if row_id != unk_id { check_itf_scldom(vars, itf_id, unk_id, row_id)?; }
        else if dom_id != itf.dom_a_id && dom_id != itf.dom_b_id { return Err(VFVChemError::OperatorScalarDomainNotOnInterface { scldom_a_id: unk_id, scldom_b_id: row_id, itf_id }) }
        check_dom_scldom(vars, dom_id, cond_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclItfSolidDiffusion::default();
        oper.itf_id = itf_id;
        oper.row_id = row_id;
        oper.unk_id = unk_id;
        oper.cond_id = cond_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclItfSolidDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf -div( -(1 - epsl) * k_s * grad(T_s) ) dA
        // sum_ngh ((1 - epsl) * k_s)_f * (T_n - T_v) / dn * dl

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_id = vars.scldom[self.unk_id].dom_id;
        let dom = &vars.dom[dom_id];
        let dom_r = &vars.dom[vars.scldom[self.row_id].dom_id];
        let tri = &vars.tri[dom_id];

        // iterate through nodes
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
            let vid_r_v = dom_r.num_pts() + dom_r.nds_nid[&gnd_v];

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
                self.add_sclscl_a(vars, a_triplet, self.row_id, vid_r_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.row_id, vid_r_v, self.unk_id, vid_v, -coeff)?;
            }
        }

        Ok(())
    }
}
