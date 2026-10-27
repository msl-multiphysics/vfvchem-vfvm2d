use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclBndDiffusion {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (c_i)
    pub diff_id: usize,  // diffusion coefficient (G_ij)
    pub drv_id: usize,   // driving scalar (c_j)
    pub epsl_id: usize,  // porosity (epsl)
    pub tort_id: usize,  // tortuosity (tau)
}

impl PoreOpSclBndDiffusion {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, diff_id: usize, drv_id: usize, epsl_id: usize, tort_id: usize) -> Result<PoreOpSclBndDiffusion, VFVChemError> {
        // adds the effective diffusion term of a porous medium to the control volumes of boundary nodes
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( -G_eff,ij * grad(c_j) - rho * D_disp * grad(c_i) ) + epsl * R_v,i + a_v * R_s,i
        //                                                ^^^^^^^^^^^^^^^^^^^^^^^^^^
        // G_eff,ij = (epsl / tau) * G_ij; G_ij includes the density (rho * D_ij) on domains with flow

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, diff_id)?;
        check_bnd_scldom(vars, bnd_id, drv_id)?;
        check_bnd_scldom(vars, bnd_id, epsl_id)?;
        check_bnd_scldom(vars, bnd_id, tort_id)?;

        // create struct
        let mut oper = PoreOpSclBndDiffusion::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.diff_id = diff_id;
        oper.drv_id = drv_id;
        oper.epsl_id = epsl_id;
        oper.tort_id = tort_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclBndDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf -div( -G_eff * grad(c) ) dA
        // sum_ngh G_eff,f * (c_n - c_v) / dn * dl
        // G_eff,v = (epsl_v / tau_v) * G_v; harmonic mean and skew correction as in OpSclDomDiffusion

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // get vertex properties
            let (diff_v, (diff_x_v, diff_y_v)) = pore_helper::diff_eff(vars, self.diff_id, self.epsl_id, self.tort_id, vid_v);

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (sx, sy) = skew_n[loc];

                // get neighbor properties
                let (diff_n, (diff_x_n, diff_y_n)) = pore_helper::diff_eff(vars, self.diff_id, self.epsl_id, self.tort_id, vid_n);

                // get face diffusion coefficient
                let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n));  // skew correction
                let diff_f = if diff_f > 0.0 { diff_f } else { diff_m };  // guard against negative diffusion

                // load into matrix
                let coeff = diff_f * face / dist;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.drv_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.drv_id, vid_v, -coeff)?;
            }
        }

        Ok(())
    }
}
