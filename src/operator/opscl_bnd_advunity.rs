use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndAdvectionUnity {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,           // unknown scalar (c_i); unknown
    pub diff_id: Option<usize>,  // self-diffusivity (D_ii); for the Peclet number
    pub den_id: usize,           // density (rho); converts the mass flux

    // vectors
    pub vel_id: usize,  // velocity (vel)

    // flow
    pub flw_id: usize,  // flow data; for the Rhie-Chow mass flux
}

impl OpSclBndAdvectionUnity {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, diff_id: Option<usize>, vel_id: usize, den_id: usize, flw_id: usize) -> Result<OpSclBndAdvectionUnity, VFVChemError> {
        // adds the advection term with unit weight to the control volumes of boundary nodes
        // dc_i/dt + vel . grad(c_i) = -div( -D_ij * grad(c_j) ) + R_i
        //           ^^^^^^^^^^^^^^^
        // must be used together with OpSclBndDiffusion; the power law scheme corrects its central difference
        // a missing self-diffusivity sets D_f to zero and falls back to pure upwinding
        // the boundary half-edges are upwinded from the node and contribute nothing

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        if let Some(diff_id) = diff_id {
            check_bnd_scldom(vars, bnd_id, diff_id)?;
        }
        check_bnd_scldom(vars, bnd_id, den_id)?;
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_sclflw(vars, bnd_id, flw_id)?;

        // create struct
        let mut oper = OpSclBndAdvectionUnity::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.diff_id = diff_id;
        oper.den_id = den_id;
        oper.vel_id = vel_id;
        oper.flw_id = flw_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndAdvectionUnity {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -vel . grad(c) dA
        // int_lin -(c - c_v) * (vel . n) dl
        // sum_ngh -(c - c_v) * (vel . n) dl + sum_bnd -(c - c_v) * (vel . n) dl
        //
        // introduce volume flux from the mass flux
        // F = vel * dl = (1 / rho) * rho * vel * dl
        // sum_ngh -F_f * (c_f - c_v) + sum_bnd -F_b * (c_v - c_v)
        // sum_ngh -F_f * (c_f - c_v)
        //
        // apply power law scheme
        // A(|Pe|) = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (D_f * (A(|Pe_f|) - 1) + max(-F_f, 0)) * (c_n - c_v)
        //
        // non-conservative form; uniform c is preserved even if continuity is not satisfied

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let diff_dom = self.diff_id.map(|diff_id| &vars.scldom[diff_id]);
        let den_dom = &vars.scldom[self.den_id];
        let flow_dom = &vars.sclflw[self.flw_id];

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
            let flux_n = &flow_dom.flux[vid_v];

            // get vertex properties
            let den_v = den_dom.value[vid_v];
            let (den_x_v, den_y_v) = den_dom.grad[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (sx, sy) = skew_n[loc];

                // get neighbor properties
                let den_n = den_dom.value[vid_n];
                let (den_x_n, den_y_n) = den_dom.grad[vid_n];

                // get face density
                let den_f = 0.5 * (den_v + den_n) + 0.5 * (sx * (den_x_v + den_x_n) + sy * (den_y_v + den_y_n));  // skew-corrected arithmetic mean

                // get face diffusion coefficient; zero if absent
                let diff_f = match diff_dom {
                    Some(diff_dom) => {
                        let diff_v = diff_dom.value[vid_v];
                        let (diff_x_v, diff_y_v) = diff_dom.grad[vid_v];
                        let diff_n = diff_dom.value[vid_n];
                        let (diff_x_n, diff_y_n) = diff_dom.grad[vid_n];
                        let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                        let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n));  // skew correction
                        if diff_f > 0.0 { diff_f } else { diff_m }  // guard against negative diffusion
                    }
                    None => 0.0,
                };

                // get power law term
                let cond = diff_f * face / dist;
                let flux = if den_f != 0.0 { flux_n[loc] / den_f } else { 0.0 };
                let peclet = if cond != 0.0 { flux / cond } else { 0.0 };
                let pow_law = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff = cond * (pow_law - 1.0) + (-flux).max(0.0);
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -coeff)?;
            }

        }

        Ok(())
    }
}
