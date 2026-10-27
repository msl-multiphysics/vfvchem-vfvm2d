use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::operator::extend::pore_helper;

const NEG_FRAC_TWO_THREE: f64 = -2.0/3.0;

#[derive(Default)]
pub struct PoreOpVecBndDilatatory {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // superficial velocity (vel); unknown

    // scalars
    pub visc_id: usize,  // viscosity (mu)
    pub epsl_id: usize,  // porosity (epsl)
}

impl PoreOpVecBndDilatatory {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, visc_id: usize, epsl_id: usize) -> Result<PoreOpVecBndDilatatory, VFVChemError> {
        // adds the dilatatory diffusion term to the control volumes of boundary nodes
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        //                                                                               ^^^^^^^^^^^^^^                              ^^^^^^^^^^^^^^^^^^
        // lam = -2/3 * mu_eff

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, visc_id)?;
        check_bnd_scldom(vars, bnd_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpVecBndDilatatory::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;
        oper.epsl_id = epsl_id;
        Ok(oper)
    }
}

impl OperatorBase for PoreOpVecBndDilatatory {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf div( lam * div(vel) * I ) dA
        // int_lin lam * div(vel) * I . n dl
        // int_lin lam * div(vel) * n dl
        // sum_ngh lam * div(vel) * n dl + sum_bnd lam * div(vel) * n dl
        // sum_ngh lam_f * div(vel)_f * n * dl + sum_bnd lam_b * div(vel)_b * n * dl
        //
        // apply face skewness correction
        // lam_f = (lam_v + lam_n) / 2 + (grad(lam_v) + grad(lam_n)) / 2 . skew
        //
        // apply half-edge interpolation
        // lam_b = 0.75 * lam_v + 0.25 * lam_n

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let face_n = &tri.face_len[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];
            let list_b = &tri.bnd_vid[vid_v];
            let half_b = &tri.bnd_half[vid_v];
            let norm_b = &tri.bnd_norm[vid_v];

            // get vertex properties
            let (visc_v, (visc_x_v, visc_y_v)) = pore_helper::visc_eff(vars, self.visc_id, self.epsl_id, vid_v);
            let (u_x_v, _, _, v_y_v) = vel_dom.grad[vid_v];
            let div_v = u_x_v + v_y_v;

            // iterate through faces
            let mut sumx = 0.0;
            let mut sumy = 0.0;
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let (nx, ny) = norm_n[loc];
                let (sx, sy) = skew_n[loc];

                // get neighbor properties
                let (visc_n, (visc_x_n, visc_y_n)) = pore_helper::visc_eff(vars, self.visc_id, self.epsl_id, vid_n);
                let (u_x_n, _, _, v_y_n) = vel_dom.grad[vid_n];
                let div_n = u_x_n + v_y_n;

                // get face second viscosity
                let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n));  // skew correction
                let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity
                let lam_f = NEG_FRAC_TWO_THREE * visc_f;

                // get face divergence
                let div_f = 0.5 * (div_v + div_n);

                // add to sum
                let coeff = lam_f * div_f * face;
                sumx += coeff * nx;
                sumy += coeff * ny;
            }

            // iterate through boundary half-edges
            for (blc, &vid_n) in list_b.iter().enumerate() {
                // get geometry
                let half = half_b[blc];
                let (nx, ny) = norm_b[blc];

                // get neighbor properties
                let (visc_n, _) = pore_helper::visc_eff(vars, self.visc_id, self.epsl_id, vid_n);
                let (u_x_n, _, _, v_y_n) = vel_dom.grad[vid_n];
                let div_n = u_x_n + v_y_n;

                // get half-edge second viscosity
                let visc_b = 0.75 * visc_v + 0.25 * visc_n;
                let lam_b = NEG_FRAC_TWO_THREE * visc_b;

                // get half-edge divergence
                let div_b = 0.75 * div_v + 0.25 * div_n;

                // add to sum
                let coeff = lam_b * div_b * half;
                sumx += coeff * nx;
                sumy += coeff * ny;
            }

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -sumx)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -sumy)?;

        }

        Ok(())
    }
}
