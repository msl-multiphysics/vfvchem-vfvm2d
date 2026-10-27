use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::operator::extend::pore_helper;

const NEG_FRAC_TWO_THREE: f64 = -2.0/3.0;

#[derive(Default)]
pub struct PoreOpVecDomDilatatory {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // superficial velocity (vel); unknown

    // scalars
    pub visc_id: usize,  // viscosity (mu)
    pub epsl_id: usize,  // porosity (epsl)
}

impl PoreOpVecDomDilatatory {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, visc_id: usize, epsl_id: usize) -> Result<PoreOpVecDomDilatatory, VFVChemError> {
        // adds the dilatatory diffusion term to the Darcy-Brinkman-Forchheimer equation
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        //                                                                               ^^^^^^^^^^^^^^                              ^^^^^^^^^^^^^^^^^^
        // lam = -2/3 * mu_eff

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpVecDomDilatatory::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;
        oper.epsl_id = epsl_id;
        Ok(oper)
    }
}

impl OperatorBase for PoreOpVecDomDilatatory {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf div( lam * div(vel) * I ) dA
        // int_lin lam * div(vel) * I . n dl
        // int_lin lam * div(vel) * n dl
        // sum_ngh lam * div(vel) * n dl
        // sum_ngh lam_f * div(vel)_f * n * dl
        //
        // apply face skewness correction
        // lam_f = (lam_v + lam_n) / 2 + (grad(lam_v) + grad(lam_n)) / 2 . skew

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let face_n = &tri.face_len[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];

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

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -sumx)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -sumy)?;

        }

        Ok(())
    }
}
