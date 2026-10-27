use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::operator::extend::pore_helper;

#[derive(Default)]
pub struct PoreOpVecDomTranspose {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // superficial velocity (vel); unknown

    // scalars
    pub visc_id: usize,  // viscosity (mu)
    pub epsl_id: usize,  // porosity (epsl)
}

impl PoreOpVecDomTranspose {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, visc_id: usize, epsl_id: usize) -> Result<PoreOpVecDomTranspose, VFVChemError> {
        // adds the transpose diffusion term to the Darcy-Brinkman-Forchheimer equation
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        //                                                                               ^^^^^^^^^^^^^^              ^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpVecDomTranspose::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;
        oper.epsl_id = epsl_id;
        Ok(oper)
    }
}

impl OperatorBase for PoreOpVecDomTranspose {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf div( mu_eff * grad(vel)^T ) dA
        // int_lin mu_eff * grad(vel)^T . n dl
        // sum_ngh mu_eff * grad(vel)^T . n dl
        // sum_ngh mu_eff_f * grad(vel)^T_f . n * dl
        //
        // apply face skewness correction
        // mu_eff_f = (mu_eff_v + mu_eff_n) / 2 + (grad(mu_eff_v) + grad(mu_eff_n)) / 2 . skew

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
            let (u_x_v, u_y_v, v_x_v, v_y_v) = vel_dom.grad[vid_v];

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
                let (u_x_n, u_y_n, v_x_n, v_y_n) = vel_dom.grad[vid_n];

                // get face viscosity
                let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n));  // skew correction
                let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity

                // get face velocity gradient; average of both sides
                let u_x_f = 0.5 * (u_x_v + u_x_n);
                let u_y_f = 0.5 * (u_y_v + u_y_n);
                let v_x_f = 0.5 * (v_x_v + v_x_n);
                let v_y_f = 0.5 * (v_y_v + v_y_n);

                // add to sum
                sumx += visc_f * (u_x_f * nx + v_x_f * ny) * face;
                sumy += visc_f * (u_y_f * nx + v_y_f * ny) * face;
            }

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -sumx)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -sumy)?;

        }

        Ok(())
    }
}
