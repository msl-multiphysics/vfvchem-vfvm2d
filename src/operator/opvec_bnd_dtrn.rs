use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecBndTranspose {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub visc_id: usize,  // viscosity (mu)
}

impl OpVecBndTranspose {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, visc_id: usize) -> Result<OpVecBndTranspose, VFVChemError> {
        // adds the transpose diffusion term to the control volumes of boundary nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                     ^^^^^^^^^^              ^^^^^^^^^^^

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, visc_id)?;

        // create struct
        let mut oper = OpVecBndTranspose::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecBndTranspose {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf div( mu * grad(vel)^T ) dA
        // int_lin mu * grad(vel)^T . n dl
        // sum_ngh mu * grad(vel)^T . n dl + sum_bnd mu * grad(vel)^T . n dl
        // sum_ngh mu_f * grad(vel)^T_f . n * dl + sum_bnd mu_b * grad(vel)^T_b . n * dl
        //
        // apply face skewness correction
        // mu_f = (mu_v + mu_n) / 2 + (grad(mu_v) + grad(mu_n)) / 2 . skew
        //
        // apply half-edge interpolation
        // mu_b = 0.75 * mu_v + 0.25 * mu_n

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let visc_dom = &vars.scldom[self.visc_id];

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
            let visc_v = visc_dom.value[vid_v];
            let (visc_x_v, visc_y_v) = visc_dom.grad[vid_v];
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
                let visc_n = visc_dom.value[vid_n];
                let (visc_x_n, visc_y_n) = visc_dom.grad[vid_n];
                let (u_x_n, u_y_n, v_x_n, v_y_n) = vel_dom.grad[vid_n];

                // get face viscosity
                let visc_m = 2.0 * visc_v * visc_n / (visc_v + visc_n);  // harmonic mean
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

            // iterate through boundary half-edges
            for (blc, &vid_n) in list_b.iter().enumerate() {
                // get geometry
                let half = half_b[blc];
                let (nx, ny) = norm_b[blc];

                // get neighbor properties
                let visc_n = visc_dom.value[vid_n];
                let (u_x_n, u_y_n, v_x_n, v_y_n) = vel_dom.grad[vid_n];

                // get half-edge viscosity
                let visc_b = 0.75 * visc_v + 0.25 * visc_n;

                // get half-edge velocity gradient
                let u_x_b = 0.75 * u_x_v + 0.25 * u_x_n;
                let u_y_b = 0.75 * u_y_v + 0.25 * u_y_n;
                let v_x_b = 0.75 * v_x_v + 0.25 * v_x_n;
                let v_y_b = 0.75 * v_y_v + 0.25 * v_y_n;

                // add to sum
                sumx += visc_b * (u_x_b * nx + v_x_b * ny) * half;
                sumy += visc_b * (u_y_b * nx + v_y_b * ny) * half;
            }

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -sumx)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -sumy)?;

        }

        Ok(())
    }
}
