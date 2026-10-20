use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecBndAdvection {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub den_id: usize,   // density (rho)
    pub visc_id: usize,  // viscosity (mu); for the Peclet number

    // flow
    pub flw_id: usize,  // flow data; for the Rhie-Chow mass flux
}

impl OpVecBndAdvection {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, den_id: usize, visc_id: usize, flw_id: usize) -> Result<OpVecBndAdvection, VFVChemError> {
        // adds the advection term to the control volumes of boundary nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                   ^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;
        check_bnd_scldom(vars, bnd_id, visc_id)?;
        check_bnd_sclflw(vars, bnd_id, flw_id)?;

        // create struct
        let mut oper = OpVecBndAdvection::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.visc_id = visc_id;
        oper.flw_id = flw_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecBndAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -div(rho * vel * vel) dA
        // int_lin -rho * vel * (vel . n) dl
        // sum_ngh -rho * vel * (vel . n) dl + sum_bnd -rho * vel * (vel . n) dl
        //
        // introduce mass flux
        // F = rho * vel * dl
        // sum_ngh -F * (vel . n) + sum_bnd -F * (vel . n)
        // sum_ngh -F_f * vel_f + sum_bnd -F_b * vel_v
        // F_b = (0.75 * rho_v * vel_v + 0.25 * rho_n * vel_n) . n * dl
        //
        // apply power law scheme
        // A(|Pe|) = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (D_f * (A(|Pe_f|) - 1) + max(-F_f, 0)) * vel_n - (D_f * (A(|Pe_f|) - 1) + max(F_f, 0)) * vel_v

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        let visc_dom = &vars.scldom[self.visc_id];
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
            let list_b = &tri.bnd_vid[vid_v];
            let half_b = &tri.bnd_half[vid_v];
            let norm_b = &tri.bnd_norm[vid_v];

            // get vertex properties
            let (u_v, v_v) = vel_dom.value[vid_v];
            let den_v = den_dom.value[vid_v];
            let visc_v = visc_dom.value[vid_v];
            let (visc_x_v, visc_y_v) = visc_dom.grad[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (sx, sy) = skew_n[loc];

                // get neighbor properties
                let visc_n = visc_dom.value[vid_n];
                let (visc_x_n, visc_y_n) = visc_dom.grad[vid_n];

                // get face viscosity
                let visc_m = 2.0 * visc_v * visc_n / (visc_v + visc_n);  // harmonic mean
                let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n));  // skew correction
                let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity

                // get power law term
                let cond = visc_f * face / dist;
                let flux = flux_n[loc];
                let peclet = if cond != 0.0 { flux / cond } else { 0.0 };
                let pow_law = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff_n = cond * (pow_law - 1.0) + (-flux).max(0.0);
                let coeff_v = cond * (pow_law - 1.0) + flux.max(0.0);
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_n,  coeff_n)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v, -coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_n,  coeff_n)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v, -coeff_v)?;
            }

            // iterate through boundary half-edges
            for (blc, &vid_n) in list_b.iter().enumerate() {
                // get geometry
                let half = half_b[blc];
                let (nx, ny) = norm_b[blc];

                // get neighbor properties
                let (u_n, v_n) = vel_dom.value[vid_n];
                let den_n = den_dom.value[vid_n];

                // get boundary mass flux
                let flux = (0.75 * den_v * (u_v * nx + v_v * ny) + 0.25 * den_n * (u_n * nx + v_n * ny)) * half;

                // load into matrix
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v, -flux)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v, -flux)?;
            }

        }

        Ok(())
    }
}
