use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::operator::extend::pore_helper;

#[derive(Default)]
pub struct PoreOpVecBndAdvection {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // superficial velocity (vel); unknown

    // scalars
    pub den_id: usize,   // density (rho)
    pub visc_id: usize,  // viscosity (mu); for the Peclet number
    pub epsl_id: usize,  // porosity (epsl)

    // flow
    pub flw_id: usize,  // flow data; for the Rhie-Chow mass flux
}

impl PoreOpVecBndAdvection {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, den_id: usize, visc_id: usize, epsl_id: usize, flw_id: usize) -> Result<PoreOpVecBndAdvection, VFVChemError> {
        // adds the advection term to the control volumes of boundary nodes
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        //                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // the boundary half-edges are upwinded from the node and contribute nothing

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;
        check_bnd_scldom(vars, bnd_id, visc_id)?;
        check_bnd_scldom(vars, bnd_id, epsl_id)?;
        check_bnd_sclflw(vars, bnd_id, flw_id)?;

        // create struct
        let mut oper = PoreOpVecBndAdvection::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.visc_id = visc_id;
        oper.epsl_id = epsl_id;
        oper.flw_id = flw_id;
        Ok(oper)
    }
}

impl OperatorBase for PoreOpVecBndAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf -(rho / epsl) * vel . grad(vel / epsl) dA
        // int_lin -(rho / epsl_v) * (vel / epsl - vel_v / epsl_v) * (vel . n) dl
        // sum_ngh -(rho / epsl_v) * (vel / epsl - vel_v / epsl_v) * (vel . n) dl + sum_bnd -(rho / epsl_v) * (vel / epsl - vel_v / epsl_v) * (vel . n) dl
        //
        // introduce mass flux
        // F = rho * vel * dl
        // sum_ngh -F_f / epsl_v * (vel_f / epsl_f - vel_v / epsl_v) + sum_bnd -F_b / epsl_v * (vel_v / epsl_v - vel_v / epsl_v)
        // sum_ngh -F_f / epsl_v * (vel_f / epsl_f - vel_v / epsl_v)
        //
        // apply power law scheme
        // A(|Pe|) = max(0, (1 - 0.1 * |Pe|)^5)
        // Pe_f = F_f / (epsl_v * epsl_f * D_f) with 1 / epsl_f = (1 / epsl_v + 1 / epsl_n) / 2
        // sum_ngh D_f * (A(|Pe_f|) - 1) * (vel_n - vel_v) + max(-F_f, 0) / epsl_v * (vel_n / epsl_n - vel_v / epsl_v)

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let epsl_dom = &vars.scldom[self.epsl_id];
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
            let (visc_v, (visc_x_v, visc_y_v)) = pore_helper::visc_eff(vars, self.visc_id, self.epsl_id, vid_v);
            let epsl_v = epsl_dom.value[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (sx, sy) = skew_n[loc];

                // get neighbor properties
                let (visc_n, (visc_x_n, visc_y_n)) = pore_helper::visc_eff(vars, self.visc_id, self.epsl_id, vid_n);
                let epsl_n = epsl_dom.value[vid_n];

                // get face viscosity and inverse porosity
                let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n));  // skew correction
                let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity
                let ieps_f = 0.5 * (1.0 / epsl_v + 1.0 / epsl_n);  // arithmetic mean of 1 / epsl

                // get power law term
                let cond = visc_f * face / dist;
                let flux = flux_n[loc];
                let peclet = if cond != 0.0 { ieps_f * flux / (epsl_v * cond) } else { 0.0 };
                let pow_law = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff_n = cond * (pow_law - 1.0) + (-flux).max(0.0) / (epsl_v * epsl_n);
                let coeff_v = cond * (pow_law - 1.0) + (-flux).max(0.0) / (epsl_v * epsl_v);
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_n,  coeff_n)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v, -coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_n,  coeff_n)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v, -coeff_v)?;
            }

        }

        Ok(())
    }
}
