use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::operator::extend::pore_helper;

#[derive(Default)]
pub struct PoreOpVecDomAdvection {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // superficial velocity (vel); unknown

    // scalars
    pub den_id: usize,   // density (rho)
    pub visc_id: usize,  // viscosity (mu); for the Peclet number
    pub epsl_id: usize,  // porosity (epsl)

    // flow
    pub flw_id: usize,  // flow data; for the Rhie-Chow mass flux
}

impl PoreOpVecDomAdvection {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, den_id: usize, visc_id: usize, epsl_id: usize, flw_id: usize) -> Result<PoreOpVecDomAdvection, VFVChemError> {
        // adds the advection term to the Darcy-Brinkman-Forchheimer equation
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        //                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;
        check_dom_sclflw(vars, dom_id, flw_id)?;

        // create struct
        let mut oper = PoreOpVecDomAdvection::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.visc_id = visc_id;
        oper.epsl_id = epsl_id;
        oper.flw_id = flw_id;
        Ok(oper)
    }
}

impl OperatorBase for PoreOpVecDomAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf -(rho / epsl) * vel . grad(vel / epsl) dA
        // int_lin -(rho / epsl_v) * (vel / epsl - vel_v / epsl_v) * (vel . n) dl
        // sum_ngh -(rho / epsl_v) * (vel / epsl - vel_v / epsl_v) * (vel . n) dl
        //
        // introduce mass flux
        // F = rho * vel * dl
        // sum_ngh -F_f / epsl_v * (vel_f / epsl_f - vel_v / epsl_v)
        //
        // apply power law scheme
        // A(|Pe|) = max(0, (1 - 0.1 * |Pe|)^5)
        // Pe_f = F_f / (epsl_v * epsl_f * D_f) with 1 / epsl_f = (1 / epsl_v + 1 / epsl_n) / 2
        // sum_ngh D_f * (A(|Pe_f|) - 1) * (vel_n - vel_v) + max(-F_f, 0) / epsl_v * (vel_n / epsl_n - vel_v / epsl_v)

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let epsl_dom = &vars.scldom[self.epsl_id];
        let flow_dom = &vars.sclflw[self.flw_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
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
