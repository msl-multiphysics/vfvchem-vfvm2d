use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::operator::extend::pore_helper;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpVecItfDiffusion {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,   // superficial velocity (vel^a); unknown
    pub visc_a_id: usize,  // viscosity (mu^a)
    pub epsl_a_id: usize,  // porosity (epsl^a)

    // domain b
    pub vel_b_id: usize,   // superficial velocity (vel^b); unknown
    pub visc_b_id: usize,  // viscosity (mu^b)
    pub epsl_b_id: usize,  // porosity (epsl^b)
}

impl PoreOpVecItfDiffusion {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, visc_a_id: usize, epsl_a_id: usize, vel_b_id: usize, visc_b_id: usize, epsl_b_id: usize) -> Result<PoreOpVecItfDiffusion, VFVChemError> {
        // adds the diffusion term to the combined control volumes of interface nodes
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        //                                                                               ^^^^^^^^^^^^^^^^^^^^^^^^^
        // each domain contributes the faces of its half of the control volume; loads into the rows of domain b
        // the interface flux leaves one half and enters the other, so it cancels

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, visc_a_id, visc_b_id)?;
        check_itf_scldom(vars, itf_id, epsl_a_id, epsl_b_id)?;

        // create struct
        let mut oper = PoreOpVecItfDiffusion::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.visc_a_id = visc_a_id;
        oper.epsl_a_id = epsl_a_id;
        oper.vel_b_id = vel_b_id;
        oper.visc_b_id = visc_b_id;
        oper.epsl_b_id = epsl_b_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpVecItfDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf div( mu_eff * grad(vel) ) dA
        // int_lin mu_eff * grad(vel) . n dl
        // sum_ngh mu_eff * grad(vel) . n dl
        //
        // apply normal derivative
        // sum_ngh mu_eff_f * (vel_n - vel_v) / dn * dl

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, visc_id, epsl_id) in [(self.vel_a_id, self.visc_a_id, self.epsl_a_id), (self.vel_b_id, self.visc_b_id, self.epsl_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let list_n = &tri.ngh_vid[vid_v];
                let dist_n = &tri.ngh_dist[vid_v];
                let face_n = &tri.face_len[vid_v];
                let skew_n = &tri.face_skew[vid_v];

                // get vertex properties
                let (visc_v, (visc_x_v, visc_y_v)) = pore_helper::visc_eff(vars, visc_id, epsl_id, vid_v);

                // iterate through faces
                for (loc, &vid_n) in list_n.iter().enumerate() {
                    // get geometry
                    let face = face_n[loc];
                    let dist = dist_n[loc];
                    let (sx, sy) = skew_n[loc];

                    // get neighbor properties
                    let (visc_n, (visc_x_n, visc_y_n)) = pore_helper::visc_eff(vars, visc_id, epsl_id, vid_n);

                    // get face viscosity
                    let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                    let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n));  // skew correction
                    let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity

                    // load into matrix
                    let coeff = visc_f * face / dist;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 0, vid_b_v, vel_id, 0, vid_n,  coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 0, vid_b_v, vel_id, 0, vid_v, -coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 1, vid_b_v, vel_id, 1, vid_n,  coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 1, vid_b_v, vel_id, 1, vid_v, -coeff)?;
                }

            }
        }

        Ok(())
    }
}
