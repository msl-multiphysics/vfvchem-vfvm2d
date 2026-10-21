use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

const NEG_FRAC_TWO_THREE: f64 = -2.0/3.0;

#[derive(Default)]
pub struct OpVecItfDilatatory {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,   // velocity (vel^a); unknown
    pub visc_a_id: usize,  // viscosity (mu^a)

    // domain b
    pub vel_b_id: usize,   // velocity (vel^b); unknown
    pub visc_b_id: usize,  // viscosity (mu^b)
}

impl OpVecItfDilatatory {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, visc_a_id: usize, vel_b_id: usize, visc_b_id: usize) -> Result<OpVecItfDilatatory, VFVChemError> {
        // adds the dilatatory diffusion term to the combined control volumes of interface nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                     ^^^^^^^^^^                              ^^^^^^^^^^^^^^^^^^
        // each domain contributes the faces of its half of the control volume; loads into the rows of domain b
        // half-edges on the interface are inside the combined control volume and are skipped

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, visc_a_id, visc_b_id)?;

        // create struct
        let mut oper = OpVecItfDilatatory::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.visc_a_id = visc_a_id;
        oper.vel_b_id = vel_b_id;
        oper.visc_b_id = visc_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecItfDilatatory {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
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
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, visc_id) in [(self.vel_a_id, self.visc_a_id), (self.vel_b_id, self.visc_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];
            let vel_dom = &vars.vecdom[vel_id];
            let visc_dom = &vars.scldom[visc_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

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
                    let visc_n = visc_dom.value[vid_n];
                    let (visc_x_n, visc_y_n) = visc_dom.grad[vid_n];
                    let (u_x_n, _, _, v_y_n) = vel_dom.grad[vid_n];
                    let div_n = u_x_n + v_y_n;

                    // get face second viscosity
                    let visc_m = 2.0 * visc_v * visc_n / (visc_v + visc_n);  // harmonic mean
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

                // iterate through boundary half-edges off the interface
                for (blc, &vid_n) in list_b.iter().enumerate() {
                    if itf.nds_gnd.contains(&dom.nds_gnd[vid_n - dom.num_pts()]) { continue; }

                    // get geometry
                    let half = half_b[blc];
                    let (nx, ny) = norm_b[blc];

                    // get neighbor properties
                    let visc_n = visc_dom.value[vid_n];
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
                self.add_vec_b(vars, b_vec, self.vel_b_id, 0, vid_b_v, -sumx)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 1, vid_b_v, -sumy)?;

            }
        }

        Ok(())
    }
}
