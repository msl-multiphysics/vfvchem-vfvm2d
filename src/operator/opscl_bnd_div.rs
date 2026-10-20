use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndDivergence {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub pres_id: usize,  // pressure (p)
    pub den_id: usize,  // density (rho)
}

impl OpSclBndDivergence {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, pres_id: usize, den_id: usize) -> Result<OpSclBndDivergence, VFVChemError> {
        // adds the divergence term to the control volumes of boundary nodes
        // d(rho)/dt = -div(rho * vel)
        //             ^^^^^^^^^^^^^^^

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, pres_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;

        // create struct
        let mut oper = OpSclBndDivergence::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;
        oper.den_id = den_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -div(rho * vel) dA
        // int_lin -rho * vel . n dl
        // sum_ngh -rho * vel . n dl + sum_bnd -rho * vel . n dl
        // sum_ngh -rho_f * vel_f . n * dl + sum_bnd -(0.75 * rho_v * vel_v + 0.25 * rho_n * vel_n) . n * dl
        //
        // apply face skewness correction
        // rho_f = (rho_v + rho_n) / 2 + (grad(rho_v) + grad(rho_n)) / 2 . skew

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let den_dom = &vars.scldom[self.den_id];

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
            let den_v = den_dom.value[vid_v];
            let (den_x_v, den_y_v) = den_dom.grad[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let (nx, ny) = norm_n[loc];
                let (sx, sy) = skew_n[loc];

                // get neighbor properties
                let den_n = den_dom.value[vid_n];
                let (den_x_n, den_y_n) = den_dom.grad[vid_n];

                // get face density
                let den_m = 0.5 * (den_v + den_n);  // arithmetic mean
                let den_f = den_m + 0.5 * (sx * (den_x_v + den_x_n) + sy * (den_y_v + den_y_n));  // skewness correction

                // load into matrix
                let coeff = 0.5 * den_f * face;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_n, -coeff * nx)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_v, -coeff * nx)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_n, -coeff * ny)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_v, -coeff * ny)?;

                // vertex gradient contribution
                let list_vm = &tri.ngh_vid[vid_v];
                let wgt_vm = &tri.grad_wgt[vid_v];
                for (loc_m, &vid_m) in list_vm.iter().enumerate() {
                    // project gradient weight onto skewness vector
                    let (w_x, w_y) = wgt_vm[loc_m];
                    let w_skew = w_x * sx + w_y * sy;

                    // get coefficients for x and y components
                    let val_x = -coeff * w_skew * nx;
                    let val_y = -coeff * w_skew * ny;

                    // (vel_m - vel_v) terms
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_m,  val_x)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_v, -val_x)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_m,  val_y)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_v, -val_y)?;
                }

                // neighbor gradient contribution
                let list_nm = &tri.ngh_vid[vid_n];
                let wgt_nm = &tri.grad_wgt[vid_n];
                for (loc_m, &vid_m) in list_nm.iter().enumerate() {
                    // project gradient weight onto skewness vector
                    let (w_x, w_y) = wgt_nm[loc_m];
                    let w_skew = w_x * sx + w_y * sy;

                    // get coefficients for x and y components
                    let val_x = -coeff * w_skew * nx;
                    let val_y = -coeff * w_skew * ny;

                    // (vel_m - vel_n) terms
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_m,  val_x)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_n, -val_x)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_m,  val_y)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_n, -val_y)?;
                }
            }

            // iterate through boundary half-edges
            for (blc, &vid_n) in list_b.iter().enumerate() {
                // get geometry
                let half = half_b[blc];
                let (nx, ny) = norm_b[blc];

                // get neighbor properties
                let den_n = den_dom.value[vid_n];

                // load into matrix
                let coeff_v = 0.75 * den_v * half;
                let coeff_n = 0.25 * den_n * half;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_v, -coeff_v * nx)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_n, -coeff_n * nx)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_v, -coeff_v * ny)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_n, -coeff_n * ny)?;
            }

        }

        Ok(())
    }
}
