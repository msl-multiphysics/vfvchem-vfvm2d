use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpSclDomDivergence {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub pres_id: usize,  // pressure (p)
    pub den_id: usize,  // density (rho)
}

impl OpSclDomDivergence {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, pres_id: usize, den_id: usize) -> Result<OpSclDomDivergence, VFVChemError> {
        // adds the divergence term to the continuity equation
        // d(rho)/dt = -div(rho * vel)
        //             ^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, pres_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = OpSclDomDivergence::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;
        oper.den_id = den_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_surf -d(rho)/dt dA
        // -(rho_v - rho_v_prev) * dA / dt

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let den_v = den_dom.value[vid_v];
            let den_v_prev = den_dom.value_prev[vid_v];

            // load into vector
            self.add_scl_b(vars, b_vec, self.pres_id, vid_v, (den_v - den_v_prev) * area / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for OpSclDomDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf -div(rho * vel) dA
        // int_lin -rho * vel . n dl
        // sum_ngh -rho * vel . n dl
        // sum_ngh -rho_f * vel_f . n * dl
        // 
        // apply face skewness correction
        // rho_f = (rho_v + rho_n) / 2 + (grad(rho_v) + grad(rho_n)) / 2 . skew

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let face_n = &tri.face_len[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];

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

        }

        Ok(())
    }
}
