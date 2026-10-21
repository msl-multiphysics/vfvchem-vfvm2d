use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfDivergence {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,   // velocity (vel^a); unknown
    pub pres_a_id: usize,  // pressure (p^a)
    pub den_a_id: usize,   // density (rho^a)

    // domain b
    pub vel_b_id: usize,   // velocity (vel^b); unknown
    pub pres_b_id: usize,  // pressure (p^b)
    pub den_b_id: usize,   // density (rho^b)
}

impl OpSclItfDivergence {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, pres_a_id: usize, den_a_id: usize, vel_b_id: usize, pres_b_id: usize, den_b_id: usize) -> Result<OpSclItfDivergence, VFVChemError> {
        // adds the divergence term to the combined control volumes of interface nodes
        // d(rho)/dt = -div(rho * vel)
        //             ^^^^^^^^^^^^^^^
        // each domain contributes the faces of its half of the control volume; loads into the rows of domain b
        // half-edges on the interface are inside the combined control volume and are skipped

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, pres_a_id, pres_b_id)?;
        check_itf_scldom(vars, itf_id, den_a_id, den_b_id)?;

        // create struct
        let mut oper = OpSclItfDivergence::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.pres_a_id = pres_a_id;
        oper.den_a_id = den_a_id;
        oper.vel_b_id = vel_b_id;
        oper.pres_b_id = pres_b_id;
        oper.den_b_id = den_b_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_surf -d(rho)/dt dA
        // -(rho_v - rho_v_prev) * dA / dt

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.pres_b_id].dom_id];

        // iterate through both sides of the interface
        for den_id in [self.den_a_id, self.den_b_id] {
            // get objects on this side
            let dom_id = vars.scldom[den_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];
            let den_dom = &vars.scldom[den_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let area = tri.cvol_area[vid_v];

                // get vertex properties
                let den_v = den_dom.value[vid_v];
                let den_v_prev = den_dom.value_prev[vid_v];

                // load into vector
                self.add_scl_b(vars, b_vec, self.pres_b_id, vid_b_v, (den_v - den_v_prev) * area / dt)?;
            }
        }

        Ok(())
    }
}

impl OperatorBase for OpSclItfDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -div(rho * vel) dA
        // int_lin -rho * vel . n dl
        // sum_ngh -rho * vel . n dl + sum_bnd -rho * vel . n dl
        // sum_ngh -rho_f * vel_f . n * dl + sum_bnd -(0.75 * rho_v * vel_v + 0.25 * rho_n * vel_n) . n * dl
        //
        // apply face skewness correction
        // rho_f = (rho_v + rho_n) / 2 + (grad(rho_v) + grad(rho_n)) / 2 . skew

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.pres_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, den_id) in [(self.vel_a_id, self.den_a_id), (self.vel_b_id, self.den_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];
            let den_dom = &vars.scldom[den_id];

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
                    self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 0, vid_n, -coeff * nx)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 0, vid_v, -coeff * nx)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 1, vid_n, -coeff * ny)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 1, vid_v, -coeff * ny)?;

                    // vertex and neighbor gradient contributions; (vel_m - vel_s) terms for s = v, n
                    for vid_s in [vid_v, vid_n] {
                        for (loc_m, &vid_m) in tri.ngh_vid[vid_s].iter().enumerate() {
                            // project gradient weight onto skewness vector
                            let (w_x, w_y) = tri.grad_wgt[vid_s][loc_m];
                            let w_skew = w_x * sx + w_y * sy;

                            // get coefficients for x and y components
                            let val_x = -coeff * w_skew * nx;
                            let val_y = -coeff * w_skew * ny;

                            // load into matrix
                            self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 0, vid_m,  val_x)?;
                            self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 0, vid_s, -val_x)?;
                            self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 1, vid_m,  val_y)?;
                            self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 1, vid_s, -val_y)?;
                        }
                    }
                }

                // iterate through boundary half-edges off the interface
                for (blc, &vid_n) in list_b.iter().enumerate() {
                    if itf.nds_gnd.contains(&dom.nds_gnd[vid_n - dom.num_pts()]) { continue; }

                    // get geometry
                    let half = half_b[blc];
                    let (nx, ny) = norm_b[blc];

                    // get neighbor properties
                    let den_n = den_dom.value[vid_n];

                    // load into matrix
                    let coeff_v = 0.75 * den_v * half;
                    let coeff_n = 0.25 * den_n * half;
                    self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 0, vid_v, -coeff_v * nx)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 0, vid_n, -coeff_n * nx)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 1, vid_v, -coeff_v * ny)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, 1, vid_n, -coeff_n * ny)?;
                }

            }
        }

        Ok(())
    }
}
