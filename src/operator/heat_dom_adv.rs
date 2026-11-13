use crate::base::prelude::*;
use crate::operator::prelude::*;

const VENKATAKRISHNAN_SCALE: f64 = 1.0;

#[derive(Default)]
pub struct HeatDomAdvection {
    // domain
    pub dom_id: usize,

    // unknowns
    pub temp_id: usize,   // temperature (T)

    // others
    pub vel_id: usize,   // velocity (vel)
    pub den_id: usize,   // density (rho)
    pub hcap_id: usize,  // heat capacity (cp)
}

impl HeatDomAdvection {
    pub fn new(vars: &Variables, dom_id: usize, temp_id: usize, vel_id: usize, den_id: usize, hcap_id: usize) -> Result<HeatDomAdvection, VFVChemError> {
        // adds the advection term to the heat transfer equation
        // rho * cp * d(T)/dt + rho * cp * vel . grad(T) = -div( -k * grad(T) ) + L * T + R
        //                      ^^^^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, temp_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, hcap_id)?;

        // create struct
        let mut oper = HeatDomAdvection::default();
        oper.dom_id = dom_id;
        oper.temp_id = temp_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.hcap_id = hcap_id;

        Ok(oper)
    }
}

impl OperatorBase for HeatDomAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the RHS of the PDE
        // -rho * cp * vel . grad(T)
        //
        // rewrite advection term
        // -rho * cp * vel . grad(T)
        // -div(rho * cp * vel * T) + T * div(rho * cp * vel)
        //
        // get surface integral of 1st term
        // int_surf -div(rho * cp * vel * T) dA
        // int_lin -(rho * cp * vel * T) . n dl
        // sum_ngh -(rho * cp * vel * T)_f . n dl
        //
        // define enthalpy flux
        // let F = rho * cp * vel . n * dl
        // sum_ngh (-F_f * T_f) -> 1st term
        //
        // get surface integral of 2nd term
        // int_surf T * div(rho * cp * vel) dA
        // T_c * div(rho * cp * vel)_c dA
        //
        // approximate divergence using divergence theorem
        // div(rho * cp * vel)_c = (1/dA) sum_ngh (rho * cp * vel)_f . n dl
        // T_c * sum_ngh (rho * cp * vel)_f . n dl
        // T_c * sum_ngh F_f
        // 
        // apply off-center correction
        // T_c = T_v + T_o 
        // (T_v + T_o) * sum_ngh F_f -> 2nd term
        //
        // add the terms
        // sum_ngh -F_f * T_f + (T_v + T_o) * F_f
        // sum_ngh -F_f * (T_f - T_v) + F_f * T_o
        //
        // face temperature T_f is computed using venkatakrishnan limiter
        // let T_min be the minimum T_v or T_n among neighbors
        // let T_max be the maximum T_v or T_n among neighbors
        // del_f = grad(T)_v . (x_f - x_v)
        // del_lim = (T_max - T_v) if del_f > 0 else (T_min - T_v) if del_f < 0
        // eps^2 = (K * h / L)^3 * del_T^2; h = sqrt(dA), L = sqrt(sum dA), del_T = range of T in domain
        // alp_f = (del_lim^2 + 2 * del_lim * del_f + eps^2) / (del_lim^2 + del_lim * del_f + 2 * del_f^2 + eps^2)
        // alp_v = min(1, alp_f among faces)
        //
        // let T_u be the upwind vertex temperature (either T_v or T_n)
        // T_f = T_u + alp_u * grad(T)_u . (x_f - x_u)

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let temp_dom = &vars.scldom[self.temp_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        let hcap_dom = &vars.scldom[self.hcap_id];

        // step 1: limiter computation

        // initialize
        let mut alp_vec = vec![1.0; dom.num_vert()];

        // get limiter scales
        let (mut temp_min, mut temp_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for &temp in temp_dom.value.iter() {
            temp_min = temp_min.min(temp);
            temp_max = temp_max.max(temp);
        }
        let temp2_scale = (temp_max - temp_min) * (temp_max - temp_min);
        let len_scale = tri.vert_area.iter().sum::<f64>().sqrt();

        // compute limiter at each vertex
        for vid_v in 0..dom.num_vert() {
            // get bounds from neighbors
            let temp_v = temp_dom.value[vid_v];
            let (gradx_v, grady_v) = temp_dom.grad[vid_v];
            let (mut temp_min, mut temp_max) = (temp_v, temp_v);
            for &vid_n in tri.ngh_vid[vid_v].iter() {
                let temp_n = temp_dom.value[vid_n];
                temp_min = temp_min.min(temp_n);
                temp_max = temp_max.max(temp_n);
            }

            // get limiter tolerance
            let eps2 = (VENKATAKRISHNAN_SCALE * tri.vert_area[vid_v].sqrt() / len_scale).powi(3) * temp2_scale;

            // iterate through faces
            let mut alp: f64 = 1.0;
            for loc in 0..tri.ngh_vid[vid_v].len() {
                // get position of face center relative to vertex
                let dist = tri.ngh_dist[vid_v][loc];
                let (nx, ny) = tri.face_norm[vid_v][loc];
                let (sx, sy) = tri.face_skew[vid_v][loc];
                let (rx, ry) = (0.5 * dist * nx + sx, 0.5 * dist * ny + sy);  // x_f - x_v

                // get limiter
                let del = gradx_v * rx + grady_v * ry;
                if del != 0.0 {
                    let del_lim = if del > 0.0 {temp_max - temp_v} else {temp_min - temp_v};
                    let alp_f = (del_lim * del_lim + 2.0 * del_lim * del + eps2) / (del_lim * del_lim + del_lim * del + 2.0 * del * del + eps2);
                    alp = alp.min(alp_f);
                }
            }

            // store data
            alp_vec[vid_v] = alp;
        }

        // step 2: matrix computation

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let list_n = &tri.ngh_vid[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let offc = tri.vert_offc[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // get centroid properties
            let temp_o = temp_dom.correct_vert(offc, vid_v);

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let dist = dist_n[loc];
                let face = face_n[loc];
                let (nx, ny) = norm_n[loc];
                let skew = skew_n[loc];

                // get face properties
                let (velx_f, vely_f) = vel_dom.interp_face(skew, vid_v, vid_n);
                let den_f = den_dom.interp_face(skew, vid_v, vid_n);
                let hcap_f = hcap_dom.interp_face(skew, vid_v, vid_n);
                let flux = den_f * hcap_f * (velx_f * nx + vely_f * ny) * face;

                // load implicit upwind term into matrix
                // sparsity pattern must always be the same
                let coeff_v = flux.min(0.0);
                self.add_sclscl_a(vars, a_triplet, self.temp_id, vid_v, self.temp_id, vid_n, -coeff_v)?;
                self.add_sclscl_a(vars, a_triplet, self.temp_id, vid_v, self.temp_id, vid_v,  coeff_v)?;
                self.add_scl_b(vars, b_vec, self.temp_id, vid_v, -flux * temp_o)?;

                // get limited linear upwind face value
                let (vid_u, sign) = if flux > 0.0 {(vid_v, 1.0)} else {(vid_n, -1.0)};
                let (sx, sy) = skew;
                let (rx, ry) = (sign * 0.5 * dist * nx + sx, sign * 0.5 * dist * ny + sy);
                let (gradx_u, grady_u) = temp_dom.grad[vid_u];
                let alp_u = alp_vec[vid_u];
                let temp_corr = alp_u * (gradx_u * rx + grady_u * ry);

                // load deferred correction
                self.add_scl_b(vars, b_vec, self.temp_id, vid_v, flux * temp_corr)?;
            }
        }

        Ok(())
    }
}
