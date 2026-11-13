use crate::base::prelude::*;
use crate::operator::prelude::*;

const VENKATAKRISHNAN_SCALE: f64 = 1.0;

#[derive(Default)]
pub struct MoleDomAdvection {
    // domain
    pub dom_id: usize,

    // unknowns
    pub conc_id: usize,   // concentration (c_i); mole solute per mass solvent

    // others
    pub vel_id: usize,  // velocity (vel)
    pub den_id: usize,  // density (rho)
}

impl MoleDomAdvection {
    pub fn new(vars: &Variables, dom_id: usize, conc_id: usize, vel_id: usize, den_id: usize) -> Result<MoleDomAdvection, VFVChemError> {
        // adds the advection term to the mass transfer equation
        // rho * d(c_i)/dt + rho * vel . grad(c_i) = -div( -rho * D_ij * grad(c_j) ) + L_i * c_i + R_i
        //                   ^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, conc_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = MoleDomAdvection::default();
        oper.dom_id = dom_id;
        oper.conc_id = conc_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;

        Ok(oper)
    }
}

impl OperatorBase for MoleDomAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the RHS of the PDE
        // -rho * vel . grad(c)
        //
        // rewrite advection term
        // -rho * vel . grad(c)
        // -div(rho * vel * c) + c * div(rho * vel)
        //
        // get surface integral of 1st term
        // int_surf -div(rho * vel * c) dA
        // int_lin -(rho * vel * c) . n dl
        // sum_ngh -(rho * vel * c)_f . n dl
        //
        // define mass flux
        // let F = rho * vel . n * dl
        // sum_ngh (-F_f * c_f) -> 1st term
        //
        // get surface integral of 2nd term
        // int_surf c * div(rho * vel) dA
        // c_c * div(rho * vel)_c dA
        //
        // approximate divergence using divergence theorem
        // div(rho * vel)_c = (1/dA) sum_ngh (rho * vel)_f . n dl
        // c_c * sum_ngh (rho * vel)_f . n dl
        // c_c * sum_ngh F_f
        // 
        // apply off-center correction
        // c_c = c_v + c_o 
        // (c_v + c_o) * sum_ngh F_f -> 2nd term
        //
        // add the terms
        // sum_ngh -F_f * c_f + (c_v + c_o) * F_f
        // sum_ngh -F_f * (c_f - c_v) + F_f * c_o
        //
        // face concentration c_f is computed using venkatakrishnan limiter
        // let c_min be the minimum c_v or c_n among neighbors
        // let c_max be the maximum c_v or c_n among neighbors
        // del_f = grad(c)_v . (x_f - x_v)
        // del_lim = (c_max - c_v) if del_f > 0 else (c_min - c_v) if del_f < 0
        // eps^2 = (K * h / L)^3 * del_c^2; h = sqrt(dA), L = sqrt(sum dA), del_c = range of c in domain
        // alp_f = (del_lim^2 + 2 * del_lim * del_f + eps^2) / (del_lim^2 + del_lim * del_f + 2 * del_f^2 + eps^2)
        // alp_v = min(1, alp_f among faces)
        //
        // let c_u be the upwind vertex concentration (either c_v or c_n)
        // c_f = c_u + alp_u * grad(c)_u . (x_f - x_u)

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let conc_dom = &vars.scldom[self.conc_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];

        // step 1: limiter computation

        // initialize
        let mut alp_vec = vec![1.0; dom.num_vert()];

        // get limiter scales
        let (mut conc_min, mut conc_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for &conc in conc_dom.value.iter() {
            conc_min = conc_min.min(conc);
            conc_max = conc_max.max(conc);
        }
        let conc2_scale = (conc_max - conc_min) * (conc_max - conc_min);
        let len_scale = tri.vert_area.iter().sum::<f64>().sqrt();

        // compute limiter at each vertex
        for vid_v in 0..dom.num_vert() {
            // get bounds from neighbors
            let conc_v = conc_dom.value[vid_v];
            let (gradx_v, grady_v) = conc_dom.grad[vid_v];
            let (mut conc_min, mut conc_max) = (conc_v, conc_v);
            for &vid_n in tri.ngh_vid[vid_v].iter() {
                let conc_n = conc_dom.value[vid_n];
                conc_min = conc_min.min(conc_n);
                conc_max = conc_max.max(conc_n);
            }

            // get limiter tolerance
            let eps2 = (VENKATAKRISHNAN_SCALE * tri.vert_area[vid_v].sqrt() / len_scale).powi(3) * conc2_scale;

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
                    let del_lim = if del > 0.0 {conc_max - conc_v} else {conc_min - conc_v};
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
            let conc_o = conc_dom.correct_vert(offc, vid_v);

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
                let flux = den_f * (velx_f * nx + vely_f * ny) * face;

                // load implicit upwind term into matrix
                // sparsity pattern must always be the same
                let coeff_v = flux.min(0.0);
                self.add_sclscl_a(vars, a_triplet, self.conc_id, vid_v, self.conc_id, vid_n, -coeff_v)?;
                self.add_sclscl_a(vars, a_triplet, self.conc_id, vid_v, self.conc_id, vid_v,  coeff_v)?;
                self.add_scl_b(vars, b_vec, self.conc_id, vid_v, -flux * conc_o)?;

                // get limited linear upwind face value
                let (vid_u, sign) = if flux > 0.0 {(vid_v, 1.0)} else {(vid_n, -1.0)};
                let (sx, sy) = skew;
                let (rx, ry) = (sign * 0.5 * dist * nx + sx, sign * 0.5 * dist * ny + sy);
                let (gradx_u, grady_u) = conc_dom.grad[vid_u];
                let alp_u = alp_vec[vid_u];
                let conc_corr = alp_u * (gradx_u * rx + grady_u * ry);

                // load deferred correction
                self.add_scl_b(vars, b_vec, self.conc_id, vid_v, flux * conc_corr)?;
            }
        }

        Ok(())
    }
}
