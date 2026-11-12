use crate::base::prelude::*;
use crate::operator::prelude::*;

const VENKATAKRISHNAN_SCALE: f64 = 1.0;

#[derive(Default)]
pub struct FlowDomAdvection {
    // domain
    pub dom_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub den_id: usize,  // density (rho)
}

impl FlowDomAdvection {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, den_id: usize, _: usize) -> Result<FlowDomAdvection, VFVChemError> {
        // adds the advection term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                   ^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = FlowDomAdvection::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;

        Ok(oper)
    }
}

impl OperatorBase for FlowDomAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the RHS of the PDE
        // -rho * vel . grad(vel)
        //
        // rewrite advection term
        // -rho * vel . grad(vel)
        // -div(rho * vel * vel) + vel * div(rho * vel)
        // -div(rho * vel * vel) + vel * div(rho * vel)
        //
        // get surface integral of 1st term
        // int_surf -div(rho * vel * vel) dA
        // int_lin -(rho * vel * vel) . n dl
        // sum_ngh -(rho * vel * vel)_f . n dl
        //
        // define mass flux
        // let F = rho * vel . n * dl
        // sum_ngh (-F_f * vel_f) -> 1st term
        //
        // get surface integral of 2nd term
        // int_surf vel * div(rho * vel) dA
        // vel_c * div(rho * vel)_c dA
        //
        // approximate divergence using divergence theorem
        // div(rho * vel)_c = (1/dA) sum_ngh (rho * vel)_f . n dl
        // vel_c * sum_ngh (rho * vel)_f . n dl
        // vel_c * sum_ngh F_f
        // 
        // apply off-center correction
        // vel_c = vel_v + vel_o 
        // (vel_v + vel_o) * sum_ngh F_f -> 2nd term
        //
        // add the terms
        // sum_ngh -F_f * vel_f + (vel_v + vel_o) * F_f
        // sum_ngh -F_f * (vel_f - vel_v) + F_f * vel_o
        //
        // face velocity vel_f is computed using venkatakrishnan limiter
        // let vel_min be the minimum vel_v or vel_n among neighbors
        // let vel_max be the maximum vel_v or vel_n among neighbors
        // del_f = grad(vel)_v . (x_f - x_v)
        // del_lim = (vel_max - vel_v) if del_f > 0 else (vel_min - vel_v) if del_f < 0
        // eps^2 = (K * h / L)^3 * del_vel^2; h = sqrt(dA), L = sqrt(sum dA), del_vel = range of vel in domain
        // alp_f = (del_lim^2 + 2 * del_lim * del_f + eps^2) / (del_lim^2 + del_lim * del_f + 2 * del_f^2 + eps^2)
        // alp_v = min(1, alp_f among faces)
        //
        // let vel_u be the upwind vertex velocity (either vel_v or vel_n)
        // vel_f = vel_u + alp_u * grad(vel)_u . (x_f - x_u)
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];

        // step 1: limiter computation

        // initialize
        let mut alp_vec = vec![(1.0, 1.0); dom.num_vert()];

        // get limiter scales
        let (mut velx_min, mut velx_max, mut vely_min, mut vely_max) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        for &(velx, vely) in vel_dom.value.iter() {
            velx_min = velx_min.min(velx);
            velx_max = velx_max.max(velx);
            vely_min = vely_min.min(vely);
            vely_max = vely_max.max(vely);
        }
        let velx2_scale = (velx_max - velx_min) * (velx_max - velx_min);
        let vely2_scale = (vely_max - vely_min) * (vely_max - vely_min);
        let len_scale = tri.vert_area.iter().sum::<f64>().sqrt();

        // compute limiter at each vertex
        for vid_v in 0..dom.num_vert() {
            // get bounds from neighbors
            let (velx_v, vely_v) = vel_dom.value[vid_v];
            let (gradux_v, graduy_v, gradvx_v, gradvy_v) = vel_dom.grad[vid_v];
            let (mut velx_min, mut velx_max, mut vely_min, mut vely_max) = (velx_v, velx_v, vely_v, vely_v);
            for &vid_n in tri.ngh_vid[vid_v].iter() {
                let (velx_n, vely_n) = vel_dom.value[vid_n];
                velx_min = velx_min.min(velx_n);
                velx_max = velx_max.max(velx_n);
                vely_min = vely_min.min(vely_n);
                vely_max = vely_max.max(vely_n);
            }

            // get limiter tolerance
            let eps2_part = (VENKATAKRISHNAN_SCALE * tri.vert_area[vid_v].sqrt() / len_scale).powi(3);
            let epsx2 = eps2_part * velx2_scale;
            let epsy2 = eps2_part * vely2_scale;

            // iterate through faces
            let (mut alpx, mut alpy): (f64, f64) = (1.0, 1.0);
            for loc in 0..tri.ngh_vid[vid_v].len() {
                // get position of face center relative to vertex
                let dist = tri.ngh_dist[vid_v][loc];
                let (nx, ny) = tri.face_norm[vid_v][loc];
                let (sx, sy) = tri.face_skew[vid_v][loc];
                let (rx, ry) = (0.5 * dist * nx + sx, 0.5 * dist * ny + sy);  // x_f - x_v

                // get limiter for each component
                let delx = gradux_v * rx + graduy_v * ry;
                let dely = gradvx_v * rx + gradvy_v * ry;
                if delx != 0.0 {
                    let del_lim = if delx > 0.0 {velx_max - velx_v} else {velx_min - velx_v};
                    let alp_f = (del_lim * del_lim + 2.0 * del_lim * delx + epsx2) / (del_lim * del_lim + del_lim * delx + 2.0 * delx * delx + epsx2);
                    alpx = alpx.min(alp_f);
                }
                if dely != 0.0 {
                    let del_lim = if dely > 0.0 {vely_max - vely_v} else {vely_min - vely_v};
                    let alp_f = (del_lim * del_lim + 2.0 * del_lim * dely + epsy2) / (del_lim * del_lim + del_lim * dely + 2.0 * dely * dely + epsy2);
                    alpy = alpy.min(alp_f);
                }
            }

            // store data
            alp_vec[vid_v] = (alpx, alpy);
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
            let (velx_o, vely_o) = vel_dom.correct_vert(offc, vid_v);

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
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_n, -coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v,  coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_n, -coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v,  coeff_v)?;
                self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -flux * velx_o)?;
                self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -flux * vely_o)?;

                // get limited linear upwind face value
                let (vid_u, sign) = if flux > 0.0 {(vid_v, 1.0)} else {(vid_n, -1.0)};
                let (sx, sy) = skew;
                let (rx, ry) = (sign * 0.5 * dist * nx + sx, sign * 0.5 * dist * ny + sy);
                let (gradux_u, graduy_u, gradvx_u, gradvy_u) = vel_dom.grad[vid_u];
                let (alpx_u, alpy_u) = alp_vec[vid_u];
                let velx_corr = alpx_u * (gradux_u * rx + graduy_u * ry);
                let vely_corr = alpy_u * (gradvx_u * rx + gradvy_u * ry);

                // load deferred correction
                self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, flux * velx_corr)?;
                self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, flux * vely_corr)?;
            }
        }

        Ok(())
    }
}
