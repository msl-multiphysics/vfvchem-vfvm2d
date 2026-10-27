use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::operator::extend::pore_helper;

#[derive(Default)]
pub struct PoreOpVecItfTranspose {
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

impl PoreOpVecItfTranspose {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, visc_a_id: usize, epsl_a_id: usize, vel_b_id: usize, visc_b_id: usize, epsl_b_id: usize) -> Result<PoreOpVecItfTranspose, VFVChemError> {
        // adds the transpose diffusion term to the combined control volumes of interface nodes
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        //                                                                               ^^^^^^^^^^^^^^              ^^^^^^^^^^^
        // each domain contributes its half of the control volume; loads into the rows of domain b
        // each half is integrated over its closed contour as in PoreOpVecBndTranspose, then the interface flux is removed
        // the one-sided vertex gradient cancels on a closed contour, which keeps the deferred iteration stable

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, visc_a_id, visc_b_id)?;
        check_itf_scldom(vars, itf_id, epsl_a_id, epsl_b_id)?;

        // create struct
        let mut oper = PoreOpVecItfTranspose::default();
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

impl OperatorBase for PoreOpVecItfTranspose {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf div( mu_eff * grad(vel)^T ) dA
        // int_lin mu_eff * grad(vel)^T . n dl
        // sum_ngh mu_eff * grad(vel)^T . n dl + sum_bnd mu_eff * grad(vel)^T . n dl - sum_itf mu_eff * grad(vel)^T . n dl
        // sum_ngh mu_eff_f * grad(vel)^T_f . n * dl + sum_bnd mu_eff_b * grad(vel)^T_b . n * dl - sum_itf mu_eff_b * grad(vel)^T_i . n * dl
        //
        // apply face skewness correction
        // mu_eff_f = (mu_eff_v + mu_eff_n) / 2 + (grad(mu_eff_v) + grad(mu_eff_n)) / 2 . skew
        //
        // apply half-edge interpolation
        // mu_eff_b = 0.75 * mu_eff_v + 0.25 * mu_eff_n
        //
        // use differences along the interface with tangent t
        // grad(vel)^T_i . n = t * d(vel . n)/dt + n * (div(vel)_b - d(vel . t)/dt)

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, visc_id, epsl_id) in [(self.vel_a_id, self.visc_a_id, self.epsl_a_id), (self.vel_b_id, self.visc_b_id, self.epsl_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];
            let vel_dom = &vars.vecdom[vel_id];

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
                let (u_v, v_v) = vel_dom.value[vid_v];
                let (visc_v, (visc_x_v, visc_y_v)) = pore_helper::visc_eff(vars, visc_id, epsl_id, vid_v);
                let (u_x_v, u_y_v, v_x_v, v_y_v) = vel_dom.grad[vid_v];

                // iterate through faces
                let mut sumx = 0.0;
                let mut sumy = 0.0;
                for (loc, &vid_n) in list_n.iter().enumerate() {
                    // get geometry
                    let face = face_n[loc];
                    let (nx, ny) = norm_n[loc];
                    let (sx, sy) = skew_n[loc];

                    // get neighbor properties
                    let (visc_n, (visc_x_n, visc_y_n)) = pore_helper::visc_eff(vars, visc_id, epsl_id, vid_n);
                    let (u_x_n, u_y_n, v_x_n, v_y_n) = vel_dom.grad[vid_n];

                    // get face viscosity
                    let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                    let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n));  // skew correction
                    let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity

                    // get face velocity gradient; average of both sides
                    let u_x_f = 0.5 * (u_x_v + u_x_n);
                    let u_y_f = 0.5 * (u_y_v + u_y_n);
                    let v_x_f = 0.5 * (v_x_v + v_x_n);
                    let v_y_f = 0.5 * (v_y_v + v_y_n);

                    // add to sum
                    sumx += visc_f * (u_x_f * nx + v_x_f * ny) * face;
                    sumy += visc_f * (u_y_f * nx + v_y_f * ny) * face;
                }

                // iterate through boundary half-edges
                for (blc, &vid_n) in list_b.iter().enumerate() {
                    // get geometry
                    let half = half_b[blc];
                    let (nx, ny) = norm_b[blc];

                    // get neighbor properties
                    let (visc_n, _) = pore_helper::visc_eff(vars, visc_id, epsl_id, vid_n);
                    let (u_x_n, u_y_n, v_x_n, v_y_n) = vel_dom.grad[vid_n];

                    // get half-edge viscosity
                    let visc_b = 0.75 * visc_v + 0.25 * visc_n;

                    // get half-edge velocity gradient
                    let u_x_b = 0.75 * u_x_v + 0.25 * u_x_n;
                    let u_y_b = 0.75 * u_y_v + 0.25 * u_y_n;
                    let v_x_b = 0.75 * v_x_v + 0.25 * v_x_n;
                    let v_y_b = 0.75 * v_y_v + 0.25 * v_y_n;

                    // add to sum
                    sumx += visc_b * (u_x_b * nx + v_x_b * ny) * half;
                    sumy += visc_b * (u_y_b * nx + v_y_b * ny) * half;

                    // remove the interface flux
                    let gnd_n = dom.nds_gnd[vid_n - dom.num_pts()];
                    if !itf.nds_gnd.contains(&gnd_n) { continue; }

                    // get tangent along the interface
                    let (x_v, y_v) = vars.nds.nds_pos[gnd_v];
                    let (x_n, y_n) = vars.nds.nds_pos[gnd_n];
                    let dist = 2.0 * half;
                    let (tx, ty) = ((x_n - x_v) / dist, (y_n - y_v) / dist);

                    // get tangential derivatives and half-edge divergence
                    let (u_n, v_n) = vel_dom.value[vid_n];
                    let dnorm = ((u_n - u_v) * nx + (v_n - v_v) * ny) / dist;
                    let dtang = ((u_n - u_v) * tx + (v_n - v_v) * ty) / dist;
                    let div_b = u_x_b + v_y_b;

                    // subtract from sum
                    sumx -= visc_b * (tx * dnorm + nx * (div_b - dtang)) * half;
                    sumy -= visc_b * (ty * dnorm + ny * (div_b - dtang)) * half;
                }

                // load into vector
                self.add_vec_b(vars, b_vec, self.vel_b_id, 0, vid_b_v, -sumx)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 1, vid_b_v, -sumy)?;

            }
        }

        Ok(())
    }
}
