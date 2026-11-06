use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpVecDomTranspose {
    // domain
    pub dom_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub visc_id: usize,  // viscosity (mu)
}

impl OpVecDomTranspose {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, visc_id: usize) -> Result<OpVecDomTranspose, VFVChemError> {
        // adds the advection term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                      ^^^^^^^^^               ^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;

        // create struct
        let mut oper = OpVecDomTranspose::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecDomTranspose {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf div( mu * grad(vel)^T ) dA
        //
        // apply divergence theorem
        // int_lin mu * grad(vel)^T . n dl
        //
        // discretize into sum over face values
        // sum_ngh mu_f * (grad(vel)^T . n)_f dl
        //
        // approximate face gradient as average of vertex gradients
        // do not apply corrections to avoid noise from higher-order terms
        // (grad(vel)^T . n)_f = 0.5 * ((grad(vel)^T . n)_v + (grad(vel)^T . n)_n)

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let visc_dom = &vars.scldom[self.visc_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let face_n = &tri.face_len[vid_v];
            let list_n = &tri.ngh_vid[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // get vertex properties
            let (velx_x_v, velx_y_v, vely_x_v, vely_y_v) = vel_dom.grad[vid_v];

            // iterate through faces
            let mut sumx = 0.0;
            let mut sumy = 0.0;
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let (nx, ny) = norm_n[loc];
                let skew = skew_n[loc];

                // get neighbor properties
                let (velx_x_n, velx_y_n, vely_x_n, vely_y_n) = vel_dom.grad[vid_n];

                // get face properties
                let velx_x_f = 0.5 * (velx_x_v + velx_x_n);
                let velx_y_f = 0.5 * (velx_y_v + velx_y_n);
                let vely_x_f = 0.5 * (vely_x_v + vely_x_n);
                let vely_y_f = 0.5 * (vely_y_v + vely_y_n);
                let visc_f = visc_dom.interp_face_harmonic(skew, vid_v, vid_n);

                // add to sum
                sumx += visc_f * (velx_x_f * nx + vely_x_f * ny) * face;
                sumy += visc_f * (velx_y_f * nx + vely_y_f * ny) * face;
            }

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -sumx)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -sumy)?;
        }

        Ok(())
    }
}
