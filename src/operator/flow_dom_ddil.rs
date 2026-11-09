use crate::base::prelude::*;
use crate::operator::prelude::*;

const NEG_FRAC_TWO_THREE: f64 = -2.0/3.0;

#[derive(Default)]
pub struct FlowDomDilatatory {
    // domain
    pub dom_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub visc_id: usize,  // viscosity (mu)
}

impl FlowDomDilatatory {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, visc_id: usize) -> Result<FlowDomDilatatory, VFVChemError> {
        // adds the advection term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                      ^^^^^^^^^                               ^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;

        // create struct
        let mut oper = FlowDomDilatatory::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for FlowDomDilatatory {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf div( lam * div(vel) * I ) dA
        //
        // apply divergence theorem
        // int_lin lam * div(vel) * I . n dl
        // int_lin lam * div(vel) * n dl
        // 
        // discretize into sum over face values
        // sum_ngh lam_f * div(vel)_f * n * dl
        //
        // approximate face divergence as average of vertex divergences
        // do not apply corrections to avoid noise from higher-order terms
        // div(vel)_f = 0.5 * (div(vel)_v + div(vel)_n)

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let visc_dom = &vars.scldom[self.visc_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let face_n = &tri.face_len[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // get vertex properties
            let (velx_x_v, _, _, vely_y_v) = vel_dom.grad[vid_v];
            let div_v = velx_x_v + vely_y_v;

            // iterate through faces
            let mut sumx = 0.0;
            let mut sumy = 0.0;
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let (nx, ny) = norm_n[loc];
                let skew = skew_n[loc];

                // get neighbor properties
                let (velx_x_n, _, _, vely_y_n) = vel_dom.grad[vid_n];
                let div_n = velx_x_n + vely_y_n;

                // get face properties
                let div_f = 0.5 * (div_v + div_n);
                let visc_f = visc_dom.interp_face_harmonic(skew, vid_v, vid_n);
                let lam_f= NEG_FRAC_TWO_THREE * visc_f;
                
                // add to sum
                let coeff = lam_f * div_f * face;
                sumx += coeff * nx;
                sumy += coeff * ny;
            }

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -sumx)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -sumy)?;
        }

        Ok(())
    }
}
