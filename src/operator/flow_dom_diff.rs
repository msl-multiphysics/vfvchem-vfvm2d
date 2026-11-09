use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct FlowDomDiffusion {
    // domain
    pub dom_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub visc_id: usize,  // viscosity (mu)
}

impl FlowDomDiffusion {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, visc_id: usize) -> Result<FlowDomDiffusion, VFVChemError> {
        // adds the advection term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                      ^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;

        // create struct
        let mut oper = FlowDomDiffusion::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for FlowDomDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf div( mu * grad(vel) ) dA
        //
        // apply divergence theorem
        // int_lin mu * grad(vel) . n dl
        //
        // discretize into sum over face values
        // sum_ngh mu_f * (grad(vel) . n)_f dl
        //
        // for velocity gradient, use normal derivative
        // (grad(vel) . n)_f = (vel_n - vel_v) / dn

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let visc_dom = &vars.scldom[self.visc_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let list_n = &tri.ngh_vid[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let dist = dist_n[loc];
                let face = face_n[loc];
                let skew = skew_n[loc];

                // get face properties
                let visc_f = visc_dom.interp_face_harmonic(skew, vid_v, vid_n);

                // load into matrix
                let coeff_v = -visc_f * face / dist;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_n, -coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v,  coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_n, -coeff_v)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v,  coeff_v)?;
            }
        }

        Ok(())
    }
}
