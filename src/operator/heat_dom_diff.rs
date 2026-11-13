use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct HeatDomDiffusion {
    // domain
    pub dom_id: usize,

    // unknowns
    pub temp_id: usize,   // temperature (T)

    // others
    pub cond_id: usize,  // thermal conductivity (k)
}

impl HeatDomDiffusion {
    pub fn new(vars: &Variables, dom_id: usize, temp_id: usize, cond_id: usize) -> Result<HeatDomDiffusion, VFVChemError> {
        // adds the diffusion term to the heat transfer equation
        // rho * cp * d(T)/dt + rho * cp * vel . grad(T) = -div( -k * grad(T) ) + L * T + S
        //                                                 ^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, temp_id)?;
        check_dom_scldom(vars, dom_id, cond_id)?;

        // create struct
        let mut oper = HeatDomDiffusion::default();
        oper.dom_id = dom_id;
        oper.temp_id = temp_id;
        oper.cond_id = cond_id;

        Ok(oper)
    }
}

impl OperatorBase for HeatDomDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf -div( -k * grad(T) ) dA
        //
        // apply divergence theorem
        // int_lin k * grad(T) . n dl
        //
        // discretize into sum over face values
        // sum_ngh k_f * (grad(T) . n)_f dl
        //
        // for temperature gradient, use normal derivative
        // (grad(T) . n)_f = (T_n - T_v) / dn
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let cond_dom = &vars.scldom[self.cond_id];

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
                let cond_f = cond_dom.interp_face_harmonic(skew, vid_v, vid_n);

                // load into matrix
                let coeff_v = -cond_f * face / dist;
                self.add_sclscl_a(vars, a_triplet, self.temp_id, vid_v, self.temp_id, vid_n, -coeff_v)?;
                self.add_sclscl_a(vars, a_triplet, self.temp_id, vid_v, self.temp_id, vid_v,  coeff_v)?;
            }
        }
    
        Ok(())
    }
}
