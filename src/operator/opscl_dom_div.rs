use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpSclDomDivergence {
    // domain
    pub dom_id: usize,

    // unknowns
    pub pres_id: usize,  // pressure (p)
    
    // others
    pub vel_id: usize,  // velocity (vel)
    pub den_id: usize,  // density (rho)
}

impl OpSclDomDivergence {
    pub fn new(vars: &Variables, dom_id: usize, pres_id: usize, vel_id: usize, den_id: usize) -> Result<OpSclDomDivergence, VFVChemError> {
        // adds the divergence term to the continuity equation
        // d(rho)/dt = -div(rho * vel)
        //             ^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, pres_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        
        // create struct
        let mut oper = OpSclDomDivergence::default();
        oper.dom_id = dom_id;
        oper.pres_id = pres_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        
        Ok(oper)
    }
}

impl OperatorBase for OpSclDomDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf -div(rho * vel) dA
        //
        // apply divergence theorem
        // int_lin -rho * vel . n dl
        //
        // discretize into sum over face values
        // sum_ngh -rho_f * (vel . n)_f * dl

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        
        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let face_n = &tri.face_len[vid_v];
            let list_n = &tri.ngh_vid[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let (nx, ny) = norm_n[loc];
                let skew = skew_n[loc];

                // get face propeties
                let (velx_s, vely_s) = vel_dom.correct_face(skew, vid_v, vid_n);
                let den_f = den_dom.interp_face(skew, vid_v, vid_n);

                // load into matrix
                let coeff_v = -0.5 * den_f * face;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_n, coeff_v * nx)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_v, coeff_v * nx)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_n, coeff_v * ny)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_v, coeff_v * ny)?;
                self.add_scl_b(vars, b_vec, self.pres_id, vid_v, den_f * (velx_s * nx + vely_s * ny) * face)?;
            }

        }

        Ok(())
    }
}
