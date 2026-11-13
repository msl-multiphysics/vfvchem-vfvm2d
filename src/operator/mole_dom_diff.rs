use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct MoleDomDiffusion {
    // domain
    pub dom_id: usize,

    // unknowns
    pub conc_id: usize,   // concentration (c_i); mole solute per mass solvent

    // others
    pub den_id: usize,   // density (rho)
    pub diff_id: usize,  // diffusion coefficient (D_ij)
    pub drv_id: usize,   // driving concentration (c_j)
}

impl MoleDomDiffusion {
    pub fn new(vars: &Variables, dom_id: usize, conc_id: usize, den_id: usize, diff_id: usize, drv_id: usize) -> Result<MoleDomDiffusion, VFVChemError> {
        // adds the diffusion term to the mass transfer equation
        // rho * d(c_i)/dt + rho * vel . grad(c_i) = -div( -rho * D_ij * grad(c_j) ) + L_i * c_i + S_i
        //                                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, conc_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, diff_id)?;
        check_dom_scldom(vars, dom_id, drv_id)?;

        // create struct
        let mut oper = MoleDomDiffusion::default();
        oper.dom_id = dom_id;
        oper.conc_id = conc_id;
        oper.den_id = den_id;
        oper.diff_id = diff_id;
        oper.drv_id = drv_id;

        Ok(oper)
    }
}

impl OperatorBase for MoleDomDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf -div( -rho * D * grad(c) ) dA
        //
        // apply divergence theorem
        // int_lin rho * D * grad(c) . n dl
        //
        // discretize into sum over face values
        // sum_ngh rho_f * D_f * (grad(c) . n)_f dl
        //
        // for concentration gradient, use normal derivative
        // (grad(c) . n)_f = (c_n - c_v) / dn
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let den_dom = &vars.scldom[self.den_id];
        let diff_dom = &vars.scldom[self.diff_id];

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
                let den_f = den_dom.interp_face(skew, vid_v, vid_n);
                let diff_f = diff_dom.interp_face_harmonic(skew, vid_v, vid_n);

                // load into matrix
                let coeff_v = -den_f * diff_f * face / dist;
                self.add_sclscl_a(vars, a_triplet, self.conc_id, vid_v, self.drv_id, vid_n, -coeff_v)?;
                self.add_sclscl_a(vars, a_triplet, self.conc_id, vid_v, self.drv_id, vid_v,  coeff_v)?;
            }
        }
    
        Ok(())
    }
}
