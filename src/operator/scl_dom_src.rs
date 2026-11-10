use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct SclDomSource {
    // domain
    pub dom_id: usize,

    // unknowns
    pub unk_id: usize,  // unknown scalar (c_i)

    // others
    pub src_id: usize,  // source (R_i)
}

impl SclDomSource {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, src_id: usize) -> Result<SclDomSource, VFVChemError> {
        // adds the source term to the scalar transport equation
        // rho * d(c_i)/dt + rho * vel . grad(c_i) = -div( -rho * D_ij * grad(c_j) ) + L_i * c_i + R_i
        //                                                                                         ^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, src_id)?;

        // create struct
        let mut oper = SclDomSource::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;
        
        Ok(oper)
    }
}

impl OperatorBase for SclDomSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf R dA
        //
        // discretize into centroid values
        // R_c * dA
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let src_dom = &vars.scldom[self.src_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let offc = tri.vert_offc[vid_v];
            
            // get centroid properties
            let src_c = src_dom.interp_vert(offc, vid_v);
            
            // load into matrix
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -src_c * area)?;
        }
    
        Ok(())
    }
}
