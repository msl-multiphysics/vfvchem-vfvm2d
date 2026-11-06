use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpSclDomLinear {
    // domain
    pub dom_id: usize,

    // unknowns
    pub unk_id: usize,  // unknown scalar (c_i)

    // others
    pub lin_id: usize,  // source coefficient (L_i)
}

impl OpSclDomLinear {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, lin_id: usize) -> Result<OpSclDomLinear, VFVChemError> {
        // adds the linear source term to scalar transport equations
        // rho * d(c_i)/dt + rho * vel . grad(c_i) = -div( -D_ij * grad(c_j) ) + L_i c_i + R_i
        //                                                                       ^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, lin_id)?;

        // create struct
        let mut oper = OpSclDomLinear::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.lin_id = lin_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclDomLinear {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf L c dA
        //
        // discretize into centroid values
        // L_c * c_c * dA

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let lin_dom = &vars.scldom[self.lin_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let offc = tri.vert_offc[vid_v];

            // get centroid properties
            let lin_c = lin_dom.interp_vert(offc, vid_v);
            let unk_o = unk_dom.correct_vert(offc, vid_v);

            // load into matrix and vector
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, lin_c * area)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -lin_c * unk_o * area)?;
        }

        Ok(())
    }
}
