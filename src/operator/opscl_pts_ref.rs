use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpSclPtsReference {
    // domain
    pub dom_id: usize,

    // unknowns
    pub unk_id: usize,  // unknown scalar (c_i)

    // others
    pub nid: usize,
    pub value: f64,
}

impl OpSclPtsReference {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, nid: usize, value: f64) -> Result<OpSclPtsReference, VFVChemError> {
        // sets the boundary value
        // c_i = c_set

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;

        // create struct
        let mut oper = OpSclPtsReference::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.nid = nid;
        oper.value = value;
        
        Ok(oper)
    }
}

impl OperatorBase for OpSclPtsReference {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get objects
        let dom = &vars.dom[self.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        
        // get solution vector index
        let start = match unk_dom.scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let xid = start + dom.num_pts() + self.nid;

        // replace the existing equation with referenc
        a_triplet.retain(|tri| tri.0 != xid);
        a_triplet.push((xid, xid, 1.0));
        b_vec[xid] = self.value;
    
        Ok(())
    }
}
