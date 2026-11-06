use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpVecBndDirichlet {
    // boundary
    pub bnd_id: usize,

    // unknowns
    pub unk_id: usize,  // velocity (vel_i)

    // others
    pub val_id: usize,  // value (vel_set)

    // dirichlet flag
    pub is_dir: Vec<bool>,  // xid -> true if dirichlet
}

impl OpVecBndDirichlet {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, val_id: usize) -> Result<OpVecBndDirichlet, VFVChemError> {
        // sets the boundary value
        // vel_i = vel_set

        // input checking
        check_bnd_vecdom(vars, bnd_id, unk_id)?;
        check_bnd_vecbnd(vars, bnd_id, val_id)?;
        
        // create struct
        let mut oper = OpVecBndDirichlet::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.val_id = val_id;
        
        Ok(oper)
    }

    pub fn apply_initial(&mut self, vars: &mut Variables, mat_size: usize) -> Result<(), VFVChemError> {
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];

        // initialize dirichlet flag
        self.is_dir = vec![false; mat_size];
        let start_u = match vars.vecdom[self.unk_id].vecdom_type {
            VecDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let start_v = start_u + dom.num_vert();

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // load into unknown
            vars.vecdom[self.unk_id].value[vid_v] = vars.vecbnd[self.val_id].value[idx_v];
        
            // set dirichlet flag
            self.is_dir[start_u + vid_v] = true;
            self.is_dir[start_v + vid_v] = true;
        }

        Ok(())
    }
}

impl OperatorBase for OpVecBndDirichlet {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let val_bnd = &vars.vecbnd[self.val_id];

        // remove existing equations along dirichlet rows
        a_triplet.retain(|&(row, _, _)| !self.is_dir[row]);

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;
            
            // get value
            let (valu, valv) = val_bnd.value[idx_v];

            // load into matrix
            self.add_vecvec_a(vars, a_triplet, self.unk_id, 0, vid_v, self.unk_id, 0, vid_v, 1.0)?;
            self.add_vecvec_a(vars, a_triplet, self.unk_id, 1, vid_v, self.unk_id, 1, vid_v, 1.0)?;
            self.set_vec_b(vars, b_vec, self.unk_id, 0, vid_v, valu)?;
            self.set_vec_b(vars, b_vec, self.unk_id, 1, vid_v, valv)?;
        }

        Ok(())
    }
}
