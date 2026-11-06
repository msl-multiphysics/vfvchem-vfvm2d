use crate::base::prelude::*;
use crate::operator::prelude::*;

const INVALID: usize = usize::MAX;

#[derive(Default)]
pub struct OpVecItfContinuity {
    // interface
    pub itf_id: usize,

    // unknowns
    pub unk_a_id: usize,  // unknown vector (vel_i^a)
    pub unk_b_id: usize,  // unknown vector (vel_i^b)

    // continuity rows
    pub row_move: Vec<usize>,  // row a -> row b
    pub row_pair: Vec<(usize, usize)>,  // (row a, row b)
}

impl OpVecItfContinuity {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, unk_b_id: usize) -> Result<OpVecItfContinuity, VFVChemError> {
        // sets a continuous value and flux across an interface for each component
        // vel_i^a = vel_i^b
        // R^a + R^b = 0

        // input checking
        check_itf_vecdom(vars, itf_id, unk_a_id, unk_b_id)?;

        // create struct
        let mut oper = OpVecItfContinuity::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.unk_b_id = unk_b_id;

        Ok(oper)
    }

    pub fn apply_initial(&mut self, vars: &mut Variables, mat_size: usize) -> Result<(), VFVChemError> {
        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_a = &vars.dom[vars.vecdom[self.unk_a_id].dom_id];
        let dom_b = &vars.dom[vars.vecdom[self.unk_b_id].dom_id];

        // get start indices
        let start_ua = match vars.vecdom[self.unk_a_id].vecdom_type {
            VecDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let start_ub = match vars.vecdom[self.unk_b_id].vecdom_type {
            VecDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let start_va = start_ua + dom_a.num_vert();
        let start_vb = start_ub + dom_b.num_vert();

        // initialize continuity rows
        self.row_move = vec![INVALID; mat_size];
        self.row_pair.clear();
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_a_v = dom_a.num_pts() + dom_a.nds_nid[&gnd_v];
            let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

            // set continuity rows
            for (start_a, start_b) in [(start_ua, start_ub), (start_va, start_vb)] {
                let row_a = start_a + vid_a_v;
                let row_b = start_b + vid_b_v;
                self.row_move[row_a] = row_b;
                self.row_pair.push((row_a, row_b));
            }
        }

        Ok(())
    }
}

impl OperatorBase for OpVecItfContinuity {
    fn apply(&self, _vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // outline:
        // add entries in row a to row b (flux continuity)
        // clear row a
        // set value continuity in row a

        // add entries in row a to row b (matrix)
        for (r, _, _) in a_triplet.iter_mut() {
            let row_b = self.row_move[*r];
            if row_b != INVALID { *r = row_b; }
        }

        // add entries in row a to row b (vector)
        for &(row_a, row_b) in self.row_pair.iter() {
            b_vec[row_b] += b_vec[row_a];
            b_vec[row_a] = 0.0;
        }

        // set value continuity in row a
        for &(row_a, row_b) in self.row_pair.iter() {
            a_triplet.push((row_a, row_a,  1.0));
            a_triplet.push((row_a, row_b, -1.0));
        }

        Ok(())
    }
}
