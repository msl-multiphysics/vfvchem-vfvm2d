use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct SclBndOutflow {
    // boundary
    pub bnd_id: usize,
}

impl SclBndOutflow {
    pub fn new(_vars: &Variables, bnd_id: usize) -> Result<SclBndOutflow, VFVChemError> {
        // sets the outflow boundary condition

        // create struct
        let mut oper = SclBndOutflow::default();
        oper.bnd_id = bnd_id;

        Ok(oper)
    }
}

impl OperatorBase for SclBndOutflow {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // this operator does nothing
        // normal diffusive flux is zero
        // advection leaves naturally
        
        Ok(())
    }
}
