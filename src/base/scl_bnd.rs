use crate::base::prelude::*;

// input: time, (x, y), slice of unknown values at a point
// output: function value at the point
pub type SclBndFunc = Box<dyn Fn(f64, (f64, f64), &[f64]) -> f64 + Send + Sync>;

#[derive(Default)]
pub enum SclBndType {
    #[default] Constant,
    Function {
        func: SclBndFunc,  
        unkn_id: Vec<usize>,
    },
}

#[derive(Default)]
pub struct ScalarBoundary {
    // ids
    pub sclbnd_id: usize,
    pub bnd_id: usize,  // boundary this is attached to

    // node data
    pub sclbnd_type: SclBndType,
    pub value: Vec<f64>,  // idx -> value
}

impl ScalarBoundary {
    pub fn new_cnst(sclbnd_id: usize, bnd: &Boundary, value: f64) -> Result<ScalarBoundary, VFVChemError> {
        // default class
        let mut sclbnd = ScalarBoundary::default();
        sclbnd.sclbnd_id = sclbnd_id;
        sclbnd.bnd_id = bnd.bnd_id;

        // node data
        sclbnd.sclbnd_type = SclBndType::Constant;
        sclbnd.value = vec![value; bnd.num_nds()];

        Ok(sclbnd)
    }

    pub fn new_func(sclbnd_id: usize, bnd: &Boundary, scldom_all: &Vec<ScalarDomain>, func: SclBndFunc, unkn_id: Vec<usize>) -> Result<ScalarBoundary, VFVChemError> {
        // check that unkn_ids refer to unknown scalars
        for &i in unkn_id.iter() {
            match scldom_all[i].scldom_type {
                SclDomType::Unknown { start: _ } => { continue; }
                _ => { return Err(VFVChemError::ScalarBoundaryFunctionDependsOnNonUnknown { scldom_id: i }) }
            }
        }
        
        // default class
        let mut sclbnd = ScalarBoundary::default();
        sclbnd.sclbnd_id = sclbnd_id;
        sclbnd.bnd_id = bnd.bnd_id;

        // node data
        sclbnd.sclbnd_type = SclBndType::Function { func, unkn_id };
        sclbnd.value = vec![f64::NAN; bnd.num_nds()];  // placeholder

        Ok(sclbnd)
    }

    pub fn new_cnst_nonuniform(sclbnd_id: usize, nds: &Nodes, bnd: &Boundary, value_func: SclPosFunc) -> Result<ScalarBoundary, VFVChemError> {
        // default class
        let mut sclbnd = ScalarBoundary::default();
        sclbnd.sclbnd_id = sclbnd_id;
        sclbnd.bnd_id = bnd.bnd_id;

        // compute values
        let mut value = Vec::with_capacity(bnd.num_nds());  // idx -> value
        for &gnd in bnd.nds_gnd.iter() {
            let (x, y) = nds.nds_pos[gnd];
            value.push(value_func(x, y));
        }

        // node data
        sclbnd.sclbnd_type = SclBndType::Constant;
        sclbnd.value = value;

        Ok(sclbnd)
    }

    pub fn update_func(vars: &mut Variables, sclbnd_id: usize, t: f64) {
        // does not use &mut self for similarity with scldom's update_func
        // assume that scldoms in unkn_id have the same dom as sclbnd_id
        
        // pseudocode:
        // collect values from unknown-type scldom
        // update the selected function-type sclbnd

        // compute updated values
        let sclbnd = &vars.sclbnd[sclbnd_id];
        match &sclbnd.sclbnd_type {
            SclBndType::Function { func, unkn_id } => {
                // get objects
                let bnd = &vars.bnd[sclbnd.bnd_id];
                let dom = &vars.dom[bnd.dom_id];

                // collect positions and values from unknown scalars
                let mut pos = Vec::with_capacity(bnd.num_nds());  // idx -> (x, y)
                let mut value_u = Vec::with_capacity(bnd.num_nds());  // idx -> values of unknowns
                for idx in 0..bnd.num_nds() {
                    // get domain index
                    let gnd = bnd.nds_gnd[idx];
                    let nid = dom.nds_nid[&gnd];
                    let vid = dom.num_pts() + nid;

                    // collect position and values
                    pos.push(vars.nds.nds_pos[gnd]);
                    let mut value_sub = Vec::with_capacity(unkn_id.len());  // values of unknowns
                    for &i in unkn_id {
                        let scldom = &vars.scldom[i];
                        value_sub.push(scldom.value[vid]);
                    }
                    value_u.push(value_sub);
                }

                // compute function value
                let mut value_f = Vec::with_capacity(bnd.num_nds());
                for idx in 0..bnd.num_nds() {
                    let value_sub = func(t, pos[idx], value_u[idx].as_slice());
                    value_f.push(value_sub)
                }

                // modify scldom values
                let sclbnd = &mut vars.sclbnd[sclbnd_id];
                sclbnd.value = value_f;

            }
            _ => { return; }  // skip if not function-type
        };

    }

}
