use crate::base::prelude::*;

// input: time, (x, y), slice of unknown values at a point
// output: function value at the point
pub type VecBndFunc = Box<dyn Fn(f64, (f64, f64), &[f64]) -> (f64, f64) + Send + Sync>;

#[derive(Default)]
pub enum VecBndType {
    #[default] Constant,
    Function {
        func: VecBndFunc,  
        unkn_id: Vec<usize>,
    },
}

#[derive(Default)]
pub struct VectorBoundary {
    // ids
    pub vecbnd_id: usize,
    pub bnd_id: usize,  // boundary this is attached to

    // node data
    pub vecbnd_type: VecBndType,
    pub value: Vec<(f64, f64)>,  // idx -> value
    pub value_prev: Vec<(f64, f64)>,
}

impl VectorBoundary {
    pub fn new_cnst(vecbnd_id: usize, bnd: &Boundary, value: (f64, f64)) -> Result<VectorBoundary, VFVChemError> {
        // default class
        let mut vecbnd = VectorBoundary::default();
        vecbnd.vecbnd_id = vecbnd_id;
        vecbnd.bnd_id = bnd.bnd_id;

        // node data
        vecbnd.vecbnd_type = VecBndType::Constant;
        vecbnd.value = vec![value; bnd.num_nds()];
        vecbnd.value_prev = vecbnd.value.clone();
        
        Ok(vecbnd)
    }

    pub fn new_func(vecbnd_id: usize, bnd: &Boundary, scldom_all: &Vec<ScalarDomain>, func: VecBndFunc, unkn_id: Vec<usize>) -> Result<VectorBoundary, VFVChemError> {
        // check that unkn_ids refer to unknown scalars
        for &i in unkn_id.iter() {
            match scldom_all[i].scldom_type {
                SclDomType::Unknown { start: _ } => { continue; }
                _ => { return Err(VFVChemError::VectorBoundaryFunctionDependsOnNonUnknown { scldom_id: i }) }
            }
        }
        
        // default class
        let mut vecbnd = VectorBoundary::default();
        vecbnd.vecbnd_id = vecbnd_id;
        vecbnd.bnd_id = bnd.bnd_id;

        // node data
        vecbnd.vecbnd_type = VecBndType::Function { func, unkn_id };
        vecbnd.value = vec![(f64::NAN, f64::NAN); bnd.num_nds()];  // placeholder
        vecbnd.value_prev = vecbnd.value.clone();

        Ok(vecbnd)
    }

    pub fn update_prev(&mut self) {
        self.value_prev = self.value.clone();
    }

    pub fn update_func(vars: &mut Variables, vecbnd_id: usize, t: f64) {
        // does not use &mut self for similarity with scldom's update_func
        // assume that scldoms in unkn_id have the same dom as vecbnd_id
        
        // pseudocode:
        // collect values from unknown-type scldom
        // update the selected function-type vecbnd

        // compute updated values
        let vecbnd = &vars.vecbnd[vecbnd_id];
        match &vecbnd.vecbnd_type {
            VecBndType::Function { func, unkn_id } => {
                // get objects
                let bnd = &vars.bnd[vecbnd.bnd_id];
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
                let vecbnd = &mut vars.vecbnd[vecbnd_id];
                vecbnd.value = value_f;

            }
            _ => { return; }  // skip if not function-type
        };

    }

}
