use crate::base::prelude::*;
use std::collections::HashMap;

// input: time
// output: velocity, origin, angular velocity
pub type BndVelFunc = Box<dyn Fn(f64) -> ((f64, f64), (f64, f64), f64) + Send + Sync>;

#[derive(Default)]
pub enum BndVelType {
    #[default] Constant,
    Function {
        func: BndVelFunc,
    }
}

#[derive(Default)]
pub struct Boundary {
    // ids
    pub bnd_id: usize,
    pub dom_id: usize,  // domain this boundary is attached to
    pub grp1d_id: String,
    
    // nodes on boundary
    pub nds_gnd: Vec<usize>,  // idx (index) -> gnd; local to global
    pub nds_nid: HashMap<usize, usize>,  // gnd -> idx; global to local
    pub nds_idx: HashMap<usize, usize>,  // gnd -> idx; only nodes on this boundary
    pub nds_norm: Vec<(f64, f64)>,  // idx -> (nx, ny)

    // boundary movement
    pub bndvel_type: BndVelType,
    pub vel: (f64, f64),
    pub org: (f64, f64),
    pub angvel: f64,
}

impl Boundary {
    pub fn new(bnd_id: usize, nds: &Nodes, dom: &Domain, grp1d_id: &str) -> Result<Boundary, VFVChemError> {
        // input checking
        if !nds.grp1d_gnd.contains_key(grp1d_id) { return Err(VFVChemError::BoundaryGroupIDNotFound { grp1d_id: grp1d_id.to_string() })}

        // default struct
        let mut bnd = Boundary::default();
        bnd.bnd_id = bnd_id;
        bnd.dom_id = dom.dom_id;
        bnd.grp1d_id = grp1d_id.to_string();

        // get nodes
        bnd.nds_gnd = nds.grp1d_gnd[grp1d_id].clone();
        bnd.nds_nid = dom.nds_gnd.iter().enumerate().map(|(idx, &gnd)| (gnd, idx)).collect();
        bnd.nds_idx = bnd.nds_gnd.iter().enumerate().map(|(idx, &gnd)| (gnd, idx)).collect();

        // get outward normals
        let nds_norm: Vec<(f64, f64)> = nds.grp2d_norm[&dom.grp2d_id].iter().flatten().copied().collect();
        bnd.nds_norm = bnd.nds_gnd.iter().map(|&gnd| nds_norm[bnd.nds_nid[&gnd]]).collect();

        Ok(bnd)
    }

    pub fn set_bndvel_cnst(&mut self, vel: (f64, f64), org: (f64, f64), angvel: f64) -> Result<(), VFVChemError> {
        self.bndvel_type = BndVelType::Constant;
        self.vel = vel;
        self.org = org;
        self.angvel = angvel;
        Ok(())
    }

    pub fn set_bndvel_func(&mut self, func: BndVelFunc) -> Result<(), VFVChemError> {
        self.bndvel_type = BndVelType::Function { func };
        Ok(())
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_gnd.len()
    }
}
