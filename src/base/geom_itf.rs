use crate::base::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct Interface {
    // ids
    pub itf_id: usize,
    pub dom_a_id: usize,
    pub dom_b_id: usize,
    pub grp1d_id: String,

    // nodes on interface
    pub nds_gnd: Vec<usize>,  // idx (index) -> gnd; local to global
    pub nds_nid: HashMap<usize, usize>,  // gnd -> idx; global to local
}

impl Interface {
    pub fn new(itf_id: usize, nds: &Nodes, dom_a: &Domain, dom_b: &Domain, grp1d_id: &str) -> Result<Interface, VFVChemError> {
        // input checking
        if !nds.grp1d_eid.contains_key(grp1d_id) { return Err(VFVChemError::InterfaceGroupIDNotFound { grp1d_id: grp1d_id.to_string() })}

        // check that grp1d_id nodes are a subset of dom_a and dom_b
        let grp1d_eid = &nds.grp1d_eid[grp1d_id];
        let not_subset = !grp1d_eid.iter().flat_map(|&eid| nds.edge_gnd[eid].iter()).all(|gnd| {dom_a.nds_nid.contains_key(gnd) && dom_b.nds_nid.contains_key(gnd)});
        if not_subset { return Err(VFVChemError::InterfaceNoCommonInterface { dom_a: dom_a.dom_id, dom_b: dom_b.dom_id, grp1d_id: grp1d_id.to_string() })}

        // default struct
        let mut itf = Interface::default();
        itf.itf_id = itf_id;
        itf.dom_a_id = dom_a.dom_id;
        itf.dom_b_id = dom_b.dom_id;
        itf.grp1d_id = grp1d_id.to_string();

        // get nodes
        for &eid in nds.grp1d_eid[grp1d_id].iter() {
        for &gnd in nds.edge_gnd[eid].iter() {
            itf.nds_gnd.push(gnd);
        }}
        itf.nds_nid = dom_a.nds_gnd.iter().enumerate().map(|(idx, &gnd)| (gnd, idx)).collect();

        Ok(itf)
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_gnd.len()
    }
}
