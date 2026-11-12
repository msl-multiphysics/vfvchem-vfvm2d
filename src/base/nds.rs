use crate::base::prelude::*;
use std::collections::HashMap;

const OFFSET_SPACE: f64 = 0.2;
const CORNER_SPACE: f64 = 5.0;
const CORNER_RATIO: f64 = 1.2;

#[derive(Default)]
pub struct Nodes {
    // nodes
    // nodes are on boundaries and interfaces
    pub nds_pos: Vec<(f64, f64)>,  // gnd (global node id)

    // edges
    // defines the geometry
    // note: vertex here (edge endpoints) is different from vertex in domain object (points and nodes)
    pub vert_pos: Vec<(f64, f64)>,  // gvt (geometric vertex id)
    pub edge_gvt: Vec<[usize; 2]>,  // eid (edge id) -> gvt; different orders/windings are given unique edges
    pub edge_gnd: Vec<Vec<usize>>,  // eid -> [gnd, ...]

    // 1d groups
    // defines boundaries and interfaces
    pub grp1d_eid: HashMap<String, Vec<usize>>, // grp1d_id -> [eid, ...]

    // 2d groups
    // defines domains (ccw winding -> interior)
    pub grp2d_eid: HashMap<String, Vec<usize>>, // grp2d_id -> [eid, ...]

    // sampling
    pub nds_space: f64,  // node spacing
    pub off_space: f64,  // corner offset length to node spacing ratio
    pub crn_space: f64,  // corner refinement length to node spacing ratio
    pub crn_ratio: f64,  // geometric growth rate from corner offset length to node spacing
}

impl Nodes {
    pub fn new_obj(path: &str, nds_space: f64) -> Result<Nodes, VFVChemError> {
        Nodes::new_obj_options(path, nds_space, OFFSET_SPACE, CORNER_SPACE, CORNER_RATIO)
    }

    pub fn new_obj_options(path: &str, nds_space: f64, off_space: f64, crn_space: f64, crn_ratio: f64) -> Result<Nodes, VFVChemError> {
        // input checking
        if nds_space <= 0.0 { return Err(VFVChemError::NodesNonPositiveSpacing { nds_space }) }
        if off_space <= 0.0 || off_space > 1.0 { return Err(VFVChemError::NodesInvalidOffset { off_space }) }
        if crn_space <= 1.0 { return Err(VFVChemError::NodesInvalidCornerSpace { crn_space }) }
        if crn_ratio <= 1.0 { return Err(VFVChemError::NodesInvalidCornerRatio { crn_ratio }) }
        
        // create struct
        let nds = read_nds_obj(path, nds_space, off_space, crn_space, crn_ratio)?;
        Ok(nds)
    }

    pub fn new_bounds(x_min: f64, y_min: f64, x_max: f64, y_max: f64, nds_space: f64) -> Result<Nodes, VFVChemError> {
        Nodes::new_bounds_options(x_min, y_min, x_max, y_max, nds_space, OFFSET_SPACE, CORNER_SPACE, CORNER_RATIO)
    }
    
    pub fn new_bounds_options(x_min: f64, y_min: f64, x_max: f64, y_max: f64, nds_space: f64, off_space: f64, crn_space: f64, crn_ratio: f64) -> Result<Nodes, VFVChemError> {
        // create nodes spaced at least h apart around rectangle
        // no corner nodes are made; create nodes offset from corner by dh
        // no corner nodes prevents incompatible conditions and simplifies bookkeeping
        
        // input checking
        if x_max < x_min { return Err(VFVChemError::NodesInvalidBoundsX { x_min, x_max }) }
        if y_max < y_min { return Err(VFVChemError::NodesInvalidBoundsY { y_min, y_max }) }
        if nds_space <= 0.0 { return Err(VFVChemError::NodesNonPositiveSpacing { nds_space }) }
        if off_space <= 0.0 || off_space > 1.0 { return Err(VFVChemError::NodesInvalidOffset { off_space }) }
        if crn_space <= 1.0 { return Err(VFVChemError::NodesInvalidCornerSpace { crn_space }) }
        if crn_ratio <= 1.0 { return Err(VFVChemError::NodesInvalidCornerRatio { crn_ratio }) }

        // create struct
        let mut nds = Nodes::default();
        nds.nds_space = nds_space;
        nds.off_space = off_space;
        nds.crn_space = crn_space;
        nds.crn_ratio = crn_ratio;

        // construct nodes
        let nds_b = _nodes_line((x_min, y_min), (x_max, y_min), nds_space, off_space, crn_space, crn_ratio);  // bottom
        let nds_r = _nodes_line((x_max, y_min), (x_max, y_max), nds_space, off_space, crn_space, crn_ratio);  // right
        let nds_t = _nodes_line((x_max, y_max), (x_min, y_max), nds_space, off_space, crn_space, crn_ratio);  // top
        let nds_l = _nodes_line((x_min, y_max), (x_min, y_min), nds_space, off_space, crn_space, crn_ratio);  // left
        let num_x = nds_b.len();
        let num_y = nds_r.len();
        nds.nds_pos.extend(nds_b);
        nds.nds_pos.extend(nds_r);
        nds.nds_pos.extend(nds_t);
        nds.nds_pos.extend(nds_l);

        // construct edges
        nds.vert_pos = vec![(x_min, y_min), (x_max, y_min), (x_max, y_max), (x_min, y_max)];
        nds.edge_gvt = vec![[0, 1], [1, 2], [2, 3], [3, 0]];
        nds.edge_gnd = vec![
            (0..num_x).collect(),
            (num_x..num_x + num_y).collect(),
            (num_x + num_y..2*num_x + num_y).collect(),
            (2*num_x + num_y..2*num_x + 2*num_y).collect(),
        ];

        // 1d group nodes
        nds.grp1d_eid.insert("bnd_b".to_string(), vec![0]);  // bottom
        nds.grp1d_eid.insert("bnd_r".to_string(), vec![1]);  // right
        nds.grp1d_eid.insert("bnd_t".to_string(), vec![2]);  // top
        nds.grp1d_eid.insert("bnd_l".to_string(), vec![3]);  // left

        // 2d group nodes
        nds.grp2d_eid.insert("dom".to_string(), vec![0, 1, 2, 3]);

        Ok(nds)
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_pos.len()
    }
}

fn _nodes_line(pos0: (f64, f64), pos1: (f64, f64), nds_space: f64, off_space: f64, crn_space: f64, crn_ratio: f64) -> Vec<(f64, f64)> {
    // outline:
    // compute node spacing in 1D
    // place points on 2D line
    
    // compute node spacing in 1D

    // get endpoints of segment
    let (x0, y0) = pos0;
    let (x1, y1) = pos1;

    // get segment length
    let (len_x, len_y) = (x1 - x0, y1 - y0);
    let len_s = (len_x * len_x + len_y * len_y).sqrt();

    // compute corner refinement node positions
    // initial spacing is the corner offset
    let mut space_cnr = Vec::new();
    let a = off_space * nds_space;
    let b = crn_space * nds_space;
    let n = ((1.0 + (b - a) * (crn_ratio - 1.0) / a).ln() / crn_ratio.ln()) as i32;
    let clear_cnr = 0.5 * len_s - a;  // available space for refinement
    for i in 0..=n {
        let c = a + (b - a) * (crn_ratio.powi(i) - 1.0) / (crn_ratio.powi(n) - 1.0);
        if c > clear_cnr { break; }
        space_cnr.push(c);
    }
    let start_cnr = space_cnr[space_cnr.len() - 1];

    // compute regularly spaced node positions
    let clear_reg = len_s - 2.0 * b;  // amount of space for regular spacing
    let num_reg = (clear_reg / nds_space) as usize;  // number of regularly spaced intervals
    let ds = clear_reg / (num_reg as f64);
    let space_reg = (1..num_reg).map(|i| start_cnr + (i as f64) * ds);

    // place nodes along segment length
    // refinement -> regular -> refinement in reverse
    let mut space_len = space_cnr.clone();
    space_len.extend(space_reg);
    for &c in space_cnr.iter().rev() {
        let s = len_s - c;
        space_len.push(s);
    }

    // place points on 2D line
    let mut node_pos = Vec::new();
    for s in space_len {
        let x = x0 + s * len_x / len_s;
        let y = y0 + s * len_y / len_s;
        node_pos.push((x, y));
    }

    // return
    node_pos
}
