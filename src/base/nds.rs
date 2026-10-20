use crate::base::prelude::*;
use std::collections::HashMap;

const OFFSET_RATIO: f64 = 0.1;

#[derive(Default)]
pub struct Nodes {
    // spacing
    pub space: f64,

    // nodes
    // nodes are on boundaries and interfaces
    pub nds_pos: Vec<(f64, f64)>,  // gnd (global node id) -> (x, y)
    pub nds_vel: Vec<(f64, f64)>,  // gnd -> (u, v)
    
    // 1d groups
    // defines boundaries and interfaces
    pub grp1d_gnd: HashMap<String, Vec<usize>>, // grp1d_id -> [gnd, ...]

    // 2d groups
    // defines domain and outward normals
    // outer vector defines loop; inner vector defines gnd and normals
    pub grp2d_gnd: HashMap<String, Vec<Vec<usize>>>, // grp2d_id, loop -> [gnd, ...]
    pub grp2d_norm: HashMap<String, Vec<Vec<(f64, f64)>>>, // grp2d_id, loop -> [(nx, ny), ...]
}

impl Nodes {
    pub fn new_from_obj(file_path: &str, h: f64) -> Result<Nodes, VFVChemError> {
        let nds = read_nds_obj(file_path, h, OFFSET_RATIO)?;
        Ok(nds)
    }

    pub fn new_from_bounds(x_min: f64, y_min: f64, x_max: f64, y_max: f64, h: f64) -> Result<Nodes, VFVChemError> {
        // creates boundary nodes around the bounding box
        // nodes are spaced by at least h
        // nodes are offset from corners by a distance dh
        // no corner nodes to prevent incompatible conditions and simplify bookkeeping

        // input checking
        if x_max < x_min { return Err(VFVChemError::NodesInvalidBoundsX { x_min, x_max }) }
        if y_max < y_min { return Err(VFVChemError::NodesInvalidBoundsY { y_min, y_max }) }
        if h <= 0.0 { return Err(VFVChemError::NodesNonPositiveSpacing { h }) }

        // default struct
        let mut nds = Nodes::default();
        nds.space = h;

        // get spacing along x and y; must be greater than h
        let num_x = ((x_max - x_min + 1e-6) / h) as usize;  // number of segments
        let num_y = ((y_max - y_min + 1e-6) / h) as usize;
        let dx = (x_max - x_min) / (num_x as f64);
        let dy = (y_max - y_min) / (num_y as f64);

        // construct nodes
        let dh = OFFSET_RATIO * h;
        let mut nds_b = vec![(x_min + dh, y_min)];  // bottom
        for i in 1..num_x {  
            let x = x_min + (i as f64) * dx;
            nds_b.push((x, y_min));
        }
        nds_b.push((x_max - dh, y_min));
        let mut nds_r = vec![(x_max, y_min + dh)];  // right
        for j in 1..num_y {
            let y = y_min + (j as f64) * dy;
            nds_r.push((x_max, y));
        }
        nds_r.push((x_max, y_max - dh));
        let mut nds_t = vec![(x_max - dh, y_max)];  // top
        for i in 1..num_x {  // top
            let x = x_max - (i as f64) * dx;
            nds_t.push((x, y_max));
        }
        nds_t.push((x_min + dh, y_max));
        let mut nds_l = vec![(x_min, y_max - dh)];  // left
        for j in 1..num_y {
            let y = y_max - (j as f64) * dy;
            nds_l.push((x_min, y));
        }
        nds_l.push((x_min, y_min + dh));

        // store nodes
        nds.nds_pos.extend(nds_b);
        nds.nds_pos.extend(nds_r);
        nds.nds_pos.extend(nds_t);
        nds.nds_pos.extend(nds_l);
        nds.nds_vel = vec![(0.0, 0.0); 2*num_x + 2*num_y + 4];

        // nodes in 1d groups
        nds.grp1d_gnd.insert("bnd_b".to_string(), (0..num_x + 1).collect());                                // bottom
        nds.grp1d_gnd.insert("bnd_r".to_string(), (num_x + 1..num_x + num_y + 2).collect());                // right
        nds.grp1d_gnd.insert("bnd_t".to_string(), (num_x + num_y + 2..2*num_x + num_y + 3).collect());      // top
        nds.grp1d_gnd.insert("bnd_l".to_string(), (2*num_x + num_y + 3..2*num_x + 2*num_y + 4).collect());  // left

        // nodes in 2d groups
        nds.grp2d_gnd.insert("dom".to_string(), vec![(0..2*num_x + 2*num_y + 4).collect()]);

        // outward normals in 2d groups
        let mut nds_norm = Vec::new();
        nds_norm.extend(vec![( 0.0, -1.0); num_x + 1]);  // bottom
        nds_norm.extend(vec![( 1.0,  0.0); num_y + 1]);  // right
        nds_norm.extend(vec![( 0.0,  1.0); num_x + 1]);  // top
        nds_norm.extend(vec![(-1.0,  0.0); num_y + 1]);  // left
        nds.grp2d_norm.insert("dom".to_string(), vec![nds_norm]);

        Ok(nds)
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_pos.len()
    }
}
