use crate::base::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;
use std::f64::consts::PI;

const BOUNDARY_RATIO: f64 = 1.2;
const MAX_TRY: usize = 100;
const SEED: u64 = 317;
const DELAUNAY_RATIO: f64 = 0.6;  // must be greater than 0.5 so that boundaries become delaunay edges

#[derive(Default)]
pub struct Points {
    // points
    // points are located inside boundaries
    pub pts_pos: Vec<(f64, f64)>,  // gpt (global point id) -> (x, y)

    // 2d groups
    // uses the same grp2d_id in nodes
    pub grp2d_gpt: HashMap<String, Vec<usize>>,  // grp2d_id -> [gpt, ...]

    // spacing
    pub pts_space: f64,  // point spacing
    pub bnd_ratio: f64,  // geometric growth rate from node spacing to point spacing
}

impl Points {
    pub fn new(nds: &Nodes, pts_space: f64) -> Result<Points, VFVChemError> {
        Points::new_options(nds, pts_space, BOUNDARY_RATIO, MAX_TRY, SEED)
    }

    pub fn new_options(nds: &Nodes, pts_space: f64, bnd_ratio: f64, max_try: usize, seed: u64) -> Result<Points, VFVChemError> {
        // input checking
        if max_try == 0 { return Err(VFVChemError::PointsNonPositiveAttempts { max_try }) }
        if !(pts_space >= nds.nds_space) { return Err(VFVChemError::PointsInvalidSpacing { pts_space }) }
        if bnd_ratio <= 1.0 { return Err(VFVChemError::PointsInvalidBoundaryRatio { bnd_ratio }) }

        // default struct
        let mut pts = Points::default();
        pts.pts_space = pts_space;  // point spacing
        pts.bnd_ratio = bnd_ratio;  // geometric growth rate from node spacing to point spacing

        // iterate through 2d groups
        for (grp2d_id, curve_eid) in nds.grp2d_eid.iter() {
            // create points
            let pts_pos = _point_sample(nds, curve_eid, pts_space, bnd_ratio, max_try, seed)?;
            
            // register points
            let mut pts_id = Vec::new();
            for pos in pts_pos.iter() {
                let gpt = pts.pts_pos.len();
                pts.pts_pos.push(*pos);
                pts_id.push(gpt);
            }
            pts.grp2d_gpt.insert(grp2d_id.clone(), pts_id);
        }

        Ok(pts)
    }

    #[inline]
    pub fn num_pts(&self) -> usize {
        self.pts_pos.len()
    }
}

fn _point_sample(nds: &Nodes, curve_eid: &Vec<usize>, pts_space: f64, bnd_ratio: f64, max_try: usize, seed: u64) -> Result<Vec<(f64, f64)>, VFVChemError> {
    // outline:
    // create bounding box and grid
    // get initial sampling point
    // get remaining sampling points
    //
    // notes:
    // point is feasible if far enough from boundaries and inside region

    // step 1: create bounding box

    // initialize bounding box
    let mut x_min = f64::INFINITY;
    let mut y_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_max = f64::NEG_INFINITY;

    // initialize vertex positions per curve
    // not used for bounding box but used later on
    let mut curve_pos = Vec::new();

    // iterate through curves
    for &eid in curve_eid.iter() {
        let mut edge_pos = [(0.0, 0.0); 2];
        for (i, &gvt) in nds.edge_gvt[eid].iter().enumerate() {
            let (x, y) = nds.vert_pos[gvt];
            x_min = x_min.min(x);
            y_min = y_min.min(y);
            x_max = x_max.max(x);
            y_max = y_max.max(y);
            edge_pos[i] = (x, y);
        }
        curve_pos.push(edge_pos);
    }
    
    // step 2: create grid

    // initialize point list
    let mut pts_pos = Vec::new();  // pid -> (x, y)
    let mut pts_h = Vec::new();  // pid -> local spacing

    // compute grid spacing
    // set grid spacing to regular node spacing
    let dh = nds.nds_space;

    // create grid to point id mapping
    // grid ids (gid) is i + j * num_x
    // point ids (pid) here are not the same as those used outside this function
    let num_x = ((x_max - x_min) / dh).ceil() as usize;
    let num_y = ((y_max - y_min) / dh).ceil() as usize;
    let mut grid = vec![vec![]; num_x * num_y];  // gid -> [pid, ...]

    // get node spacing parameters
    let nds_space = nds.nds_space;
    let off_space = nds.off_space;
    let crn_space = nds.crn_space;
    let crn_ratio = nds.crn_ratio;

    // load nodes into grid
    for &eid in curve_eid.iter() {
    for &gnd in nds.edge_gnd[eid].iter() {
        // get point and grid ids
        let (x, y) = nds.nds_pos[gnd];
        let i = (((x - x_min) / dh) as usize).min(num_x - 1);
        let j = (((y - y_min) / dh) as usize).min(num_y - 1);
        let pid = pts_pos.len();
        let gid = i + j * num_x;
        
        // get local spacing (includes corner refinement)
        let (h, _) = _point_space(&curve_pos, nds_space, off_space, crn_space, crn_ratio, pts_space, bnd_ratio, x, y);

        // store point
        grid[gid].push(pid);
        pts_pos.push((x, y));
        pts_h.push(h);
    }}

    // number of grid cells occupied by nodes
    let num_nds = pts_pos.len();

    // step 3: get initial sampling point

    // seed rng
    let mut rng = StdRng::seed_from_u64(seed);

    // initialize active points
    let mut pts_act: Vec<usize> = (0..num_nds).collect();

    // insert initial point
    for _ in 0..max_try {
        // candidate position
        let xc = rng.random_range(x_min..x_max);
        let yc = rng.random_range(y_min..y_max);

        // test point
        // rejects point if not feasible
        if !_point_inside(&curve_pos, xc, yc) { continue; };
        let (hc, dc) = _point_space(&curve_pos, nds_space, off_space, crn_space, crn_ratio, pts_space, bnd_ratio, xc, yc);
        if dc < DELAUNAY_RATIO * hc { continue; };  // reject if spacing is greater than distance to boundary
        if !_point_far(&pts_pos, &grid, x_min, y_min, num_x, num_y, dh, hc, xc, yc) { continue; };

        // register point if feasible
        // get point and grid ids
        let ic = (((xc - x_min) / dh) as usize).min(num_x - 1);
        let jc = (((yc - y_min) / dh) as usize).min(num_y - 1);
        let pid = pts_pos.len();
        let gid = ic + jc * num_x;

        // store point
        grid[gid].push(pid);
        pts_pos.push((xc, yc));
        pts_h.push(hc);
        pts_act.push(pid);  // set initial point as active

        // stop search for initial point
        break;
    }

    // step 4: get remaining sampling points

    // continue until no more active points
    while !pts_act.is_empty() {
        // choose active point
        let idx_a = rng.random_range(0..pts_act.len());
        let pid_a = pts_act[idx_a];
        let (xa, ya) = pts_pos[pid_a];
        let ha = pts_h[pid_a];

        // insert neighbor point
        let mut is_insert = false;  // true if inserted
        for _ in 0..max_try {
            // candidate point
            let ang = rng.random_range(0.0..(2.0 * PI));
            let rad = rng.random_range(0.0..(2.0 * ha));
            let xc = xa + rad * ang.cos();
            let yc = ya + rad * ang.sin();

            // test point
            // rejects point if not feasible
            if !_point_inside(&curve_pos, xc, yc) { continue; };
            let (hc, dc) = _point_space(&curve_pos, nds_space, off_space, crn_space, crn_ratio, pts_space, bnd_ratio, xc, yc);
            if dc < DELAUNAY_RATIO * hc { continue; };
            if !_point_far(&pts_pos, &grid, x_min, y_min, num_x, num_y, dh, hc, xc, yc) { continue; };

            // register point if feasible
            // get point and grid ids
            let ic = (((xc - x_min) / dh) as usize).min(num_x - 1);
            let jc = (((yc - y_min) / dh) as usize).min(num_y - 1);
            let pid = pts_pos.len();
            let gid = ic + jc * num_x;

            // store point
            grid[gid].push(pid);
            pts_pos.push((xc, yc));
            pts_h.push(hc);
            pts_act.push(pid);  // set new point as active

            // stop search for neighbor point
            is_insert = true;
            break;
        }

        // flag as inactive if neighbors cannot be inserted
        if !is_insert {
            pts_act.swap_remove(idx_a);
        }

    }
    
    // skip nodes (first num_nds points) to get internal points only
    pts_pos.drain(0..num_nds);  

    Ok(pts_pos)
}

fn _point_inside(curve_pos: &Vec<[(f64, f64); 2]>, xc: f64, yc: f64) -> bool {
    // winding number algorithm
    // +1 for ccw curves and -1 for cw curves
    // point is inside if winding number is positive

    // initialize winding number
    let mut wind = 0;

    // iterate through curves
    for &[(x0, y0), (x1, y1)] in curve_pos {
        // cross product of (xc - x0) and (x1 - x0)
        // positive if xc on left of segment
        let cross = (x1 - x0) * (yc - y0) - (xc - x0) * (y1 - y0); 

        // CCW - point on left; segment goes up
        // CW - point on right; segment goes down
        if y0 <= yc && y1 > yc && cross > 0.0 { wind += 1; }
        if y0 > yc && y1 <= yc && cross < 0.0 { wind -= 1; }
    }

    // return
    wind > 0
}

fn _point_space(curve_pos: &Vec<[(f64, f64); 2]>, nds_space: f64, off_space: f64, crn_space: f64, crn_ratio: f64, pts_space: f64, bnd_ratio: f64, xc: f64, yc: f64) -> (f64, f64) {
    // outline
    // find shortest perpendicular distance to segments
    // find shortest distance to corners if within corner distance
    // if corner distance is short enough, use corner refinement length
    // otherwise, use boundary refinement length
    // return spacing and distance from segments

    // initialize minimum distances
    let mut dist2_bnd = f64::INFINITY;  // from boundary
    let mut dist2_crn = f64::INFINITY;  // from corner

    // find shortest perpendicular distance to segments
    // iterate through curves
    for &[(x0, y0), (x1, y1)] in curve_pos {
        // distance from segment
        
        // project (xc - x0) onto (x1 - x0)
        // clamp to segment endpoints
        let len_x = x1 - x0;
        let len_y = y1 - y0;
        let len_s2 = len_x * len_x + len_y * len_y;
        let t = (((xc - x0) * len_x + (yc - y0) * len_y) / len_s2).clamp(0.0, 1.0);
        let xp = x0 + t * len_x;
        let yp = y0 + t * len_y;

        // get distance to projection point
        let dist_x = xc - xp;
        let dist_y = yc - yp;
        let dist2 = dist_x * dist_x + dist_y * dist_y;
        dist2_bnd = dist2_bnd.min(dist2);

        // distance from corner
        
        // compute distance to point 0
        let dist_x = xc - x0;
        let dist_y = yc - y0;
        let dist2 = dist_x * dist_x + dist_y * dist_y;
        dist2_crn = dist2_crn.min(dist2);

        // compute distance to point 1
        let dist_x = xc - x1;
        let dist_y = yc - y1;
        let dist2 = dist_x * dist_x + dist_y * dist_y;
        dist2_crn = dist2_crn.min(dist2);

    }

    // compute spacing using corner refinement
    let a = off_space * nds_space;
    let b = crn_space * nds_space;
    let n = ((1.0 + (b - a) * (crn_ratio - 1.0) / a).ln() / crn_ratio.ln()) as i32;
    let hc_crn = crn_ratio.ln() * ((b - a) / (crn_ratio.powi(n) - 1.0) + dist2_crn.sqrt() - a);
    
    // compute spacing using boundary refinement
    let p = nds_space;
    let q = pts_space;
    let m = ((1.0 + (q - p) * (bnd_ratio - 1.0) / p).ln() / bnd_ratio.ln()) as i32;
    let dist_bnd = dist2_bnd.sqrt();
    let hc_bnd = bnd_ratio.ln() * ((q - p) / (bnd_ratio.powi(m) - 1.0) + dist_bnd - p);
    let hc_bnd = hc_bnd.min(q);
    
    // return
    (hc_crn.min(hc_bnd), dist_bnd)
}

fn _point_far(pts_pos: &Vec<(f64, f64)>, grid: &Vec<Vec<usize>>, x_min: f64, y_min: f64, num_x: usize, num_y: usize, dh: f64, hc: f64, xc: f64, yc: f64) -> bool {
    // search only within hc of (xc, yc)

    // get bounds of cells to check
    let ic = (((xc - x_min) / dh) as usize).min(num_x - 1);
    let jc = (((yc - y_min) / dh) as usize).min(num_y - 1);
    let num_cell = (hc / dh).ceil() as usize;
    let i_lb = ic.saturating_sub(num_cell);  // (ic - num_cell).max(0) without underflow
    let i_ub = (ic + num_cell).min(num_x - 1);
    let j_lb = jc.saturating_sub(num_cell);
    let j_ub = (jc + num_cell).min(num_y - 1);

    // compute distance from candidate to neighbors
    let h2 = hc * hc;
    for i in i_lb..=i_ub {
    for j in j_lb..=j_ub {
        // skip if cell is empty
        let gid = i + j * num_x;
        let pid_sub = &grid[gid];
        if pid_sub.is_empty() { continue; }

        // compute distance to neighbors
        // return if neighbor is too close
        for &pid in pid_sub.iter() {
            let (xn, yn) = pts_pos[pid];
            let dist_x = xn - xc;
            let dist_y = yn - yc;
            let dist2 = dist_x * dist_x + dist_y * dist_y;        
            if dist2 < h2 { return false; }
        }
    }}

    // return
    true
}
