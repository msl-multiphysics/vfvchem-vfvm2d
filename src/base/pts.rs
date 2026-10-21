use crate::base::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;
use std::f64::consts::PI;

pub const INVALID: usize = usize::MAX;

#[derive(Default)]
pub struct Points {
    // sampling parameters
    pub max_try: usize,
    pub seed: u64,

    // points
    // points are located inside boundaries
    pub pts_pos: Vec<(f64, f64)>,  // gpt (global point id) -> (x, y)
    pub pts_vel: Vec<(f64, f64)>,  // gpt -> (u, v)

    // 2d groups
    // uses the same grp2d_id in nodes
    pub grp2d_gpt: HashMap<String, Vec<usize>>,  // grp2d_id -> [gpt, ...]
}

impl Points {
    pub fn new(nds: &Nodes, max_try: usize, seed: u64) -> Result<Points, VFVChemError> {
        // input checking
        if max_try == 0 { return Err(VFVChemError::PointsNonPositiveAttempts { max_try }) }

        // default struct
        let mut pts = Points::default();
        pts.max_try = max_try;
        pts.seed = seed;

        // iterate through 2d groups
        for (grp2d_id, loop_gnd) in nds.grp2d_gnd.iter() {
            // create internal points
            let pts_pos = _point_sample(nds, loop_gnd, max_try, seed);
            let mut pts_id = Vec::new();

            // register internal points
            for pos in pts_pos.iter() {
                let gpt = pts.pts_pos.len();
                pts.pts_pos.push(*pos);
                pts.pts_vel.push((0.0, 0.0));
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

fn _point_sample(nds: &Nodes, loop_gnd: &Vec<Vec<usize>>, max_try: usize, seed: u64) -> Vec<(f64, f64)> {
    // pseudocode:
    // for one 2d group, find the loops and bounding box
    // use loop winding to test if candidate point is inside
    // use bounding box to create background grid
    // seed every ccw loop, then use bridson poisson-disk sampling to find points
    // any number of ccw regions and cw holes are allowed

    // step 1: get polygon and bounding box

    // initialize
    let mut loop_pos = Vec::new();
    let mut x_min = f64::INFINITY;
    let mut y_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_max = f64::NEG_INFINITY;

    // iterate through nodes of each loop
    for nds_id in loop_gnd {
        let mut nds_pos = Vec::new();
        for &gnd in nds_id {
            let (x, y) = nds.nds_pos[gnd];
            nds_pos.push((x, y));
            x_min = x_min.min(x);
            y_min = y_min.min(y);
            x_max = x_max.max(x);
            y_max = y_max.max(y);
        }
        loop_pos.push(nds_pos);
    }
    
    // step 2: create grid
    
    // initialize grid
    // vid (vertex id) refers to the combined node and point ids
    let h = nds.space;
    let cell = h / f64::sqrt(2.0);
    let nx = ((x_max - x_min) / cell).ceil() as usize;
    let ny = ((y_max - y_min) / cell).ceil() as usize;
    let mut grid = vec![INVALID; nx * ny];  // vid of point in each cell; (i + j*nx) -> vid
   
    // initialize point positions and active points
    let mut pts_pos = Vec::new();  // vid -> (x, y)
    let mut pts_act = Vec::new();  // point indices where neighbors are to be added

    // load node points into grid
    for nds_pos in &loop_pos {
    for &(xc, yc) in nds_pos {
        // get grid cell
        let vid = pts_pos.len();
        let ic = (((xc - x_min) / cell) as usize).min(nx - 1);
        let jc = (((yc - y_min) / cell) as usize).min(ny - 1);
        
        // register point if no point was added yet
        if grid[ic + jc * nx] != INVALID { continue; }
        grid[ic + jc * nx] = vid;
        pts_pos.push((xc, yc));
    }}

    // number of grid cells with at least one node
    // corner nodes may be in the same cell
    let num_nds = pts_pos.len();

    // step 3: insert one starting point per ccw loop

    // seed rng
    let mut rng = StdRng::seed_from_u64(seed);

    // growth only reaches 2h from existing points, so separate regions need their own starting point
    // candidates are the loop centroid, then a distance h left of each edge midpoint
    let h2 = h * h;
    for nds_pos in &loop_pos {
        // skip degenerate loops
        let num = nds_pos.len();
        if num < 3 { continue; }

        // skip cw loops (holes) using the signed area
        let mut area = 0.0;
        for nid0 in 0..num {
            let (x0, y0) = nds_pos[nid0];
            let (x1, y1) = nds_pos[(nid0 + 1) % num];
            area += 0.5 * (x0 * y1 - x1 * y0);
        }
        if area <= 0.0 { continue; }

        // candidate starting points
        let xg = nds_pos.iter().map(|&(x, _)| x).sum::<f64>() / (num as f64);
        let yg = nds_pos.iter().map(|&(_, y)| y).sum::<f64>() / (num as f64);
        let mut cand = vec![(xg, yg)];
        for nid0 in 0..num {
            let (x0, y0) = nds_pos[nid0];
            let (x1, y1) = nds_pos[(nid0 + 1) % num];
            let len = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
            if len == 0.0 { continue; }
            cand.push((0.5 * (x0 + x1) - h * (y1 - y0) / len, 0.5 * (y0 + y1) + h * (x1 - x0) / len));
        }

        // register the first valid candidate
        for (x0, y0) in cand {
            // check if point is (1) in the domain and (2) far enough from other points
            let is_poly = _point_in_poly(x0, y0, &loop_pos);
            let is_far = _point_is_far(x0, y0, &pts_pos, &grid, x_min, y_min, nx, ny, h2, cell);

            // register initial point
            if is_poly && is_far {
                // get grid cell
                let vid_c = pts_pos.len();
                let i0 = ((x0 - x_min) / cell) as usize;
                let j0 = ((y0 - y_min) / cell) as usize;

                // register point
                grid[i0 + j0 * nx] = vid_c;
                pts_pos.push((x0, y0));
                pts_act.push(vid_c);
                break;
            }
        }
    }

    // step 4: insert remaining points

    // insert points in the rest of the grid
    while !pts_act.is_empty() {
        // choose random active point
        let a = rng.random_range(0..pts_act.len());
        let vid_a = pts_act[a];
        let (xa, ya) = pts_pos[vid_a];

        // try to insert point next to it
        let mut is_insert = false;  // true if successfully inserted
        for _ in 0..max_try {
            // candidate point a distance [h, 2h] from active point
            let ang = rng.random_range(0.0..(2.0 * PI));
            let rad = rng.random_range(h..(2.0 * h));
            let xc = xa + rad * ang.cos();
            let yc = ya + rad * ang.sin();

            // check if point is (1) in the polygon and (2) far enough from other points
            let is_poly = _point_in_poly(xc, yc, &loop_pos);
            let is_far = _point_is_far(xc, yc, &pts_pos, &grid, x_min, y_min, nx, ny, h2, cell);
            
            // register candidate point 
            if is_poly && is_far {
                // get grid cell
                let vid_c = pts_pos.len();
                let ic = ((xc - x_min) / cell) as usize;
                let jc = ((yc - y_min) / cell) as usize;
                
                // register point
                grid[ic + jc * nx] = vid_c;
                pts_pos.push((xc, yc));
                pts_act.push(vid_c);
                is_insert = true;
                break;
            }
        }

        // flag active point as inactive if candidates cannot be inserted
        if !is_insert {
            pts_act.swap_remove(a);
        }
    }

    // return internal points
    // skip first num_nds points
    pts_pos.split_off(num_nds)

}

fn _point_in_poly(xc: f64, yc: f64, loop_pos: &Vec<Vec<(f64, f64)>>) -> bool {
    // winding number over each closed loop
    // ccw loops add +1 and cw loops add -1, so the point is inside if the winding is positive
    
    // initialize
    let mut wind = 0;
    
    // iterate through loops
    for nds_pos in loop_pos {
        // skip degenerate loops
        if nds_pos.len() < 3 { continue; }

        // iterate through each segment in traversal order
        for nid0 in 0..nds_pos.len() {
            // get endpoints of segment
            let (x0, y0) = nds_pos[nid0];
            let (x1, y1) = nds_pos[(nid0 + 1) % nds_pos.len()];

            // side of test point; positive if on the left of the segment
            let side = (x1 - x0) * (yc - y0) - (xc - x0) * (y1 - y0);

            // upward crossing with point on the left, or downward crossing with point on the right
            if y0 <= yc && y1 > yc && side > 0.0 { wind += 1; }
            if y0 > yc && y1 <= yc && side < 0.0 { wind -= 1; }
        }
    }

    wind > 0
}

fn _point_is_far(xc: f64, yc: f64, pts_pos: &[(f64, f64)], grid: &[usize], x_min: f64, y_min: f64, nx: usize, ny: usize, h2: f64, cell: f64) -> bool{
    // get upper and lower bounds of cells to check
    // get points +/- 2 cells away from candidate
    let ic = ((xc - x_min) / cell) as usize;
    let jc = ((yc - y_min) / cell) as usize;
    let i_lb = ic.saturating_sub(2); // (ic - 2).max(0) without underflow
    let i_ub = (ic + 2).min(nx - 1);
    let j_lb = jc.saturating_sub(2);
    let j_ub = (jc + 2).min(ny - 1);

    // compute distance of existing points from candidate
    // reject candidate if too close to any existing point
    let mut is_far = true;
    for j in j_lb..=j_ub {
    for i in i_lb..=i_ub {
        let pid_n = grid[i + j * nx];
        if pid_n == INVALID { continue; }  // empty cell
        let (xn, yn) = pts_pos[pid_n];
        let dx = xc - xn;
        let dy = yc - yn;
        let ds2 = dx * dx + dy * dy;
        if ds2 < h2 {
            is_far = false;
            break;
        }
    }}

    is_far
}
