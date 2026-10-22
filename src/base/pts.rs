use crate::base::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;
use std::f64::consts::PI;

pub const INVALID: usize = usize::MAX;

#[derive(Default)]
pub struct Points {
    // sampling parameters
    // spacing grows linearly from h_min at the boundary to h_max at a distance ref_len
    pub h_min: f64,
    pub h_max: f64,
    pub ref_len: f64,
    pub max_try: usize,
    pub seed: u64,

    // points
    // points are located inside boundaries
    pub pts_pos: Vec<(f64, f64)>,  // gpt (global point id) -> (x, y)

    // 2d groups
    // uses the same grp2d_id in nodes
    pub grp2d_gpt: HashMap<String, Vec<usize>>,  // grp2d_id -> [gpt, ...]
}

impl Points {
    pub fn new(nds: &Nodes, h_min: f64, h_max: f64, ref_len: f64, max_try: usize, seed: u64) -> Result<Points, VFVChemError> {
        // input checking
        if max_try == 0 { return Err(VFVChemError::PointsNonPositiveAttempts { max_try }) }
        if !(h_min > 0.0 && h_min.is_finite()) { return Err(VFVChemError::PointsInvalidMinSpacing { h_min }) }
        if !(h_max >= h_min && h_max.is_finite()) { return Err(VFVChemError::PointsInvalidMaxSpacing { h_min, h_max }) }
        if !(ref_len > 0.0 && ref_len.is_finite()) { return Err(VFVChemError::PointsInvalidRefinementLength { ref_len }) }

        // default struct
        let mut pts = Points::default();
        pts.h_min = h_min;
        pts.h_max = h_max;
        pts.ref_len = ref_len;
        pts.max_try = max_try;
        pts.seed = seed;

        // iterate through 2d groups
        for (grp2d_id, loop_gnd) in nds.grp2d_gnd.iter() {
            // create internal points
            let pts_pos = _point_sample(nds, loop_gnd, h_min, h_max, ref_len, max_try, seed);
            let mut pts_id = Vec::new();

            // register internal points
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

fn _point_sample(nds: &Nodes, loop_gnd: &Vec<Vec<usize>>, h_min: f64, h_max: f64, ref_len: f64, max_try: usize, seed: u64) -> Vec<(f64, f64)> {
    // pseudocode:
    // for one 2d group, find the loops and bounding box
    // use loop winding to test if candidate point is inside
    // use distance to the loops to get the local spacing of each candidate
    // use bounding box to create background grid
    // seed every ccw loop, then use variable-radius bridson poisson-disk sampling to find points
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
    // cells are sized by the smallest spacing; each cell keeps a list since nodes near corners can share a cell
    let cell = h_min / f64::sqrt(2.0);
    let nx = (((x_max - x_min) / cell).ceil() as usize).max(1);
    let ny = (((y_max - y_min) / cell).ceil() as usize).max(1);
    let mut grid = vec![Vec::new(); nx * ny];  // (i + j*nx) -> [vid, ...]
   
    // initialize point positions, spacings, and active points
    let mut pts_pos = Vec::new();  // vid -> (x, y)
    let mut pts_h = Vec::new();  // vid -> local spacing
    let mut pts_act = Vec::new();  // point indices where neighbors are to be added

    // load node points into grid
    // nodes are on the loops, so their spacing is h_min
    for nds_pos in &loop_pos {
    for &(xc, yc) in nds_pos {
        // get grid cell
        let vid = pts_pos.len();
        let ic = (((xc - x_min) / cell) as usize).min(nx - 1);
        let jc = (((yc - y_min) / cell) as usize).min(ny - 1);
        
        // register point
        grid[ic + jc * nx].push(vid);
        pts_pos.push((xc, yc));
        pts_h.push(h_min);
    }}

    // number of node points
    let num_nds = pts_pos.len();

    // step 3: insert one starting point per ccw loop

    // seed rng
    let mut rng = StdRng::seed_from_u64(seed);

    // growth only reaches twice the local spacing from existing points, so separate regions need their own starting point
    // candidates are the loop centroid, then a distance h_min or h_max left of each edge midpoint
    // the h_max offset always clears the loop nodes, while the h_min offset fits in thin regions
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
            for h in [h_min, h_max] {
                cand.push((0.5 * (x0 + x1) - h * (y1 - y0) / len, 0.5 * (y0 + y1) + h * (x1 - x0) / len));
            }
        }

        // register the first valid candidate
        for (x0, y0) in cand {
            // check if point is (1) in the domain and (2) far enough from other points at its local spacing
            if !_point_in_poly(x0, y0, &loop_pos) { continue; }
            let h0 = _point_space(x0, y0, &loop_pos, h_min, h_max, ref_len);
            if !_point_is_far(x0, y0, h0, &pts_pos, &pts_h, &grid, x_min, y_min, nx, ny, cell, h_max) { continue; }

            // get grid cell
            let vid_c = pts_pos.len();
            let i0 = (((x0 - x_min) / cell) as usize).min(nx - 1);
            let j0 = (((y0 - y_min) / cell) as usize).min(ny - 1);

            // register initial point
            grid[i0 + j0 * nx].push(vid_c);
            pts_pos.push((x0, y0));
            pts_h.push(h0);
            pts_act.push(vid_c);
            break;
        }
    }

    // step 4: insert remaining points

    // insert points in the rest of the grid
    while !pts_act.is_empty() {
        // choose random active point
        let a = rng.random_range(0..pts_act.len());
        let vid_a = pts_act[a];
        let (xa, ya) = pts_pos[vid_a];
        let ha = pts_h[vid_a];

        // try to insert point next to it
        let mut is_insert = false;  // true if successfully inserted
        for _ in 0..max_try {
            // candidate point a distance [ha, 2ha] from active point
            let ang = rng.random_range(0.0..(2.0 * PI));
            let rad = rng.random_range(ha..(2.0 * ha));
            let xc = xa + rad * ang.cos();
            let yc = ya + rad * ang.sin();

            // check if point is (1) in the polygon and (2) far enough from other points at its local spacing
            if !_point_in_poly(xc, yc, &loop_pos) { continue; }
            let hc = _point_space(xc, yc, &loop_pos, h_min, h_max, ref_len);
            if !_point_is_far(xc, yc, hc, &pts_pos, &pts_h, &grid, x_min, y_min, nx, ny, cell, h_max) { continue; }

            // get grid cell
            let vid_c = pts_pos.len();
            let ic = (((xc - x_min) / cell) as usize).min(nx - 1);
            let jc = (((yc - y_min) / cell) as usize).min(ny - 1);
            
            // register candidate point
            grid[ic + jc * nx].push(vid_c);
            pts_pos.push((xc, yc));
            pts_h.push(hc);
            pts_act.push(vid_c);
            is_insert = true;
            break;
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

fn _point_space(xc: f64, yc: f64, loop_pos: &Vec<Vec<(f64, f64)>>, h_min: f64, h_max: f64, ref_len: f64) -> f64 {
    // local spacing h(d) = h_min + (h_max - h_min) * min(d / ref_len, 1)
    // d is the distance to the nearest loop segment, so holes and interfaces are also refined

    // uniform spacing needs no distance
    if h_max == h_min { return h_min; }

    // iterate through each segment of each loop
    let mut d2 = f64::INFINITY;
    for nds_pos in loop_pos {
        for nid0 in 0..nds_pos.len() {
            // get endpoints of segment
            let (x0, y0) = nds_pos[nid0];
            let (x1, y1) = nds_pos[(nid0 + 1) % nds_pos.len()];

            // project test point onto the segment and clamp to its endpoints
            let (ex, ey) = (x1 - x0, y1 - y0);
            let len2 = ex * ex + ey * ey;
            let t = if len2 > 0.0 { (((xc - x0) * ex + (yc - y0) * ey) / len2).clamp(0.0, 1.0) } else { 0.0 };
            let dx = x0 + t * ex - xc;
            let dy = y0 + t * ey - yc;
            d2 = d2.min(dx * dx + dy * dy);
        }
    }

    h_min + (h_max - h_min) * (d2.sqrt() / ref_len).min(1.0)
}

fn _point_is_far(xc: f64, yc: f64, hc: f64, pts_pos: &[(f64, f64)], pts_h: &[f64], grid: &[Vec<usize>], x_min: f64, y_min: f64, nx: usize, ny: usize, cell: f64, h_max: f64) -> bool {
    // two points are far enough apart if their distance is at least the mean of their local spacings
    // existing spacings are at most h_max, so only cells within (hc + h_max) / 2 of the candidate are checked

    // get upper and lower bounds of cells to check
    let ic = (((xc - x_min) / cell) as usize).min(nx - 1);
    let jc = (((yc - y_min) / cell) as usize).min(ny - 1);
    let num = (0.5 * (hc + h_max) / cell).ceil() as usize;
    let i_lb = ic.saturating_sub(num);  // (ic - num).max(0) without underflow
    let i_ub = (ic + num).min(nx - 1);
    let j_lb = jc.saturating_sub(num);
    let j_ub = (jc + num).min(ny - 1);

    // compute distance of existing points from candidate
    // reject candidate if too close to any existing point
    for j in j_lb..=j_ub {
    for i in i_lb..=i_ub {
        for &vid_n in &grid[i + j * nx] {
            let (xn, yn) = pts_pos[vid_n];
            let dx = xc - xn;
            let dy = yc - yn;
            let hs = 0.5 * (hc + pts_h[vid_n]);
            if dx * dx + dy * dy < hs * hs { return false; }
        }
    }}

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    const H_MIN: f64 = 0.02;
    const H_MAX: f64 = 0.08;
    const REF_LEN: f64 = 0.2;

    fn _square(seed: u64) -> (Nodes, Points) {
        let nds = Nodes::new_from_bounds(0.0, 0.0, 1.0, 1.0, H_MIN).unwrap();
        let pts = Points::new(&nds, H_MIN, H_MAX, REF_LEN, 30, seed).unwrap();
        (nds, pts)
    }

    #[test]
    fn space_is_linear_in_distance() {
        let loop_pos = vec![vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]];
        let h_mid = 0.5 * (H_MIN + H_MAX);
        assert!((_point_space(0.5, 0.0, &loop_pos, H_MIN, H_MAX, REF_LEN) - H_MIN).abs() < 1e-12);
        assert!((_point_space(0.5, 0.5 * REF_LEN, &loop_pos, H_MIN, H_MAX, REF_LEN) - h_mid).abs() < 1e-12);
        assert!((_point_space(0.5, REF_LEN, &loop_pos, H_MIN, H_MAX, REF_LEN) - H_MAX).abs() < 1e-12);
        assert!((_point_space(0.5, 0.5, &loop_pos, H_MIN, H_MAX, REF_LEN) - H_MAX).abs() < 1e-12);
        assert!((_point_space(0.5, 0.5, &loop_pos, H_MIN, H_MIN, REF_LEN) - H_MIN).abs() < 1e-12);
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let nds = Nodes::new_from_bounds(0.0, 0.0, 1.0, 1.0, H_MIN).unwrap();
        assert!(matches!(Points::new(&nds, H_MIN, H_MAX, REF_LEN, 0, 42), Err(VFVChemError::PointsNonPositiveAttempts { .. })));
        assert!(matches!(Points::new(&nds, 0.0, H_MAX, REF_LEN, 30, 42), Err(VFVChemError::PointsInvalidMinSpacing { .. })));
        assert!(matches!(Points::new(&nds, f64::NAN, H_MAX, REF_LEN, 30, 42), Err(VFVChemError::PointsInvalidMinSpacing { .. })));
        assert!(matches!(Points::new(&nds, H_MIN, 0.5 * H_MIN, REF_LEN, 30, 42), Err(VFVChemError::PointsInvalidMaxSpacing { .. })));
        assert!(matches!(Points::new(&nds, H_MIN, f64::INFINITY, REF_LEN, 30, 42), Err(VFVChemError::PointsInvalidMaxSpacing { .. })));
        assert!(matches!(Points::new(&nds, H_MIN, H_MAX, 0.0, 30, 42), Err(VFVChemError::PointsInvalidRefinementLength { .. })));
        assert!(matches!(Points::new(&nds, H_MIN, H_MAX, f64::NAN, 30, 42), Err(VFVChemError::PointsInvalidRefinementLength { .. })));
    }

    #[test]
    fn sampling_is_deterministic() {
        let (_, pts_a) = _square(7);
        let (_, pts_b) = _square(7);
        assert!(pts_a.num_pts() > 0);
        assert_eq!(pts_a.pts_pos, pts_b.pts_pos);
    }

    #[test]
    fn points_respect_variable_spacing() {
        let (nds, pts) = _square(42);
        let loop_pos: Vec<Vec<(f64, f64)>> = nds.grp2d_gnd["dom"].iter().map(|gnd| gnd.iter().map(|&g| nds.nds_pos[g]).collect()).collect();
        let pts_h: Vec<f64> = pts.pts_pos.iter().map(|&(x, y)| _point_space(x, y, &loop_pos, H_MIN, H_MAX, REF_LEN)).collect();
        for (p, &(xp, yp)) in pts.pts_pos.iter().enumerate() {
            assert!(_point_in_poly(xp, yp, &loop_pos));
            for (q, &(xq, yq)) in pts.pts_pos.iter().enumerate() {
                if q == p { continue; }
                let dist = ((xp - xq).powi(2) + (yp - yq).powi(2)).sqrt();
                assert!(dist >= 0.5 * (pts_h[p] + pts_h[q]) - 1e-12);
            }
            for &(xn, yn) in &nds.nds_pos {
                let dist = ((xp - xn).powi(2) + (yp - yn).powi(2)).sqrt();
                assert!(dist >= 0.5 * (pts_h[p] + H_MIN) - 1e-12);
            }
        }
    }

    #[test]
    fn points_are_denser_near_boundary() {
        // band within 0.1 of the boundary has area 0.36; center square farther than 0.3 has area 0.16
        let (_, pts) = _square(42);
        let dist = |&(x, y): &(f64, f64)| x.min(1.0 - x).min(y).min(1.0 - y);
        let num_near = pts.pts_pos.iter().filter(|p| dist(p) < 0.1).count() as f64;
        let num_far = pts.pts_pos.iter().filter(|p| dist(p) > 0.3).count() as f64;
        assert!(num_near / 0.36 > 3.0 * num_far / 0.16);
    }
}
