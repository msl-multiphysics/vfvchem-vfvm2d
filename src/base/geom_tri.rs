use crate::base::prelude::*;
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct Triangle {
    // ids
    pub tri_id: usize,
    pub dom_id: usize,  // domain this tri refers to

    // triangulation data
    // vid is the index of combined point and node vectors
    // i.e., vid = if vid < dom.num_pts() {pid} else {dom.num_pts() + nid}
    pub tri_vid: Vec<[usize; 3]>,  // [vid, vid, vid] counterclockwise

    // control volume data
    pub cvol_area: Vec<f64>,  // vid (vertex id) -> CV (control volume) area
    pub cvol_skew: Vec<(f64, f64)>,  // vid -> CV centroid minus vertex position

    // face data
    pub face_len: Vec<Vec<f64>>,  // vid, loc (local neighbor index) -> signed face lengths; negative if circumcenters cross the edge
    pub face_norm: Vec<Vec<(f64, f64)>>,  // vid, loc -> outward normal on face
    pub face_skew: Vec<Vec<(f64, f64)>>,  // vid, loc -> face center minus edge midpoint

    // neighbor vertex data
    pub ngh_vid: Vec<Vec<usize>>,  // vid -> neighbor vid
    pub ngh_dist: Vec<Vec<f64>>,  // vid -> edge length to neighbor

    // boundary vertex data
    // boundary vertices own half of the boundary edges on either side
    pub bnd_len: Vec<f64>,  // vid -> CV length along boundary; zero for points
    pub bnd_vid: Vec<Vec<usize>>,  // vid, blc (boundary loc) -> neighbor vid along boundary; empty for points
    pub bnd_half: Vec<Vec<f64>>,  // vid, blc -> half-edge length along boundary
    pub bnd_norm: Vec<Vec<(f64, f64)>>,  // vid, blc -> outward normal on half-edge

    // gradient weight data
    // delta_u = [u1 - u0, u2 - u0, ...]; u0 and uj are point and neighbor values
    // du/dx = wgt_x . delta_u; du/dy = wgt_y . delta_u
    pub grad_wgt: Vec<Vec<(f64, f64)>>,  // vid, loc -> least-squares weights (wgt_x, wgt_y)
}

impl Triangle {
    pub fn new(tri_id: usize, dom: &Domain) -> Result<Triangle, VFVChemError> {
        // default struct
        let mut tri = Triangle::default();
        tri.tri_id = tri_id;
        tri.dom_id = dom.dom_id;

        Ok(tri)
    }

    pub fn update_tri(&mut self, nds: &Nodes, pts: &Points, dom: &Domain) -> Result<(), VFVChemError> {
        // get loops
        let loops = &nds.grp2d_gnd[&dom.grp2d_id];

        // clear tri data
        self.tri_vid.clear();
        self.cvol_area = vec![0.0; dom.num_tot()];
        self.cvol_skew = vec![(0.0, 0.0); dom.num_tot()];
        self.face_len = vec![Vec::new(); dom.num_tot()];
        self.face_norm = vec![Vec::new(); dom.num_tot()];
        self.face_skew = vec![Vec::new(); dom.num_tot()];
        self.ngh_vid = vec![Vec::new(); dom.num_tot()];
        self.ngh_dist = vec![Vec::new(); dom.num_tot()];
        self.bnd_len = vec![0.0; dom.num_tot()];
        self.bnd_vid = vec![Vec::new(); dom.num_tot()];
        self.bnd_half = vec![Vec::new(); dom.num_tot()];
        self.bnd_norm = vec![Vec::new(); dom.num_tot()];
        self.grad_wgt = vec![Vec::new(); dom.num_tot()];

        // collect points in domain
        let mut pts_pos = Vec::with_capacity(dom.num_tot());  // vid -> (x, y)
        for &gpt in &dom.pts_gpt {
            pts_pos.push(pts.pts_pos[gpt]);
        }
        for &gnd in &dom.nds_gnd {
            pts_pos.push(nds.nds_pos[gnd]);
        }

        // step 1: constrained Delaunay triangulation

        // insert points
        let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::new();
        let mut hdl_vid = vec![INVALID; dom.num_tot()];  // vertex handle index -> vid
        let mut vid_hdl = Vec::with_capacity(dom.num_tot());  // vid -> vertex handle
        for (vid, &(x, y)) in pts_pos.iter().enumerate() {
            let hdl = cdt.insert(Point2::new(x, y)).map_err(|source| VFVChemError::TriangleInsertionError { source })?;
            hdl_vid[hdl.index()] = vid;
            vid_hdl.push(hdl);
        }

        // constrain boundary edges and get boundary lengths
        // each node gets half of the edges on either side of it
        let mut edge_loop = HashSet::new();  // directed boundary edges (vid_0, vid_1)
        for loop_gnd in loops {
            for (idx_0, &gnd_0) in loop_gnd.iter().enumerate() {
                let gnd_1 = loop_gnd[(idx_0 + 1) % loop_gnd.len()];
                let vid_0 = dom.num_pts() + dom.nds_nid[&gnd_0];
                let vid_1 = dom.num_pts() + dom.nds_nid[&gnd_1];

                // add constraint
                if !cdt.can_add_constraint(vid_hdl[vid_0], vid_hdl[vid_1]) { return Err(VFVChemError::TriangleConstraintIntersection { dom_id: dom.dom_id }) }
                cdt.add_constraint(vid_hdl[vid_0], vid_hdl[vid_1]);
                edge_loop.insert((vid_0, vid_1));

                // add boundary half-edges
                // loops are oriented so that the right of the traversal is out of the domain
                let (x0, y0) = pts_pos[vid_0];
                let (x1, y1) = pts_pos[vid_1];
                let len = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
                let norm = ((y1 - y0) / len, -(x1 - x0) / len);
                self.bnd_len[vid_0] += 0.5 * len;
                self.bnd_len[vid_1] += 0.5 * len;
                self.bnd_vid[vid_0].push(vid_1);
                self.bnd_vid[vid_1].push(vid_0);
                self.bnd_half[vid_0].push(0.5 * len);
                self.bnd_half[vid_1].push(0.5 * len);
                self.bnd_norm[vid_0].push(norm);
                self.bnd_norm[vid_1].push(norm);
            }
        }

        // keep triangles inside the domain
        // removes triangles outside the boundary and inside holes
        // loops keep the domain on the left, so a counterclockwise triangle containing a loop edge is inside
        // triangles on the convex hull are outside, and labels spread across unconstrained edges
        // no geometric test is used since sliver triangles between nearly collinear boundary nodes have centroids on the boundary
        let tri_all: Vec<[usize; 3]> = cdt.inner_faces().map(|face| {
            let [hdl_0, hdl_1, hdl_2] = face.vertices();
            [hdl_vid[hdl_0.fix().index()], hdl_vid[hdl_1.fix().index()], hdl_vid[hdl_2.fix().index()]]
        }).collect();
        let mut edge_tri = HashMap::new();  // directed edge (vid_i, vid_j) -> triangle index
        for (idx, tri) in tri_all.iter().enumerate() {
            for loc in 0..3 { edge_tri.insert((tri[loc], tri[(loc + 1) % 3]), idx); }
        }

        // label triangles on the boundary and convex hull
        let mut is_inside = vec![None; tri_all.len()];
        let mut stack = Vec::new();
        for (idx, tri) in tri_all.iter().enumerate() {
            for loc in 0..3 {
                let (vid_i, vid_j) = (tri[loc], tri[(loc + 1) % 3]);
                if edge_loop.contains(&(vid_i, vid_j)) { is_inside[idx] = Some(true); }
                else if edge_loop.contains(&(vid_j, vid_i)) || !edge_tri.contains_key(&(vid_j, vid_i)) { is_inside[idx] = Some(false); }
            }
            if is_inside[idx].is_some() { stack.push(idx); }
        }

        // spread labels across unconstrained edges
        while let Some(idx) = stack.pop() {
            let tri = tri_all[idx];
            for loc in 0..3 {
                let (vid_i, vid_j) = (tri[loc], tri[(loc + 1) % 3]);
                if edge_loop.contains(&(vid_i, vid_j)) || edge_loop.contains(&(vid_j, vid_i)) { continue; }
                if let Some(&idx_ngh) = edge_tri.get(&(vid_j, vid_i)) {
                    if is_inside[idx_ngh].is_none() {
                        is_inside[idx_ngh] = is_inside[idx];
                        stack.push(idx_ngh);
                    }
                }
            }
        }
        self.tri_vid = tri_all.into_iter().zip(is_inside).filter(|(_, is_in)| *is_in == Some(true)).map(|(tri, _)| tri).collect();

        // step 2: Voronoi faces

        // iterate through edges of each triangle
        // the face piece of edge (i, j) runs from the edge midpoint m to the circumcenter c
        // its signed length is 0.5 * cot(theta_k) * d_ij, where k is the opposite vertex
        // the signed sub-triangles (x_i, m, c) and (x_j, m, c) tile the triangle
        let mut edge_face = HashMap::new();  // (vid_min, vid_max) -> (signed face length, sum of length * (piece center - m))
        let mut cent_sum = vec![(0.0, 0.0); dom.num_tot()];  // vid -> sum of area * (sub-triangle centroid - node)
        for tri in &self.tri_vid {
            // get circumcenter relative to the first vertex
            let (x0, y0) = pts_pos[tri[0]];
            let (bx, by) = (pts_pos[tri[1]].0 - x0, pts_pos[tri[1]].1 - y0);
            let (cx, cy) = (pts_pos[tri[2]].0 - x0, pts_pos[tri[2]].1 - y0);
            let den = 2.0 * (bx * cy - by * cx);
            let xc = x0 + (cy * (bx * bx + by * by) - by * (cx * cx + cy * cy)) / den;
            let yc = y0 + (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by)) / den;

            for loc in 0..3 {
                // get edge (i, j) and opposite vertex k
                let vid_i = tri[loc];
                let vid_j = tri[(loc + 1) % 3];
                let vid_k = tri[(loc + 2) % 3];
                let (xi, yi) = pts_pos[vid_i];
                let (xj, yj) = pts_pos[vid_j];
                let (xk, yk) = pts_pos[vid_k];

                // cot(theta_k) = (e_ki . e_kj) / |e_ki x e_kj|
                let dot = (xi - xk) * (xj - xk) + (yi - yk) * (yj - yk);
                let crs = ((xi - xk) * (yj - yk) - (yi - yk) * (xj - xk)).abs();
                let dist = ((xj - xi).powi(2) + (yj - yi).powi(2)).sqrt();
                let face = 0.5 * dot / crs * dist;

                // accumulate face pieces; the piece center is (m + c) / 2
                let (xm, ym) = (0.5 * (xi + xj), 0.5 * (yi + yj));
                let entry = edge_face.entry((vid_i.min(vid_j), vid_i.max(vid_j))).or_insert((0.0, 0.0, 0.0));
                entry.0 += face;
                entry.1 += face * 0.5 * (xc - xm);
                entry.2 += face * 0.5 * (yc - ym);

                // accumulate sub-triangle centroids; each has signed area 0.25 * dist * face
                let area = 0.25 * dist * face;
                cent_sum[vid_i].0 += area * ((xi + xm + xc) / 3.0 - xi);
                cent_sum[vid_i].1 += area * ((yi + ym + yc) / 3.0 - yi);
                cent_sum[vid_j].0 += area * ((xj + xm + xc) / 3.0 - xj);
                cent_sum[vid_j].1 += area * ((yj + ym + yc) / 3.0 - yj);
            }
        }

        // step 3: neighbor, face, and control volume data
        for (&(vid_a, vid_b), &(face, skew_x_sum, skew_y_sum)) in edge_face.iter() {
            // get edge geometry
            let (xa, ya) = pts_pos[vid_a];
            let (xb, yb) = pts_pos[vid_b];
            let dist = ((xb - xa).powi(2) + (yb - ya).powi(2)).sqrt();
            let (nx, ny) = ((xb - xa) / dist, (yb - ya) / dist);

            // get face skewness; zero-length faces carry no flux
            let skew = if face.abs() > 1e-12 * dist { (skew_x_sum / face, skew_y_sum / face) } else { (0.0, 0.0) };

            // store for both directions; normals point out of each control volume
            self.ngh_vid[vid_a].push(vid_b);
            self.ngh_dist[vid_a].push(dist);
            self.face_len[vid_a].push(face);
            self.face_norm[vid_a].push((nx, ny));
            self.face_skew[vid_a].push(skew);
            self.ngh_vid[vid_b].push(vid_a);
            self.ngh_dist[vid_b].push(dist);
            self.face_len[vid_b].push(face);
            self.face_norm[vid_b].push((-nx, -ny));
            self.face_skew[vid_b].push(skew);

            // control volume area: each side of the face covers a triangle of area 0.5 * face * (0.5 * dist)
            self.cvol_area[vid_a] += 0.25 * dist * face;
            self.cvol_area[vid_b] += 0.25 * dist * face;
        }

        // control volume skewness
        for vid in 0..dom.num_tot() {
            let (sum_x, sum_y) = cent_sum[vid];
            self.cvol_skew[vid] = (sum_x / self.cvol_area[vid], sum_y / self.cvol_area[vid]);
        }

        // step 4: least-squares gradient weights
        // minimize sum(omega_n * (u_c + grad(u) . dx_n - u_n)^2) with omega_n = 1 / d_n^2
        // grad(u) = M^-1 sum(omega_n * dx_n * (u_n - u_c)), where M = sum(omega_n * dx_n dx_n^T)
        for vid_c in 0..dom.num_tot() {
            // get center point
            let (xc, yc) = pts_pos[vid_c];

            // construct moment matrix M
            let mut m_xx = 0.0;
            let mut m_xy = 0.0;
            let mut m_yy = 0.0;
            for &vid_n in &self.ngh_vid[vid_c] {
                let (dx, dy) = (pts_pos[vid_n].0 - xc, pts_pos[vid_n].1 - yc);
                let omega = 1.0 / (dx * dx + dy * dy);
                m_xx += omega * dx * dx;
                m_xy += omega * dx * dy;
                m_yy += omega * dy * dy;
            }

            // invert M; omega * dx^2 <= 1, so M is dimensionless
            let det = m_xx * m_yy - m_xy * m_xy;
            if !(det > 1e-12 * (m_xx + m_yy).powi(2)) { return Err(VFVChemError::TriangleSingularGradient { dom_id: dom.dom_id, vid: vid_c }) }
            let (i_xx, i_xy, i_yy) = (m_yy / det, -m_xy / det, m_xx / det);

            // get weights w_n = omega_n * M^-1 dx_n
            for &vid_n in &self.ngh_vid[vid_c] {
                let (dx, dy) = (pts_pos[vid_n].0 - xc, pts_pos[vid_n].1 - yc);
                let omega = 1.0 / (dx * dx + dy * dy);
                self.grad_wgt[vid_c].push((omega * (i_xx * dx + i_xy * dy), omega * (i_xy * dx + i_yy * dy)));
            }
        }

        Ok(())
    }

}
