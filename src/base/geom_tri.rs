use crate::base::prelude::*;
use anyhow::Error;
use simple_delaunay_lib::delaunay_2d::delaunay_struct_2d::DelaunayStructure2D;
use simple_delaunay_lib::delaunay_2d::simplicial_struct_2d::Node;
use std::collections::HashMap;

const EPSILON: f64 = 1e-10;

#[derive(Default)]
pub struct Triangle {
    // ids
    pub tri_id: usize,
    pub dom_id: usize,  // domain this tri refers to

    // triangulation data
    // vid is the index of combined point and node vectors
    // i.e., vid = if vid < dom.num_pts() {pid} else {dom.num_pts() + nid}
    pub tri_vid: Vec<[usize; 3]>,  // vid counterclockwise

    // connectivity data
    pub ngh_vid: Vec<Vec<usize>>,  // vid -> neighbor vid
    pub ngh_dist: Vec<Vec<f64>>,  // vid -> edge length to neighbor

    // vertex data
    pub vert_area: Vec<f64>,  // vid (vertex id) -> CV (control volume) area
    pub vert_offc: Vec<(f64, f64)>,  // vid -> CV centroid minus vertex position (off-center)

    // face data
    pub face_len: Vec<Vec<f64>>,  // vid, loc (local neighbor index) -> signed face lengths; negative if circumcenters cross the edge
    pub face_norm: Vec<Vec<(f64, f64)>>,  // vid, loc -> outward normal on face
    pub face_skew: Vec<Vec<(f64, f64)>>,  // vid, loc -> face center minus edge midpoint (skewness)

    // boundary half-edge data
    // boundary vertices own half of the boundary edges on either side
    pub bnd_vid: Vec<Vec<usize>>,  // vid, blc (boundary loc) -> neighbor vid along boundary; empty for points
    pub bnd_half: Vec<Vec<f64>>,  // vid, blc -> half-edge length along boundary
    pub bnd_norm: Vec<Vec<(f64, f64)>>,  // vid, blc -> outward normal on half-edge
    pub bnd_tang: Vec<Vec<(f64, f64)>>,  // vid, blc -> unit tangent from vertex toward neighbor along boundary

    // gradient weight data
    // delta_u = [u1 - u0, u2 - u0, ...]; u0 and uj are point and neighbor values
    // du/dx = wgt_x . delta_u; du/dy = wgt_y . delta_u
    pub grad_wgt: Vec<Vec<(f64, f64)>>,  // vid, loc -> least-squares weights (wgt_x, wgt_y)
}

impl Triangle {
    pub fn new(tri_id: usize, pts: &Points, nds: &Nodes, dom: &Domain) -> Result<Triangle, VFVChemError> {
        // outline:
        // get triangulation and neighbors
        // compute per-triangle data
        // compute per-vertex data
        // compute gradient data
        
        // default struct
        let mut tri = Triangle::default();
        tri.tri_id = tri_id;
        tri.dom_id = dom.dom_id;

        // step 1: get triangulation

        // collect vertex positions
        let num_vert = dom.num_vert();
        let mut vert_pos = Vec::with_capacity(num_vert);
        for &gpt in dom.pts_gpt.iter() {
            let (x, y) = pts.pts_pos[gpt];
            vert_pos.push([x, y]);  // delaunay needs slice not tuple
        }
        for &gnd in dom.nds_gnd.iter() {
            let (x, y) = nds.nds_pos[gnd];
            vert_pos.push([x, y]);
        }

        // collect curve positions
        let mut curve_pos = Vec::new();
        for &eid in nds.grp2d_eid[&dom.grp2d_id].iter(){
            let mut edge_pos = [(0.0, 0.0); 2];
            for (i, &gvt) in nds.edge_gvt[eid].iter().enumerate() {
                let pos = nds.vert_pos[gvt];
                edge_pos[i] = pos;
            }
            curve_pos.push(edge_pos);
        }

        // get triangulation
        let tri_vid = _delaunay(&vert_pos, &curve_pos, dom.num_pts())?;

        // step 2: compute per-triangle data

        // initialize edge to geometry maps
        let mut fl_map = HashMap::new();  // (vida, vidb) -> face length
        let mut fm_map = HashMap::new();  // (vida, vidb) -> face length moment
        let mut ec_map = HashMap::new(); // (vida, vidb) -> (number of triangles, opposite vertex)

        // initialize vertex geometry
        let mut vert_area = vec![0.0; num_vert];  // vid -> vertex area (control volume)
        let mut vert_mom = vec![(0.0, 0.0); num_vert];  // vid -> vertex area moment about vertex

        // iterate through triangles
        for &vid_list in tri_vid.iter() {
            // get vertex positions
            let [vid0, vid1, vid2] = vid_list;
            let [x0, y0] = vert_pos[vid0];
            let [x1, y1] = vert_pos[vid1];
            let [x2, y2] = vert_pos[vid2];

            // get circumcenter
            let (ax, ay) = (x1 - x0, y1 - y0);
            let (bx, by) = (x2 - x0, y2 - y0);
            let det = ax * by - ay * bx;
            let amag2 = ax * ax + ay * ay;
            let bmag2 = bx * bx + by * by;
            let cx = x0 + 0.5 * (by * amag2 - ay * bmag2) / det;  // circumcenter
            let cy = y0 + 0.5 * (ax * bmag2 - bx * amag2) / det;

            // compute face and area pieces
            // iterate through edges
            for i in 0..3 {
                // get endpoint positions
                let vidi = vid_list[i];  // endpoint 1
                let vidj = vid_list[(i + 1) % 3];  // endpoint 2
                let vidk = vid_list[(i + 2) % 3];  // opposite
                let [xi, yi] = vert_pos[vidi];
                let [xj, yj] = vert_pos[vidj];
                let [xk, yk] = vert_pos[vidk];

                // compute midpoint
                let (mx, my) = (0.5 * (xi + xj), 0.5 * (yi + yj));

                // compute face pieces
                let dot = (xi - xk) * (xj - xk) + (yi - yk) * (yj - yk);
                let crs = ((xi - xk) * (yj - yk) - (yi - yk) * (xj - xk)).abs();
                let dist = ((xj - xi) * (xj - xi) + (yj - yi) * (yj - yi)).sqrt();
                let len = 0.5 * dot / crs * dist;  // length of face within the triangle
                let (fwx, fwy) = (0.5 * (cx - mx), 0.5 * (cy - my));  // weights
                let (fmx, fmy) = (len * fwx, len * fwy);  // face moment

                // add face pieces
                let pair = (vidi.min(vidj), vidi.max(vidj));
                let fl_val = fl_map.entry(pair).or_insert(0.0);
                *fl_val += len;
                let fm_val = fm_map.entry(pair).or_insert((0.0, 0.0));
                fm_val.0 += fmx;
                fm_val.1 += fmy;

                // compute vertex pieces
                // pieces on either side of the face have the same area
                // moments are taken about each endpoint
                let area = 0.25 * dist * len;  // area of triangle defined by len and dist
                let (vwix, vwiy) = ((xi + mx + cx) / 3.0 - xi, (yi + my + cy) / 3.0 - yi);  // weights for endpoint 1
                let (vwjx, vwjy) = ((xj + mx + cx) / 3.0 - xj, (yj + my + cy) / 3.0 - yj);  // weights for endpoint 2

                // add vertex pieces
                vert_area[vidi] += area;
                vert_area[vidj] += area;
                vert_mom[vidi].0 += area * vwix;
                vert_mom[vidi].1 += area * vwiy;
                vert_mom[vidj].0 += area * vwjx;
                vert_mom[vidj].1 += area * vwjy;

                // count triangles sharing the edge
                let ec_val = ec_map.entry(pair).or_insert((0, vidk));
                ec_val.0 += 1;
            }
        }

        // compute skewness from moments
        let mut fs_map = HashMap::new();  // (vida, vidb) -> face skewness
        for (pair, len) in fl_map.iter() {
            let (fmx, fmy) = fm_map[&pair];
            let skew = if len.abs() > 1e-10 {(fmx/len, fmy/len)} else {(0.0, 0.0)};
            fs_map.insert(pair, skew);
        }

        // compute off-center from moments
        let mut vert_offc = Vec::with_capacity(num_vert);
        for (vid, &area) in vert_area.iter().enumerate() {
            let (vmx, vmy) = vert_mom[vid];
            let offc = if area.abs() > 1e-10 {(vmx/area, vmy/area)} else {(0.0, 0.0)};
            vert_offc.push(offc);
        }

        // initialize vertex to geometry maps
        let mut ngh_vid = vec![Vec::new(); num_vert];
        let mut face_len = vec![Vec::new(); num_vert];
        let mut face_skew = vec![Vec::new(); num_vert];

        // iterate through triangles and edges
        for &vid_list in tri_vid.iter() {
        for i in 0..3 {
            // get endpoint positions
            let vida = vid_list[i];
            let vidb = vid_list[(i + 1) % 3];
            let pair = (vida.min(vidb), vida.max(vidb));

            // assign vidb as neighbor of vida and vice versa
            if !ngh_vid[vida].contains(&vidb) {
                ngh_vid[vida].push(vidb);
                face_len[vida].push(fl_map[&pair]);
                face_skew[vida].push(fs_map[&pair]);
            }
            if !ngh_vid[vidb].contains(&vida) {
                ngh_vid[vidb].push(vida);
                face_len[vidb].push(fl_map[&pair]);
                face_skew[vidb].push(fs_map[&pair]);
            }
        }}

        // initialize boundary half-edge maps
        let mut bnd_vid = vec![Vec::new(); num_vert];
        let mut bnd_half = vec![Vec::new(); num_vert];
        let mut bnd_norm = vec![Vec::new(); num_vert];
        let mut bnd_tang = vec![Vec::new(); num_vert];

        // iterate through edges belonging to only one triangle
        for (&(vida, vidb), &(num_tri, vidk)) in ec_map.iter() {
            // skip interior edges
            if num_tri != 1 { continue; }

            // get endpoint and opposite positions
            let [xa, ya] = vert_pos[vida];
            let [xb, yb] = vert_pos[vidb];
            let [xk, yk] = vert_pos[vidk];

            // compute outward normal
            // flip if pointing toward the opposite vertex
            let (dx, dy) = (xb - xa, yb - ya);
            let dist = (dx * dx + dy * dy).sqrt();
            let (mut nx, mut ny) = (dy / dist, -dx / dist);
            if nx * (xk - xa) + ny * (yk - ya) > 0.0 {
                nx = -nx;
                ny = -ny;
            }

            // assign half of edge to each endpoint
            bnd_vid[vida].push(vidb);
            bnd_half[vida].push(0.5 * dist);
            bnd_norm[vida].push((nx, ny));
            bnd_tang[vida].push((dx / dist, dy / dist));
            bnd_vid[vidb].push(vida);
            bnd_half[vidb].push(0.5 * dist);
            bnd_norm[vidb].push((nx, ny));
            bnd_tang[vidb].push((-dx / dist, -dy / dist));
        }

        // step 3: compute per-vertex data

        // initialize
        let mut face_norm = Vec::with_capacity(num_vert);
        let mut ngh_dist = Vec::with_capacity(num_vert);

        // iterate through vertices
        for (vid_v, list_n) in ngh_vid.iter().enumerate() {
            // get vertex coordinates
            let [xv, yv] = vert_pos[vid_v];
            
            // initialize subset
            let mut fn_list = Vec::new();
            let mut nd_list = Vec::new();

            // iterate through neighbors
            for &vid_n in list_n.iter() {
                // get neighbor coordinates
                let [xn, yn] = vert_pos[vid_n];   

                // compute distance and normals
                let (dx, dy) = (xn - xv, yn - yv);
                let dist = (dx * dx + dy * dy).sqrt();
                fn_list.push((dx/dist, dy/dist));
                nd_list.push(dist);
            }

            // store subset
            face_norm.push(fn_list);
            ngh_dist.push(nd_list);
        }

        // step 4: compute gradient data
        
        // initialize
        let mut grad_wgt = Vec::with_capacity(num_vert);

        // iterate through vertices
        for (vid_v, list_n) in ngh_vid.iter().enumerate() {
            // initialize subset
            let mut gw_list = Vec::new();
            
            // get vertex coordinates
            let [xv, yv] = vert_pos[vid_v];

            // assemble moment matrix (M)
            // M g = sum_n ((u_n - (u_v + g . dist_vn))^2 / dist_vn^2)
            // g - gradient; u_v + g . dist_vn - estimate of u_n from gradient
            let mut m_xx = 0.0;
            let mut m_xy = 0.0;
            let mut m_yy = 0.0;
            for &vid_n in list_n.iter() {
                let [xn, yn] = vert_pos[vid_n];
                let (dx, dy) = (xn - xv, yn - yv);
                let dist2 = dx * dx + dy * dy;
                m_xx += dx * dx / dist2;
                m_xy += dx * dy / dist2;
                m_yy += dy * dy / dist2;
            }

            // get inverse of M
            let det = m_xx * m_yy - m_xy * m_xy;
            let i_xx = m_yy / det;
            let i_xy = -m_xy / det;
            let i_yy = m_xx / det;

            // get gradient weights (w_x, w_y)
            // g_x = sum_n (w_x * (u_n - u_v))
            for &vid_n in list_n.iter() {
                let [xn, yn] = vert_pos[vid_n];
                let (dx, dy) = (xn - xv, yn - yv);
                let dist2 = dx * dx + dy * dy;
                let wx = (i_xx * dx + i_xy * dy) / dist2;
                let wy = (i_xy * dx + i_yy * dy) / dist2;
                gw_list.push((wx, wy));
            }

            // store subset
            grad_wgt.push(gw_list);
        }

        // store data
        tri.tri_vid = tri_vid;
        tri.ngh_vid = ngh_vid;
        tri.ngh_dist = ngh_dist;
        tri.vert_area = vert_area;
        tri.vert_offc = vert_offc;
        tri.face_len = face_len;
        tri.face_norm = face_norm;
        tri.face_skew = face_skew;
        tri.bnd_vid = bnd_vid;
        tri.bnd_half = bnd_half;
        tri.bnd_norm = bnd_norm;
        tri.bnd_tang = bnd_tang;
        tri.grad_wgt = grad_wgt;

        Ok(tri)
    }
}

fn _delaunay(vert_pos: &Vec<[f64; 2]>, curve_pos: &Vec<[(f64, f64); 2]>, num_pts: usize) -> Result<Vec<[usize; 3]>, VFVChemError> {
    // outline:
    // simple delaunay triangulation
    // remove triangles outside domain
    //
    // notes:
    // assume that circle around boundary segments has nothing inside it
    // simple delaunay triangulation should line the boundary segments 
    // nodes are offset from corners; triangles can only cross the boundary at concave corners
    // triangles with points are inside since points are inside the domain
    // triangles with only nodes are outside if two or more edges are outside
    // keeping triangles that cut concave corners ensures that the domain boundary only has nodes

    // step 1: delaunay triangulation

    // compute delaunay triangulation
    let mut dela = DelaunayStructure2D::new();
    dela.insert_vertices(&vert_pos, true).map_err(|e| _dela_err(e))?;
    let mesh = dela.get_simplicial();

    // step 2: remove triangles outside boundary

    // iterate through triangles
    let mut tri_vec = Vec::new();
    for tid in 0..mesh.get_nb_triangles() {
        // get triangle vertex ids
        let tri = mesh.get_triangle(tid).map_err(|e| _dela_err(e))?;
        if tri.contains_infinity() { continue; }  // skip ghost vertices
        let [nds0, nds1, nds2] = tri.nodes();
        let (Node::Value(vid0), Node::Value(vid1), Node::Value(vid2)) = (nds0, nds1, nds2) else { continue };
        let tri_vid = [vid0, vid1, vid2];

        // insert triangle if any vertex is a point
        if tri_vid.iter().any(|&vid| vid < num_pts) {
            tri_vec.push([vid0, vid1, vid2]);
            continue;
        }

        // get triangle vertex positions
        let mut tri_pos = [[0.0; 2]; 3];
        for i in 0..3 {
            let vid = tri_vid[i];
            tri_pos[i] = vert_pos[vid];
        }
        
        // iterate through edges and compute midpoint
        // count edges outside domain
        let mut num_out = 0;
        for i in 0..3 {
            // compute midpoint
            let [xa, ya] = tri_pos[i];
            let [xb, yb] = tri_pos[(i + 1) % 3];
            let xv = 0.5 * (xa + xb);
            let yv = 0.5 * (ya + yb);

            // check if inside
            if !_vertex_inside(curve_pos, xv, yv) {
                num_out += 1;
            }
        }

        // insert triangle if inside or if it cuts a concave corner
        if num_out <= 1 {
            tri_vec.push([vid0, vid1, vid2]);
        }
    }

    Ok(tri_vec)
}

fn _vertex_inside(curve_pos: &Vec<[(f64, f64); 2]>, xv: f64, yv: f64) -> bool {
    // winding number algorithm
    // +1 for ccw curves and -1 for cw curves
    // point is inside if winding number is positive or if point is on a segment

    // initialize winding number
    let mut wind = 0;

    // iterate through curves
    for &[(x0, y0), (x1, y1)] in curve_pos {
        // cross product of (xv - x0) and (x1 - x0)
        // positive if xv on left of segment
        let cross = (x1 - x0) * (yv - y0) - (xv - x0) * (y1 - y0); 

        // point on segment is inside
        let dot = (x1 - x0) * (xv - x0) + (y1 - y0) * (yv - y0);
        let len2 = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
        if cross.abs() <= EPSILON * len2 && dot >= 0.0 && dot <= len2 { return true; }

        // CCW - point on left; segment goes up
        // CW - point on right; segment goes down
        if y0 <= yv && y1 > yv && cross > 0.0 { wind += 1; }
        if y0 > yv && y1 <= yv && cross < 0.0 { wind -= 1; }
    }

    // return
    wind > 0
}

fn _dela_err(source: Error) -> VFVChemError {
    VFVChemError::TriangleError { source }
}
