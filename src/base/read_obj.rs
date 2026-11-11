use crate::base::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::SplitWhitespace;

pub fn read_nds_obj(path: &str, nds_space: f64, off_space: f64, ref_space: f64, ref_ratio: f64) -> Result<Nodes, VFVChemError> {
    // outline:
    // collect vertices, edges, and groupings
    // create nodes for each edge
    
    // step 1: collect vertices, edges, and groupings

    // initialize
    let mut vert_pos = Vec::new();
    let mut edge_gvt = Vec::new();
    let mut gvt_edge = HashMap::new();
    let mut grp1d_eid = HashMap::new();
    let mut grp2d_eid = HashMap::new();
    let mut group_str = "".to_string();

    // file reader
    let file = File::open(path).map_err(|e| _obj_err(path, e))?;
    let r = BufReader::new(file);

    // iterate over lines
    for line in r.lines() {
        // remove comments and skip blank lines
        let line = line.map_err(|e| _obj_err(path, e))?;
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() { continue; }

        // split at whitespace and iterate
        // process based on tag (first entry)
        let mut iter = line.split_whitespace();
        let tag = iter.next().unwrap_or("");
        match tag {
            "v" => { _get_vertex(&mut iter, &mut vert_pos, path)?; }
            "g" => { _get_group(&mut iter, &mut group_str, path)?; }
            "l" => { _get_edge(&mut iter, &mut edge_gvt, &mut gvt_edge, &mut grp1d_eid, &mut grp2d_eid, &group_str, path)?; }
            _ => { return Err(_obj_parse(path, "Encountered invalid tag.")); }
        }
    }

    // step 2: create nodes for each edge

    // initialize
    let mut nds_pos = Vec::new();
    let mut edge_gnd: Vec<Vec<usize>> = Vec::new();

    // iterate through edges
    for (eid, &[vid0, vid1]) in edge_gvt.iter().enumerate() {
        // if this edge is the reverse of an exisitng one, load nodes in reverse
        let pair_rev = [vid1, vid0];
        if gvt_edge.contains_key(&pair_rev) && gvt_edge[&pair_rev] < eid {
            let eid_rev = gvt_edge[&pair_rev];
            let node_id = edge_gnd[eid_rev].iter().rev().copied().collect();
            edge_gnd.push(node_id);
            continue;
        }

        // create nodes along edge

        // construct nodes
        let node_add = _nodes_line(vert_pos[vid0], vert_pos[vid1], nds_space, off_space, ref_space, ref_ratio);
        let mut node_id = Vec::new();
        for &pos in node_add.iter() {
            node_id.push(nds_pos.len());
            nds_pos.push(pos);
        }
        edge_gnd.push(node_id);
    }

    // store in struct
    let mut nds = Nodes::default();
    nds.nds_pos = nds_pos;
    nds.nds_space = nds_space;
    nds.vert_pos = vert_pos;
    nds.edge_gvt = edge_gvt;
    nds.edge_gnd = edge_gnd;
    nds.grp1d_eid = grp1d_eid;
    nds.grp2d_eid = grp2d_eid;

    Ok(nds)
}

fn _get_vertex(iter: &mut SplitWhitespace, vert_pos: &mut Vec<(f64, f64)>, path: &str) -> Result<(), VFVChemError> {
    // get coordinates
    let x_str = iter.next().ok_or_else(|| _obj_parse(path, "Could not find x coordinate."))?;
    let y_str = iter.next().ok_or_else(|| _obj_parse(path, "Could not find y coordinate."))?;
    let x: f64 = x_str.parse::<f64>().map_err(|_| _obj_parse(path, "Could not parse x coordinate into number."))?;
    let y: f64 = y_str.parse::<f64>().map_err(|_| _obj_parse(path, "Could not parse y coordinate into number."))?;
    vert_pos.push((x, y));

    Ok(())
}

fn _get_group(iter: &mut SplitWhitespace, group_str: &mut String, path: &str) -> Result<(), VFVChemError> {
    // reject group id if not starting with dom, bnd, or itf
    let grp_id: String = iter.collect();
    if !grp_id.starts_with("dom") && !grp_id.starts_with("bnd") && !grp_id.starts_with("itf") {
        return Err(_obj_parse(path, "Group name must start with dom, bnd, or itf."))
    }
    *group_str = grp_id;

    Ok(())
}

fn _get_edge(iter: &mut SplitWhitespace, edge_gvt: &mut Vec<[usize; 2]>, gvt_edge: &mut HashMap<[usize; 2], usize>, grp1d_eid: &mut HashMap<String, Vec<usize>>, grp2d_eid: &mut HashMap<String, Vec<usize>>, group_str: &String, path: &str) -> Result<(), VFVChemError> {
    // collect vertex ids
    let curve_str: Vec<&str> = iter.collect();
    let mut curve_vid = Vec::new();
    for &vid_str in curve_str.iter() {
        let vid_p1 = vid_str.parse::<usize>().map_err(|_| _obj_parse(path, "Could not parse vertex ID into number."))?;
        curve_vid.push(vid_p1 - 1);  // obj files start with 1
    }

    // assign edges to groups
    let mut curve_eid = Vec::new();
    for i in 0..curve_vid.len() - 1 {
        // get vertex index
        let vid0 = curve_vid[i];
        let vid1 = curve_vid[i + 1];
        let pair = [vid0, vid1];

        // get edge index corresponding to [vid0, vid1]
        // use index in edge_gvt if edge already exists, otherwise insert
        // note: [vid0, vid1] is different from [vid1, vid0]
        let eid = *gvt_edge.entry(pair).or_insert_with(|| {
            edge_gvt.push(pair);
            edge_gvt.len() - 1
        });

        // store edge index
        curve_eid.push(eid);
    }

    // add new group to maps
    if group_str.starts_with("dom") {
        grp2d_eid.entry(group_str.clone()).or_insert_with(Vec::new).extend(curve_eid);
    } else {
        grp1d_eid.entry(group_str.clone()).or_insert_with(Vec::new).extend(curve_eid);
    }
    
    Ok(())
}

fn _nodes_line(pos0: (f64, f64), pos1: (f64, f64), nds_space: f64, off_space: f64, ref_space: f64, ref_ratio: f64) -> Vec<(f64, f64)> {
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
    let b = ref_space * nds_space;
    let n = ((1.0 + (b - a) * (ref_ratio - 1.0) / a).ln() / ref_ratio.ln()) as i32;
    let clear_cnr = 0.5 * len_s - a;  // available space for refinement
    for i in 0..=n {
        let c = a + (b - a) * (ref_ratio.powi(i) - 1.0) / (ref_ratio.powi(n) - 1.0);
        if c > clear_cnr { break; }
        space_cnr.push(c);
    }
    let start_cnr = space_cnr[space_cnr.len() - 1];

    // compute regularly spaced node positions
    let clear_reg = len_s - 2.0 * b;  // amount of space for regular spacing
    let num_reg = ((clear_reg / nds_space) as usize).saturating_sub(1);  // number of regularly spaced points
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

fn _obj_err(path: &str, source: std::io::Error) -> VFVChemError {
    VFVChemError::ReadObjError { path: path.to_string(), source }
}

fn _obj_parse(path: &str, message: &str) -> VFVChemError {
    VFVChemError::ReadObjParse { path: path.to_string(), message: message.to_string() }
}
