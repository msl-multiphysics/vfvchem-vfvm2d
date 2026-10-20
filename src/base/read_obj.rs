use crate::base::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

struct ObjGroup {
    id: String,
    is_dom: bool,
    lines: Vec<(usize, Vec<usize>)>,  // (line number, vertex ids)
}

pub fn read_nds_obj(file_path: &str, h: f64, offset_ratio: f64) -> Result<Nodes, VFVChemError> {
    // read boundary, interface, and domain polylines from an obj file
    // every l record must belong to a group named bnd*, itf*, or dom*
    // nodes are spaced like Nodes::new_from_bounds and inset from vertices by dh
    // a segment shared by groups reuses the same node ids
    // domain normals point to the right of the traversal, out of a ccw domain
    // each domain line is one closed loop; a later loop is a hole, traversed clockwise

    // default struct
    let mut nds = Nodes::default();
    let dh = h * offset_ratio;
    nds.space = h;
    
    // file reader
    let file = File::open(file_path).map_err(|e| _obj_file(file_path, e))?;
    let reader = BufReader::new(file);

    // obj contents
    let mut verts = Vec::new();
    let mut groups: Vec<ObjGroup> = Vec::new();
    let mut cur: Option<usize> = None;

    // iterate over lines
    for (line_idx, line) in reader.lines().enumerate() {
        let line_no = line_idx + 1;
        let line = line.map_err(|e| _obj_file(file_path, e))?;
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() { continue; }

        let mut tok = line.split_whitespace();
        let tag = tok.next().unwrap();
        match tag {
            "v" => {
                let coords: Vec<&str> = tok.collect();
                if !(2..=4).contains(&coords.len()) {
                    return Err(_obj_parse(file_path, format!("line {line_no}: vertex must have x and y coordinates")));
                }
                let x = _obj_f64(file_path, line_no, coords[0])?;
                let y = _obj_f64(file_path, line_no, coords[1])?;
                verts.push((x, y));
            }
            "g" => {
                let Some(name) = tok.next() else {
                    return Err(_obj_parse(file_path, format!("line {line_no}: group must have one name")));
                };
                if tok.next().is_some() {
                    return Err(_obj_parse(file_path, format!("line {line_no}: group must have one name")));
                }
                let is_dom = if name.starts_with("dom") {
                    true
                } else if name.starts_with("bnd") || name.starts_with("itf") {
                    false
                } else {
                    return Err(_obj_parse(file_path, format!("line {line_no}: group name '{name}' must begin with 'bnd', 'itf', or 'dom'")));
                };

                // repeated declarations continue the same group
                if let Some(idx) = groups.iter().position(|grp| grp.id == name) {
                    cur = Some(idx);
                } else {
                    groups.push(ObjGroup { id: name.to_string(), is_dom, lines: Vec::new() });
                    cur = Some(groups.len() - 1);
                }
            }
            "l" => {
                let Some(idx) = cur else {
                    return Err(_obj_parse(file_path, format!("line {line_no}: line is not in a group")));
                };
                let mut vids = Vec::new();
                for tok in tok {
                    let vid: usize = tok.parse().map_err(|_| {
                        _obj_parse(file_path, format!("line {line_no}: invalid vertex index '{tok}'"))
                    })?;
                    if vid == 0 {
                        return Err(_obj_parse(file_path, format!("line {line_no}: vertex index must be at least 1")));
                    }
                    vids.push(vid);
                }
                if vids.len() < 2 {
                    return Err(_obj_parse(file_path, format!("line {line_no}: line must contain at least two vertices")));
                }
                groups[idx].lines.push((line_no, vids));
            }
            _ => {
                return Err(_obj_parse(file_path, format!("line {line_no}: unsupported record '{tag}'")));
            }
        }
    }

    // resolve 1-based indices after every vertex has been read
    for grp in &mut groups {
        for (line_no, vids) in &mut grp.lines {
            for vid in vids.iter_mut() {
                if *vid > verts.len() {
                    return Err(_obj_parse(file_path, format!("line {line_no}: vertex index {} is out of range", *vid)));
                }
                *vid -= 1;
            }
        }
    }

    // node ids for an unordered vertex pair, ordered from the lower index to the higher index
    let mut seg_gnd: HashMap<(usize, usize), Vec<usize>> = HashMap::new();

    // add nodes group by group, in file order
    for grp in &groups {
        let mut gnds = Vec::new();
        let mut loops = Vec::new();
        let mut norms = Vec::new();
        for (line_no, vids) in &grp.lines {
            let mut segs: Vec<(usize, usize)> = vids.windows(2).map(|w| (w[0], w[1])).collect();
            let is_closed = vids.first() == vids.last();
            if grp.is_dom && !is_closed {
                segs.push((*vids.last().unwrap(), *vids.first().unwrap()));
            }

            let mut ring = Vec::new();
            let mut ring_norm = Vec::new();
            for (a, b) in segs {
                let ids = _segment_gnd(&mut nds, &mut seg_gnd, &verts, a, b, h, dh, file_path, *line_no)?;
                if grp.is_dom {
                    let (nx, ny) = _right_normal(&verts, a, b);
                    ring_norm.resize(ring_norm.len() + ids.len(), (nx, ny));
                }
                ring.extend(ids);
            }
            if grp.is_dom {
                loops.push(ring);
                norms.push(ring_norm);
            } else {
                gnds.extend(ring);
            }
        }

        if grp.is_dom {
            if loops.is_empty() { continue; }
            nds.grp2d_gnd.insert(grp.id.clone(), loops);
            nds.grp2d_norm.insert(grp.id.clone(), norms);
        } else {
            if gnds.is_empty() { continue; }
            nds.grp1d_gnd.insert(grp.id.clone(), gnds);
        }
    }

    Ok(nds)
}

fn _segment_gnd(nds: &mut Nodes, seg_gnd: &mut HashMap<(usize, usize), Vec<usize>>, verts: &[(f64, f64)], a: usize, b: usize, h: f64, dh: f64, file_path: &str, line_no: usize) -> Result<Vec<usize>, VFVChemError> {
    // reject segments that cannot hold two inset endpoints
    let (x0, y0) = verts[a];
    let (x1, y1) = verts[b];
    let dx = x1 - x0;
    let dy = y1 - y0;
    let len = (dx * dx + dy * dy).sqrt();
    if len == 0.0 {
        return Err(_obj_parse(file_path, format!("line {line_no}: segment has zero length")));
    }
    if len <= 2.0 * dh {
        return Err(_obj_parse(file_path, format!("line {line_no}: segment length must be greater than two corner offsets. Got length = {len} and 2*dh = {}", 2.0 * dh)));
    }

    // reuse nodes when another group traverses the same vertices
    let key = (a.min(b), a.max(b));
    if let Some(gnd) = seg_gnd.get(&key) {
        if a <= b { return Ok(gnd.clone()); }
        return Ok(gnd.iter().rev().copied().collect());
    }

    // place endpoints dh from the vertices and space the interior like new_from_bounds
    // a short segment still gets one interval, so it has the two inset endpoints
    let (x0, y0) = verts[key.0];
    let (x1, y1) = verts[key.1];
    let dx = x1 - x0;
    let dy = y1 - y0;
    let tx = dx / len;
    let ty = dy / len;
    let num = (((len + 1e-6) / h) as usize).max(1);
    let ds = len / (num as f64);

    let mut pos = Vec::with_capacity(num + 1);
    pos.push((x0 + dh * tx, y0 + dh * ty));
    for i in 1..num {
        let s = (i as f64) * ds;
        pos.push((x0 + s * tx, y0 + s * ty));
    }
    pos.push((x1 - dh * tx, y1 - dh * ty));

    let g0 = nds.nds_pos.len();
    nds.nds_vel.extend(vec![(0.0, 0.0); pos.len()]);
    nds.nds_pos.extend(pos);
    let gnd: Vec<usize> = (g0..nds.nds_pos.len()).collect();
    let ordered = if a <= b { gnd.clone() } else { gnd.iter().rev().copied().collect() };
    seg_gnd.insert(key, gnd);
    Ok(ordered)
}

fn _right_normal(verts: &[(f64, f64)], a: usize, b: usize) -> (f64, f64) {
    // unit normal to the right of the directed segment
    let (x0, y0) = verts[a];
    let (x1, y1) = verts[b];
    let dx = x1 - x0;
    let dy = y1 - y0;
    let len = (dx * dx + dy * dy).sqrt();
    (dy / len, -dx / len)
}

fn _obj_f64(file_path: &str, line_no: usize, tok: &str) -> Result<f64, VFVChemError> {
    let val: f64 = tok.parse().map_err(|_| {
        _obj_parse(file_path, format!("line {line_no}: invalid number '{tok}'"))
    })?;
    if !val.is_finite() {
        return Err(_obj_parse(file_path, format!("line {line_no}: non-finite number '{tok}'")));
    }
    Ok(val)
}

fn _obj_file(file_path: &str, source: std::io::Error) -> VFVChemError {
    VFVChemError::ReadObjFile { file_path: file_path.to_string(), source }
}

fn _obj_parse(file_path: &str, reason: String) -> VFVChemError {
    VFVChemError::ReadObjParse { file_path: file_path.to_string(), reason }
}
