use crate::base::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn write_group_vtu(vars: &Variables) -> Result<(), VFVChemError> {
    // initialize file writer
    let path = Path::new(vars.folder.as_str()).join("group.vtu");
    let file = File::create(&path).map_err(|e| _vtu_err(&path, e))?;
    let mut w = BufWriter::new(file);
 
    // compute number of points (vertices) and cells (triangles)
    // double-count points along interfaces
    let mut num_vert = 0;
    let mut num_tri = 0;
    let mut offset_vec = Vec::new();  // dom_id -> cumulative vertex count
    for dom in vars.dom.iter() {
        offset_vec.push(num_vert);
        num_vert += dom.num_vert();
    }
    for tri in vars.tri.iter() {
        num_tri += tri.tri_vid.len();
    }

    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", num_vert, num_tri).map_err(|e| _vtu_err(&path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&path, e))?;
    for dom in vars.dom.iter() {
        for &gpt in &dom.pts_gpt {
            let (x, y) = vars.pts.pts_pos[gpt];
            writeln!(w, "{:.9} {:.9} 0.0", x, y).map_err(|e| _vtu_err(&path, e))?;
        }
        for &gnd in &dom.nds_gnd {
            let (x, y) = vars.nds.nds_pos[gnd];
            writeln!(w, "{:.9} {:.9} 0.0", x, y).map_err(|e| _vtu_err(&path, e))?;
        }
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(&path, e))?;
   
    // write cells
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(&path, e))?;
    for tri in vars.tri.iter() {
        let offset = offset_vec[tri.tri_id];
        for &[vid0, vid1, vid2] in tri.tri_vid.iter() {
            writeln!(w, "{} {} {}", vid0 + offset, vid1 + offset, vid2 + offset).map_err(|e| _vtu_err(&path, e))?;
        }
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(&path, e))?;
    for i in 0..num_tri {
        writeln!(w, "{}", 3 * (i + 1)).map_err(|e| _vtu_err(&path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(&path, e))?;
    for _ in 0..num_tri {
        writeln!(w, "5").map_err(|e| _vtu_err(&path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(&path, e))?;
    
    // write point data
    writeln!(w, "<PointData>").map_err(|e| _vtu_err(&path, e))?;
    for dom in vars.dom.iter() {
        // create flag
        let offset = offset_vec[dom.dom_id];
        let mut flag_vec = vec![0; num_vert];
        for vid in 0..dom.num_vert() {
            flag_vec[offset + vid] = 1;
        }

        // write flag
        writeln!(w, "<DataArray type=\"Int32\" Name=\"{}\" format=\"ascii\">", dom.grp2d_id).map_err(|e| _vtu_err(&path, e))?;
        for flag in flag_vec.iter() {
            writeln!(w, "{}", flag).map_err(|e| _vtu_err(&path, e))?;
        }
        writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?; 
    }
    for bnd in vars.bnd.iter() {
        // create flag
        let dom = &vars.dom[bnd.dom_id];
        let offset = offset_vec[dom.dom_id];
        let mut flag_vec = vec![0; num_vert];
        for gnd in bnd.nds_gnd.iter() {
            let vid = dom.num_pts() + dom.nds_nid[gnd];
            flag_vec[offset + vid] = 1;
        }
        
        // write flag
        writeln!(w, "<DataArray type=\"Int32\" Name=\"{}\" format=\"ascii\">", bnd.grp1d_id).map_err(|e| _vtu_err(&path, e))?;
        for flag in flag_vec.iter() {
            writeln!(w, "{}", flag).map_err(|e| _vtu_err(&path, e))?;
        }
        writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?; 
    }
    for itf in vars.itf.iter() {
        // create flag
        let dom_a = &vars.dom[itf.dom_a_id];
        let dom_b = &vars.dom[itf.dom_b_id];
        let offset_a = offset_vec[dom_a.dom_id];
        let offset_b = offset_vec[dom_b.dom_id];
        let mut flag_vec = vec![0; num_vert];
        for gnd in itf.nds_gnd.iter() {
            let vid_a = dom_a.num_pts() + dom_a.nds_nid[gnd];
            let vid_b = dom_b.num_pts() + dom_b.nds_nid[gnd];
            flag_vec[offset_a + vid_a] = 1;
            flag_vec[offset_b + vid_b] = 1;
        }
        
        // write flag
        writeln!(w, "<DataArray type=\"Int32\" Name=\"{}\" format=\"ascii\">", itf.grp1d_id).map_err(|e| _vtu_err(&path, e))?;
        for flag in flag_vec.iter() {
            writeln!(w, "{}", flag).map_err(|e| _vtu_err(&path, e))?;
        }
        writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?; 
    }
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(&path, e))?;

    // write cell data
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(&path, e))?;    
    
    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(&path, e))?;

    Ok(())
}

pub fn write_vars_vtu(vars: &Variables, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let path = Path::new(vars.folder.as_str()).join(format!("vars_{:04}.vtu", ts));
    let file = File::create(&path).map_err(|e| _vtu_err(&path, e))?;
    let mut w = BufWriter::new(file);
 
    // compute number of points (vertices) and cells (triangles)
    // double-count points along interfaces
    let mut num_vert = 0;
    let mut num_tri = 0;
    let mut offset_vec = Vec::new();  // dom_id -> cumulative vertex count
    for dom in vars.dom.iter() {
        offset_vec.push(num_vert);
        num_vert += dom.num_vert();
    }
    for tri in vars.tri.iter() {
        num_tri += tri.tri_vid.len();
    }

    // group variables with similar names
    let mut scldom_grp = HashMap::new();  // name -> scldom_ids
    let mut vecdom_grp = HashMap::new();  // name -> scldom_ids
    for scldom in vars.scldom.iter() {
        let scldom_vec = scldom_grp.entry(scldom.name.clone()).or_insert(vec![]);
        scldom_vec.push(scldom.scldom_id);
    }
    for vecdom in vars.vecdom.iter() {
        let vecdom_vec = vecdom_grp.entry(vecdom.name.clone()).or_insert(vec![]);
        vecdom_vec.push(vecdom.vecdom_id);
    }

    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", num_vert, num_tri).map_err(|e| _vtu_err(&path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&path, e))?;
    for dom in vars.dom.iter() {
        for &gpt in &dom.pts_gpt {
            let (x, y) = vars.pts.pts_pos[gpt];
            writeln!(w, "{:.9} {:.9} 0.0", x, y).map_err(|e| _vtu_err(&path, e))?;
        }
        for &gnd in &dom.nds_gnd {
            let (x, y) = vars.nds.nds_pos[gnd];
            writeln!(w, "{:.9} {:.9} 0.0", x, y).map_err(|e| _vtu_err(&path, e))?;
        }
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(&path, e))?;
   
    // write cells
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(&path, e))?;
    for tri in vars.tri.iter() {
        let offset = offset_vec[tri.tri_id];
        for &[vid0, vid1, vid2] in tri.tri_vid.iter() {
            writeln!(w, "{} {} {}", vid0 + offset, vid1 + offset, vid2 + offset).map_err(|e| _vtu_err(&path, e))?;
        }
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(&path, e))?;
    for i in 0..num_tri {
        writeln!(w, "{}", 3 * (i + 1)).map_err(|e| _vtu_err(&path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(&path, e))?;
    for _ in 0..num_tri {
        writeln!(w, "5").map_err(|e| _vtu_err(&path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(&path, e))?;
    
    // write point data
    writeln!(w, "<PointData>").map_err(|e| _vtu_err(&path, e))?;
    for (name, scldom_list) in scldom_grp.iter() {
        // create vector with values
        let mut value_vec = vec![f64::NAN; num_vert];
        for &scldom_id in scldom_list.iter() {
            let scldom = &vars.scldom[scldom_id];
            let offset = offset_vec[scldom.dom_id];
            for (vid, &value) in scldom.value.iter().enumerate() {
                value_vec[vid + offset] = value;
            }
        }

        // write values
        writeln!(w, "<DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">", name).map_err(|e| _vtu_err(&path, e))?;
        for value in value_vec.iter() {
            writeln!(w, "{}", value).map_err(|e| _vtu_err(&path, e))?;
        }
        writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?; 
    }
    for (name, vecdom_list) in vecdom_grp.iter() {
        // create vector with values
        let mut value_vec = vec![(f64::NAN, f64::NAN); num_vert];
        for &vecdom_id in vecdom_list.iter() {
            let vecdom = &vars.vecdom[vecdom_id];
            let offset = offset_vec[vecdom.dom_id];
            for (vid, &value) in vecdom.value.iter().enumerate() {
                value_vec[vid + offset] = value;
            }
        }

        // write values
        writeln!(w, "<DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"3\" format=\"ascii\">", name).map_err(|e| _vtu_err(&path, e))?;
        for &(value_x, value_y) in value_vec.iter() {
            writeln!(w, "{} {} 0.0", value_x, value_y).map_err(|e| _vtu_err(&path, e))?;
        }
        writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&path, e))?; 
    }
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(&path, e))?;

    // write cell data
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(&path, e))?;    
    
    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(&path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(&path, e))?;

    Ok(())
}

fn _vtu_err(path: &PathBuf, source: std::io::Error) -> VFVChemError {
    VFVChemError::WriteVtuError {path: path.to_path_buf(), source}
}
