#![allow(dead_code)]
use vfvchem_vfvm2d::*;

// shared tools for the benchmarks in examples/benchmark
// run with: `cargo run --release --example <name>`

pub fn vert_pos(vars: &Variables, dom_id: usize) -> Vec<(f64, f64)> {
    // vid -> (x, y); points first, then nodes
    let dom = &vars.dom[dom_id];
    let mut pos = Vec::with_capacity(dom.num_vert());
    for &gpt in dom.pts_gpt.iter() {
        pos.push(vars.pts.pts_pos[gpt]);
    }
    for &gnd in dom.nds_gnd.iter() {
        pos.push(vars.nds.nds_pos[gnd]);
    }
    pos
}

pub fn bnd_vid(vars: &Variables, bnd_id: usize) -> Vec<usize> {
    // boundary index -> vid in the domain of the boundary
    let bnd = &vars.bnd[bnd_id];
    let dom = &vars.dom[bnd.dom_id];
    bnd.nds_gnd.iter().map(|gnd| dom.num_pts() + dom.nds_nid[gnd]).collect()
}

pub fn itf_vid(vars: &Variables, itf_id: usize, dom_id: usize) -> Vec<usize> {
    // interface index -> vid in the selected domain
    let itf = &vars.itf[itf_id];
    let dom = &vars.dom[dom_id];
    itf.nds_gnd.iter().map(|gnd| dom.num_pts() + dom.nds_nid[gnd]).collect()
}

pub fn err_norm(err: &[f64], area: &[f64]) -> (f64, f64, f64) {
    // area-weighted L1 and L2 norms (divided by total area) and Linf norm
    let mut sum_a = 0.0;
    let mut sum_1 = 0.0;
    let mut sum_2 = 0.0;
    let mut max_e: f64 = 0.0;
    for (e, a) in err.iter().zip(area.iter()) {
        sum_a += a;
        sum_1 += e.abs() * a;
        sum_2 += e * e * a;
        max_e = max_e.max(e.abs());
    }
    (sum_1 / sum_a, (sum_2 / sum_a).sqrt(), max_e)
}

pub fn get_order(err_c: f64, err_f: f64, h_c: f64, h_f: f64) -> f64 {
    // observed order between a coarse and fine mesh
    (err_c / err_f).ln() / (h_c / h_f).ln()
}

pub fn fit_order(err: &[f64], h: &[f64]) -> f64 {
    // least-squares slope of log(err) against log(h)
    let n = err.len() as f64;
    let lx: Vec<f64> = h.iter().map(|v| v.ln()).collect();
    let ly: Vec<f64> = err.iter().map(|v| v.ln()).collect();
    let mx = lx.iter().sum::<f64>() / n;
    let my = ly.iter().sum::<f64>() / n;
    let sxy: f64 = lx.iter().zip(ly.iter()).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = lx.iter().map(|x| (x - mx) * (x - mx)).sum();
    sxy / sxx
}

pub const EPS: f64 = f64::EPSILON;  // machine epsilon for round-off tolerances

pub fn mat_size(vars: &Variables) -> usize {
    // size of the linear system from the unknown start indices
    let mut size = 0;
    for scldom in vars.scldom.iter() {
        if let SclDomType::Unknown { start } = scldom.scldom_type { size = size.max(start + scldom.value.len()); }
    }
    for vecdom in vars.vecdom.iter() {
        if let VecDomType::Unknown { start } = vecdom.vecdom_type { size = size.max(start + 2 * vecdom.value.len()); }
    }
    size
}

pub fn get_res(vars: &Variables, oper_list: &[&dyn OperatorBase]) -> (Vec<f64>, Vec<f64>) {
    // residual r = A x - b of the selected operators at the current solution
    // r of domain operators at a boundary node is minus the inward boundary flux
    // also returns the round-off scale s = sum_j |a_ij x_j| + |b_i| of each row
    let size = mat_size(vars);
    let mut a_triplet = Vec::new();
    let mut b_vec = vec![0.0; size];
    for oper in oper_list.iter() {
        oper.apply(vars, &mut a_triplet, &mut b_vec).unwrap();
    }

    // collect solution vector
    let mut x_vec = vec![0.0; size];
    for scldom in vars.scldom.iter() {
        if let SclDomType::Unknown { start } = scldom.scldom_type {
            for (vid, &value) in scldom.value.iter().enumerate() { x_vec[start + vid] = value; }
        }
    }
    for vecdom in vars.vecdom.iter() {
        if let VecDomType::Unknown { start } = vecdom.vecdom_type {
            let num = vecdom.value.len();
            for (vid, &(vx, vy)) in vecdom.value.iter().enumerate() { x_vec[start + vid] = vx; x_vec[start + num + vid] = vy; }
        }
    }

    // compute residual
    let mut r_vec: Vec<f64> = b_vec.iter().map(|b| -b).collect();
    let mut s_vec: Vec<f64> = b_vec.iter().map(|b| b.abs()).collect();
    for &(r, c, v) in a_triplet.iter() {
        r_vec[r] += v * x_vec[c];
        s_vec[r] += (v * x_vec[c]).abs();
    }
    (r_vec, s_vec)
}

pub fn scl_start(vars: &Variables, scldom_id: usize) -> usize {
    match vars.scldom[scldom_id].scldom_type {
        SclDomType::Unknown { start } => start,
        _ => panic!("not an unknown"),
    }
}

pub fn vec_start(vars: &Variables, vecdom_id: usize) -> usize {
    match vars.vecdom[vecdom_id].vecdom_type {
        VecDomType::Unknown { start } => start,
        _ => panic!("not an unknown"),
    }
}

pub fn check(name: &str, pass: bool, detail: String, fail_list: &mut Vec<String>) {
    // print and record a pass/fail criterion
    let flag = if pass { "PASS" } else { "FAIL" };
    println!("[{flag}] {name}: {detail}");
    if !pass { fail_list.push(format!("{name}: {detail}")); }
}
