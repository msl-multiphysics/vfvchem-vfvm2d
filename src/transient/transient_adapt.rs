use crate::base::prelude::*;
use crate::linear::prelude::*;
use crate::transient::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

const GROWTH_MIN: f64 = 0.25;
const GROWTH_MAX: f64 = 4.00;
const GROWTH_FAIL: f64 = 0.50;
const TOL_TIME: f64 = 100.0;

pub trait TransientAdaptive: TransientBase {
    // adaptive solver with default options
    fn solve_adaptive(
        &mut self, vars: &mut Variables, lin: &mut Box<dyn LinearBase>,
        max_iter: usize, tol: f64, damp: f64, dt_init: f64, dt_min: f64, dt_max: f64,
        t_end: f64, num_step: usize,
    ) -> Result<(), VFVChemError> {
        // create output times
        let dt = t_end / (num_step as f64);
        let t_out: Vec<f64> = (0..num_step + 1).map(|i| dt * (i as f64)).collect();

        // run adaptive solver
        self.solve_adaptive_options(
            vars, lin, max_iter, tol, damp, TOL_TIME * tol, dt_init, dt_min, dt_max,
            t_out, GROWTH_MIN, GROWTH_MAX, GROWTH_FAIL,
        )?;

        Ok(())
    }
    
    // adaptive solver with complete options
    // tol_iter is the nonlinear tolerance; tol_time is the local error tolerance
    // tol_iter should be well below tol_time; otherwise the nonlinear error is mistaken for local error
    fn solve_adaptive_options(
        &mut self, vars: &mut Variables, lin: &mut Box<dyn LinearBase>,
        max_iter: usize, tol_iter: f64, damp: f64, tol_time: f64, dt_init: f64, dt_min: f64, dt_max: f64,
        t_out: Vec<f64>, grow_min: f64, grow_max: f64, grow_fail: f64
    ) -> Result<(), VFVChemError> {
        // error handling
        if max_iter == 0 { return Err(VFVChemError::SolverNonPositiveMaxIter { max_iter }); }
        if tol_iter <= 0.0 { return Err(VFVChemError::SolverNonPositiveTolerance { tol: tol_iter }); }
        if tol_time <= 0.0 { return Err(VFVChemError::SolverNonPositiveTolerance { tol: tol_time }); }
        if damp <= 0.0 || damp > 1.0 { return Err(VFVChemError::SolverInvalidDamping { damp }); }
        if dt_init <= 0.0 { return Err(VFVChemError::SolverNonPositiveInitialStepSize { dt_init }) };
        if dt_init < dt_min { return Err(VFVChemError::SolverInvalidMinimumStepSize { dt_min }) };
        if dt_max <= dt_min { return Err(VFVChemError::SolverInvalidMaximumStepSize { dt_max }) };
        if grow_min <= 0.0 || grow_min > 1.0 { return Err(VFVChemError::SolverInvalidGrowthMinimum { grow_min }) };
        if grow_max <= 1.0 { return Err(VFVChemError::SolverInvalidGrowthMaximum { grow_max }) };
        if grow_fail <= 0.0 || grow_fail >= 1.0 { return Err(VFVChemError::SolverInvalidGrowthFail { grow_fail }) };

        // error handling of t_out
        if t_out.is_empty() { return Err(VFVChemError::SolverAdaptiveInvalidOutputTime {}); }  // t_out empty
        for i in 0..t_out.len() - 1 {
            if !(t_out[i] >= 0.0) || (t_out[i] >= t_out[i + 1]) {  // t_out positive and monotonically increasing
                return Err(VFVChemError::SolverAdaptiveInvalidOutputTime {});
            }
        }
        let t_start = t_out.iter().copied().find(|&t| t > 0.0).unwrap_or(0.0);
        if t_start < dt_init { return Err(VFVChemError::SolverAdaptiveInvalidOutputTime {}); };  // first t_out less than dt_init

        // initialize time measurement
        println!("Starting transient solver.");
        let time_start = Instant::now();
        let mut time_asm = Duration::ZERO;  // matrix assembly, value updating, etc.
        let mut time_lin = Duration::ZERO;  // linear solve
        let mut time_wrt = Duration::ZERO;  // file writing

        // write groups file
        let time_0 = Instant::now();
        write_group_vtu(vars)?;
        let time_1 = Instant::now();
        time_wrt += time_1 - time_0;

        // initialize matrices
        let mat_size = self.initial_start(vars)?;
        let mut a_triplet: Vec<(usize, usize, f64)> = Vec::new();
        let mut b_vec: Vec<f64> = vec![0.0; mat_size];
        let mut x_iter_vec: Vec<f64> = vec![0.0; mat_size];  // damped result

        // initialize operators
        self.initial_operator(vars)?;

        // load dirichlet values into unknown fields
        _update_func(vars, 0.0);  // update function-type variables
        self.initial_dirichlet(vars, mat_size)?;  // set dirichlet values

        // initial assembly of A and b
        // used for sparsity pattern
        let time_2 = Instant::now();
        _unkn_to_solu(vars, &mut x_iter_vec);  // unknown fields -> solution vector
        _update_func(vars, 0.0);
        _update_grad(vars);
        _update_prev(vars);
        self.assemble_matrix(vars, &mut a_triplet, &mut b_vec, dt_init)?;
        lin.set_sparse(&a_triplet, mat_size)?;
        let num_triplet = a_triplet.len();
        let time_3 = Instant::now();
        time_asm += time_3 - time_2;

        // write initial state if needed
        let mut num_write = 0;  // number of written output files
        if t_out[0] == 0.0 {
            let time_t0 = Instant::now();
            write_vars_vtu(vars, 0)?;
            let time_t1 = Instant::now();
            time_wrt += time_t1 - time_t0;
            num_write += 1;
        }

        // initialize for adaptive time stepping
        let mut t = 0.0;  // current time
        let mut dt = dt_init;  // current step size
        let t_end = t_out[t_out.len() - 1];  // final time
        let mut dt_prev=  0.0;  // step size between x_n and x_(n+1)
        let mut x_n1_vec = x_iter_vec.clone();  // x_n
        let mut x_n2_vec = x_iter_vec.clone();  // x_(n-1)
        let mut err_prev: f64 = 1.0;  // local error of the previous accepted step
        let mut t_write = t_out.iter().copied().find(|&ti| ti > t).unwrap(); // next time to write output

        // iterate until final time
        while t < t_end {
            // repeat until next solution is acceptable
            let mut is_accept = false;
            while !is_accept {
                // adjust dt based on t_out
                let mut t_try = t + dt;
                if t_try > t_write {
                    dt = t_write - t;
                    t_try = t_write;
                }

                // attempt a solution with step size dt
                
                // assembly of A and b at time t
                // used for first iteration
                let time_t0 = Instant::now();
                a_triplet = Vec::with_capacity(num_triplet);  // reset A
                b_vec = vec![0.0; mat_size];  // reset b
                _update_func(vars, t_try);
                _update_grad(vars);
                self.assemble_matrix(vars, &mut a_triplet, &mut b_vec, dt)?;
                let time_t1 = Instant::now();
                time_asm += time_t1 - time_t0;

                // iterate to convergence
                let mut is_conv = false;
                for iter in 0..max_iter {
                    // solve for undamped solution
                    // x_udmp in A_k * x_udmp = b_k
                    let time_i0 = Instant::now();
                    let x_udmp_vec = lin.solve(&a_triplet, &b_vec, &x_iter_vec)?;
                    let time_i1 = Instant::now();
                    time_lin += time_i1 - time_i0;

                    // solve for damped solution
                    // x_(k+1) in x_damp = (1 - damp) * x_k + damp * x_udmp
                    let x_damp_vec: Vec<f64> = (0..mat_size).map(|i| (1.0 - damp) * x_iter_vec[i] + damp * x_udmp_vec[i]).collect();
                    
                    // assemble A_(k+1) and b_(k+1) using x_(k+1)
                    let time_i2 = Instant::now();
                    a_triplet = Vec::with_capacity(num_triplet);  // reset A
                    b_vec = vec![0.0; mat_size];  // reset b
                    _solu_to_unkn(vars, &x_damp_vec);
                    _update_func(vars, t_try);
                    _update_grad(vars);
                    self.assemble_matrix(vars, &mut a_triplet, &mut b_vec, dt)?;
                    let time_i3 = Instant::now();
                    time_asm += time_i3 - time_i2;

                    // compute tolerance
                    let res = _solu_res(vars, &a_triplet, &b_vec, &x_damp_vec);
                    println!("Time: {t:.4e}; Step size: {dt:.4e}; Iteration: {iter}; Residual: {res:.4e}");

                    // reject if diverged
                    if !res.is_finite() {
                        is_conv = false;
                        break;
                    }

                    // update iteration vector
                    x_iter_vec = x_damp_vec;

                    // tolerance check
                    if res < tol_iter {
                        is_conv = true;
                        break;
                    }
                }

                // adjust dt based on result

                // if not converged at dt_min, then raise error
                if !is_conv && dt <= dt_min {
                    return Err(VFVChemError::SolverAdaptiveDidNotConverge { dt_min });
                }

                // if not converged, then shrink step size
                if !is_conv {
                    dt = (dt * grow_fail).max(dt_min);
                    _solu_to_unkn(vars, &x_n1_vec);
                    x_iter_vec = x_n1_vec.clone();
                    println!("Nonlinear solver did not converge. Reducing step size to {dt:.4e}.");
                    continue;
                }

                // if converged, adjust step size based on predictor

                // local error from the linear extrapolation predictor
                // x_pred = x_n + (dt / dt_prev) * (x_n - x_(n-1))
                // trunc_err = dt / (2 * dt + dt_prev) * |x_iter - x_pred| / (tol_time * scale)
                // scale is the field scale over x_n and x_iter
                // the first step has no previous step size, so accept it
                let mut err: f64 = 1.0;
                if dt_prev > 0.0 {
                    let err_fac = dt / (2.0 * dt + dt_prev);
                    let scale_vec = _solu_scale(vars, &x_n1_vec, &x_iter_vec);
                    err = 0.0;
                    for xid in 0..mat_size {
                        let x_pred = x_n1_vec[xid] + (dt / dt_prev) * (x_n1_vec[xid] - x_n2_vec[xid]);
                        let w = 1e-10 + tol_time * scale_vec[xid];
                        err = err.max(err_fac * (x_iter_vec[xid] - x_pred).abs() / w);
                    }
                }

                // if error is too large, then shrink step size
                // accept anyway at dt_min, since retrying would repeat the same step
                if err > 1.0 && dt > dt_min {
                    is_accept = false;
                    let grow = 0.9 * err.powf(-0.5);
                    let grow = grow.max(grow_min);  // limit to grow_min
                    dt = (dt * grow).max(dt_min);
                    _solu_to_unkn(vars, &x_n1_vec);
                    x_iter_vec = x_n1_vec.clone();
                    println!("Predictor error too large. Reducing step size to {dt:.4e}.");
                    continue;
                }

                // update previous data
                is_accept = true;
                x_n2_vec = x_n1_vec;
                x_n1_vec = x_iter_vec.clone();
                dt_prev = dt;
                t += dt;

                // update step size
                let grow = 0.9 * err.max(1e-10).powf(-0.35) * err_prev.max(1e-10).powf(0.2);
                let grow = grow.min(grow_max);  // limit to grow_max
                dt = (dt * grow).clamp(dt_min, dt_max);
                err_prev = err.max(1e-10);

                // write variables
                if t_try == t_write {
                    let time_w0 = Instant::now();
                    write_vars_vtu(vars, num_write)?;
                    let time_w1 = Instant::now();
                    time_wrt += time_w1 - time_w0;
                    num_write += 1;
                    if let Some(ti) = t_out.iter().copied().find(|&ti| ti > t_try) {
                        t_write = ti;
                    }
                }

            }

            // update for next step
            _update_prev(vars);

        }

        // write list of times
        _write_time_csv(&vars.folder, &t_out)?;

        // output time measurement
        let time_end = Instant::now();
        let time_tot = time_end - time_start;
        println!("Solution completed!");
        println!("Total: {:.6} s", time_tot.as_secs_f64());
        println!("  Matrix Assembly : {:.6} s", time_asm.as_secs_f64());
        println!("  Linear Solve    : {:.6} s", time_lin.as_secs_f64());
        println!("  File Output     : {:.6} s", time_wrt.as_secs_f64());

        Ok(())
    }
}

impl<T: TransientBase> TransientAdaptive for T {}

fn _unkn_to_solu(vars: &Variables, x_vec: &mut Vec<f64>) {
    // iterate through unknown scalars
    for scldom in &vars.scldom {
        match scldom.scldom_type {
            SclDomType::Unknown { start } => {
                // iterate through points
                for (vid, &value) in scldom.value.iter().enumerate() {
                    let xid = start + vid;
                    x_vec[xid] = value;
                }
            }
            _ => {continue;}  // skip if not unknown
        }
    }

    // iterate through unknown vectors
    for vecdom in &vars.vecdom {
        match vecdom.vecdom_type {
            VecDomType::Unknown { start } => {
                // iterate through points
                let num_tot = vecdom.value.len();
                for (vid, &value) in vecdom.value.iter().enumerate() {
                    let xid_x = start + vid;
                    let xid_y = start + num_tot + vid;
                    let (value_x, value_y) = value;
                    x_vec[xid_x] = value_x;
                    x_vec[xid_y] = value_y;
                }
            }
            _ => {continue;}  // skip if not unknown
        }
    }

}

fn _solu_scale(vars: &Variables, x_a_vec: &Vec<f64>, x_b_vec: &Vec<f64>) -> Vec<f64> {
    // field scale of every unknown: max magnitude of its field over both solution vectors
    // vector unknowns share one scale, the max of |(v_x, v_y)|
    let mut scale_vec: Vec<f64> = vec![0.0; x_a_vec.len()];

    for scldom in &vars.scldom {
        match scldom.scldom_type {
            SclDomType::Unknown { start } => {
                let num_tot = scldom.value.len();
                let mut scale: f64 = 0.0;
                for x_vec in [x_a_vec, x_b_vec] {
                    for vid in 0..num_tot {
                        scale = scale.max(x_vec[start + vid].abs());
                    }
                }
                for vid in 0..num_tot {
                    scale_vec[start + vid] = scale;
                }
            }
            _ => { continue; }
        }
    }

    for vecdom in &vars.vecdom {
        match vecdom.vecdom_type {
            VecDomType::Unknown { start } => {
                let num_tot = vecdom.value.len();
                let mut scale: f64 = 0.0;
                for x_vec in [x_a_vec, x_b_vec] {
                    for vid in 0..num_tot {
                        scale = scale.max(x_vec[start + vid].hypot(x_vec[start + num_tot + vid]));
                    }
                }
                for vid in 0..2 * num_tot {
                    scale_vec[start + vid] = scale;
                }
            }
            _ => { continue; }
        }
    }

    scale_vec
}

fn _solu_res(vars: &Variables, a_triplet: &Vec<(usize, usize, f64)>, b_vec: &Vec<f64>, x_vec: &Vec<f64>) -> f64 {
    // computes the normalized residual of each field
    // res_fld = |b - A * x|_L1 / (|A * x - A * x_ave|_L1 + |b - A * x_ave|_L1);
    // res = max res_fld
    //
    // x_ave is the average of each unknown field (scalar or vector components)
    // overall residual is the maximum res_fld

    // get index ranges of fields
    let mut xid_range = Vec::new();
    for scldom in vars.scldom.iter() {
        match scldom.scldom_type {
            SclDomType::Unknown { start } => {
                let num_vert = scldom.value.len();
                xid_range.push((start, start + num_vert));
            }
            _ => { continue; }  // skip if not unknown
        }
    }
    for vecdom in vars.vecdom.iter() {
        match vecdom.vecdom_type {
            VecDomType::Unknown { start } => {
                let num_vert = vecdom.value.len();
                xid_range.push((start, start + num_vert));
                xid_range.push((start + num_vert, start + 2*num_vert));
            }
            _ => { continue; }  // skip if not unknown
        }
    }

    // compute x_ave
    let mat_size = x_vec.len();
    let mut x_ave: Vec<f64> = vec![0.0; mat_size];
    for &(xid_lb, xid_ub) in xid_range.iter() {
        let ave = x_vec[xid_lb..xid_ub].iter().sum::<f64>() / ((xid_ub - xid_lb) as f64);
        for xid in xid_lb..xid_ub {
            x_ave[xid] = ave;
        }
    }

    // compute A * x and A * x_ave
    let mut ax_vec: Vec<f64> = vec![0.0; mat_size];
    let mut ax_ave: Vec<f64> = vec![0.0; mat_size];
    for &(r, c, v) in a_triplet.iter() {
        ax_vec[r] += v * x_vec[c];
        ax_ave[r] += v * x_ave[c];
    }

    // compute residual of each field
    let mut res: f64 = 0.0;
    for &(xid_lb, xid_ub) in xid_range.iter() {
        let mut sum_res = 0.0;
        let mut sum_norm = 0.0;
        for xid in xid_lb..xid_ub {
            sum_res += (b_vec[xid] - ax_vec[xid]).abs();
            sum_norm += (ax_vec[xid] - ax_ave[xid]).abs() + (b_vec[xid] - ax_ave[xid]).abs();
        }
        let res_fld = sum_res / (sum_norm + 1e-10);
        if res_fld > res || res_fld.is_nan() {
            res = res_fld;  // keep NaN so that a diverged field fails the tolerance check
        }
    }

    res
}

fn _solu_to_unkn(vars: &mut Variables, x_vec: &Vec<f64>) {
    // iterate through unknown scalars
    for scldom in &mut vars.scldom {
        match scldom.scldom_type {
            SclDomType::Unknown { start } => {
                // iterate through points
                let num_tot = scldom.value.len();
                for vid in 0..num_tot {
                    let xid = start + vid;
                    let value = x_vec[xid];
                    scldom.value[vid] = value;
                }
            }
            _ => {continue;}  // skip if not unknown
        }
    }

    // iterate through unknown vectors
    for vecdom in &mut vars.vecdom {
        match vecdom.vecdom_type {
            VecDomType::Unknown { start } => {
                // iterate through points
                let num_tot = vecdom.value.len();
                for vid in 0..num_tot {
                    let xid_x = start + vid;
                    let xid_y = start + num_tot + vid;
                    let value_x = x_vec[xid_x];
                    let value_y = x_vec[xid_y];
                    vecdom.value[vid] = (value_x, value_y);
                }
            }
            _ => {continue;}  // skip if not unknown
        }
    }

}

fn _update_func(vars: &mut Variables, t: f64) {
    for scldom_id in 0..vars.scldom.len() {
        ScalarDomain::update_func(vars, scldom_id, t);
    }
    for sclbnd_id in 0..vars.sclbnd.len() {
        ScalarBoundary::update_func(vars, sclbnd_id, t);
    }
    for sclitf_id in 0..vars.sclitf.len() {
        ScalarInterface::update_func(vars, sclitf_id, t);
    }
    for vecdom_id in 0..vars.vecdom.len() {
        VectorDomain::update_func(vars, vecdom_id, t);
    }
    for vecbnd_id in 0..vars.vecbnd.len() {
        VectorBoundary::update_func(vars, vecbnd_id, t);
    }
    for vecitf_id in 0..vars.vecitf.len() {
        VectorInterface::update_func(vars, vecitf_id, t);
    }
}

fn _update_grad(vars: &mut Variables) {
    for scldom in vars.scldom.iter_mut() {
        scldom.update_grad(&vars.tri[scldom.dom_id]);
    }
    for vecdom in vars.vecdom.iter_mut() {
        vecdom.update_grad(&vars.tri[vecdom.dom_id]);
    }
}

fn _update_prev(vars: &mut Variables) {
    for scldom in &mut vars.scldom {
        scldom.update_prev();
    }
    for sclbnd in &mut vars.sclbnd {
        sclbnd.update_prev();
    }
    for sclitf in &mut vars.sclitf {
        sclitf.update_prev();
    }
    for vecdom in &mut vars.vecdom {
        vecdom.update_prev();
    }
    for vecbnd in &mut vars.vecbnd {
        vecbnd.update_prev();
    }
    for vecitf in &mut vars.vecitf {
        vecitf.update_prev();
    }
}

fn _write_time_csv(folder: &str, t_out: &[f64]) -> Result<(), VFVChemError> {
    // initialize writer
    let path = Path::new(folder).join("time.csv");
    let file = File::create(&path).map_err(|e| _csv_err(&path, e))?;
    let mut w = BufWriter::new(file);

    // write header
    writeln!(w, "step,time").map_err(|e| _csv_err(&path, e))?;

    // write remaining times
    for (step, &time) in t_out.iter().enumerate() {
        writeln!(w, "{step},{time}").map_err(|e| _csv_err(&path, e))?;
    }

    // flush
    w.flush().map_err(|e| _csv_err(&path, e))?;

    Ok(())
}

fn _csv_err(path: &PathBuf, source: std::io::Error) -> VFVChemError {
    VFVChemError::WriteCsvError { path: path.to_path_buf(), source }
}
