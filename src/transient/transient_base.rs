use crate::base::prelude::*;
use crate::linear::prelude::*;
use std::time::{Duration, Instant};

pub trait TransientBase {
    // to implement in specific solvers
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError>;  // set start indices and compute matrix sizes
    fn initial_dirichlet(&mut self, vars: &mut Variables, mat_size: usize) -> Result<(), VFVChemError>;  // load dirichlet values into unknown-type variables
    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError>;  // initialize operators
    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, twgt: (f64, f64, f64)) -> Result<(), VFVChemError>;  // assemble global matrix

    // default solver loop
    fn solve(
        &mut self, vars: &mut Variables, lin: &mut Box<dyn LinearBase>,
        max_iter: usize, tol: f64, damp: f64, dt: f64, num_ts: usize, num_ts_write: usize,
    ) -> Result<(), VFVChemError> {
        // error handling
        if max_iter == 0 { return Err(VFVChemError::SolverNonPositiveMaxIter { max_iter }); }
        if tol <= 0.0 { return Err(VFVChemError::SolverNonPositiveTolerance { tol }); }
        if damp <= 0.0 || damp > 1.0 { return Err(VFVChemError::SolverInvalidDamping { damp }); }
        if dt <= 0.0 { return Err(VFVChemError::SolverNonPositiveTimeStep { dt }); }
        if num_ts == 0 { return Err(VFVChemError::SolverNonPositiveNumTimeSteps { num_ts }); }
        if num_ts_write == 0 { return Err(VFVChemError::SolverNonPositiveWriteFrequency { num_ts_write }); }

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

        // time weights
        // du/dt = tw0 * u + tw1 * u_prv1 + tw2 * u_prv2
        // bdf1 for the first timestep; bdf2 afterward
        let twgt_bdf1 = (1.0 / dt, -1.0 / dt, 0.0);
        let twgt_bdf2 = (1.5 / dt, -2.0 / dt, 0.5 / dt);

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
        _update_prev(vars);  // twice to fill both previous timesteps
        self.assemble_matrix(vars, &mut a_triplet, &mut b_vec, twgt_bdf1)?;
        lin.set_sparse(&a_triplet, mat_size)?;
        let num_triplet = a_triplet.len();
        let time_3 = Instant::now();
        time_asm += time_3 - time_2;

        // write initial state
        write_vars_vtu(vars, 0)?;
        let time_4 = Instant::now();
        time_wrt += time_4 - time_3;

        // iterate over time steps
        for ts in 0..num_ts {
            // assembly of A and b at timestep ts
            // used for first iteration
            let time_t0 = Instant::now();
            let t = dt * ((ts + 1) as f64);
            let twgt = if ts == 0 { twgt_bdf1 } else { twgt_bdf2 };
            a_triplet = Vec::with_capacity(num_triplet);  // reset A
            b_vec = vec![0.0; mat_size];  // reset b
            _update_func(vars, t);
            _update_grad(vars);
            self.assemble_matrix(vars, &mut a_triplet, &mut b_vec, twgt)?;
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
                _update_func(vars, t);
                _update_grad(vars);
                self.assemble_matrix(vars, &mut a_triplet, &mut b_vec, twgt)?;
                let time_i3 = Instant::now();
                time_asm += time_i3 - time_i2;

                // compute tolerance
                let res = _solu_res(vars, &a_triplet, &b_vec, &x_damp_vec);
                println!("Timestep: {ts}; Iteration: {iter}; Residual: {res:.4e}");

                // error if diverged
                if !res.is_finite() {
                    return Err(VFVChemError::SolverNonFiniteResidual { });
                }

                // update iteration vector
                x_iter_vec = x_damp_vec;

                // tolerance check
                if res < tol {
                    is_conv = true;
                    break;
                }
            }

            // error if not converged
            if !is_conv {
                return Err(VFVChemError::SolverDidNotConverge { max_iter });
            }

            // write variables
            if (ts + 1) % num_ts_write == 0 {
                let time_w0 = Instant::now();
                write_vars_vtu(vars, (ts + 1) / num_ts_write)?;
                let time_w1 = Instant::now();
                time_wrt += time_w1 - time_w0;
            }

            // update for next timestep
            _update_prev(vars);
        }

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
    for vecdom in &mut vars.vecdom {
        vecdom.update_prev();
    }
}
