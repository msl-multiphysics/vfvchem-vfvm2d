use crate::base::prelude::*;
use crate::linear::prelude::*;
use crate::transient::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

const GROWTH_MIN: f64 = 0.25;
const GROWTH_MAX: f64 = 2.00;  // below 1 + sqrt(2) for zero-stability of variable-step bdf2
const GROWTH_FAIL: f64 = 0.50;
const TOL_TIME: f64 = 100.0;
const SAFETY: f64 = 0.9;
const MAX_REJECT: usize = 30;  // consecutive rejections before giving up

pub trait TransientAdaptive: TransientBase {
    // adaptive solver with default options
    fn solve_adaptive(
        &mut self, vars: &mut Variables, lin: &mut Box<dyn LinearBase>,
        max_iter: usize, tol: f64, damp: f64, dt_init: f64, dt_min: f64, dt_max: f64,
        t_end: f64, num_step: usize,
    ) -> Result<(), VFVChemError> {
        // error handling
        if num_step == 0 { return Err(VFVChemError::SolverNonPositiveNumTimeSteps { num_ts: num_step }); }

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
    //
    // integration
    // the first step compares one backward euler step of dt with two backward euler steps of dt / 2
    // the two half steps are kept; the initial state is not used for prediction, since algebraic unknowns
    // (e.g., pressure) and constrained fields (e.g., velocity) may be inconsistent at t = 0
    // variable-step bdf2 afterward; w = dt / dt_prv1
    //   du/dt = ((1 + 2w) / (1 + w) * u - (1 + w) * u_prv1 + w^2 / (1 + w) * u_prv2) / dt
    //
    // local error estimate; no extra solves after the first step
    // first step: err = |x - x_full|, the difference between the half steps and the full step
    // second step: linear predictor through the two half steps; conservative backward euler factor
    //   x_pred = x_n + (dt / dt_prv1) * (x_n - x_(n-1))
    //   err = dt / (2 * dt + dt_prv1) * |x - x_pred|
    // later steps: quadratic predictor through the last three accepted states
    //   x - x_pred = (c_cor + c_prd) * x'''
    //   c_cor = dt^2 * (dt + dt_prv1)^2 / (6 * (2 * dt + dt_prv1))  (bdf2)
    //   c_prd = dt * (dt + dt_prv1) * (dt + dt_prv1 + dt_prv2) / 6  (predictor)
    //   err = c_cor / (c_cor + c_prd) * |x - x_pred|
    // err is divided by tol_time times the field scale (max magnitude of the field over x_n and x)
    //
    // step control; err ~ dt^ord
    // rejected: dt *= max(0.9 * err^(-1 / ord), grow_min), or dt *= grow_fail if not converged
    // accepted: dt *= 0.9 * err^(-0.7 / ord) * err_prv^(0.4 / ord), limited to [grow_min, grow_max]
    // no growth right after a rejection; grow_max < 1 + sqrt(2) keeps variable-step bdf2 zero-stable
    // steps are shortened to land on the output times; time-dependent functions are evaluated at the end of each step
    fn solve_adaptive_options(
        &mut self, vars: &mut Variables, lin: &mut Box<dyn LinearBase>,
        max_iter: usize, tol_iter: f64, damp: f64, tol_time: f64, dt_init: f64, dt_min: f64, dt_max: f64,
        t_out: Vec<f64>, grow_min: f64, grow_max: f64, grow_fail: f64
    ) -> Result<(), VFVChemError> {
        // error handling
        if max_iter == 0 { return Err(VFVChemError::SolverNonPositiveMaxIter { max_iter }); }
        if !(tol_iter > 0.0) { return Err(VFVChemError::SolverNonPositiveTolerance { tol: tol_iter }); }
        if !(tol_time > 0.0) { return Err(VFVChemError::SolverNonPositiveTolerance { tol: tol_time }); }
        if !(damp > 0.0 && damp <= 1.0) { return Err(VFVChemError::SolverInvalidDamping { damp }); }
        if !(dt_init > 0.0) { return Err(VFVChemError::SolverNonPositiveInitialStepSize { dt_init }) };
        if !(dt_min > 0.0 && dt_min <= dt_init) { return Err(VFVChemError::SolverInvalidMinimumStepSize { dt_min }) };
        if !(dt_max >= dt_init) { return Err(VFVChemError::SolverInvalidMaximumStepSize { dt_max }) };
        if !(grow_min > 0.0 && grow_min <= 1.0) { return Err(VFVChemError::SolverInvalidGrowthMinimum { grow_min }) };
        if !(grow_max > 1.0 && grow_max < 1.0 + 2.0_f64.sqrt()) { return Err(VFVChemError::SolverInvalidGrowthMaximum { grow_max }) };
        if !(grow_fail > 0.0 && grow_fail < 1.0) { return Err(VFVChemError::SolverInvalidGrowthFail { grow_fail }) };

        // error handling of t_out
        if t_out.is_empty() || !(t_out[t_out.len() - 1] > 0.0) { return Err(VFVChemError::SolverAdaptiveInvalidOutputTime {}); }  // positive final time
        for i in 0..t_out.len() - 1 {
            if !(t_out[i] >= 0.0) || (t_out[i] >= t_out[i + 1]) {  // t_out non-negative and strictly increasing
                return Err(VFVChemError::SolverAdaptiveInvalidOutputTime {});
            }
        }

        // initialize time measurement
        println!("Starting adaptive transient solver.");
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
        let mut x_n1_vec: Vec<f64> = vec![0.0; mat_size];  // last accepted state x_n

        // initialize operators
        self.initial_operator(vars)?;

        // load dirichlet values into unknown fields
        _update_func(vars, 0.0);  // update function-type variables
        self.initial_dirichlet(vars, mat_size)?;  // set dirichlet values

        // initial assembly of A and b
        // used for sparsity pattern
        let time_2 = Instant::now();
        _unkn_to_solu(vars, &mut x_n1_vec);  // unknown fields -> solution vector
        _update_func(vars, 0.0);
        _update_grad(vars);
        _update_prev(vars);
        _update_prev(vars);  // twice to fill both previous timesteps
        self.assemble_matrix(vars, &mut a_triplet, &mut b_vec, (1.0 / dt_init, -1.0 / dt_init, 0.0))?;
        lin.set_sparse(&a_triplet, mat_size)?;
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
        let t_end = t_out[t_out.len() - 1];  // final time
        let t_eps = 1e-12 * t_end;  // tolerance for landing on output times
        let mut t = 0.0;  // current time
        let mut t_write = t_out.iter().copied().find(|&ti| ti > 0.0).unwrap();  // next time to write output
        let mut dt = dt_init;  // step size of the next attempt
        let mut dt_prv1 = 0.0;  // step size between x_(n-1) and x_n
        let mut dt_prv2 = 0.0;  // step size between x_(n-2) and x_(n-1)
        let mut x_n2_vec = x_n1_vec.clone();  // x_(n-1)
        let mut x_n3_vec = x_n1_vec.clone();  // x_(n-2)
        let mut err_prv: f64 = 1.0;  // local error of the previous accepted step

        // statistics
        let mut num_acc = 0;  // accepted steps
        let mut num_rej_time = 0;  // rejected steps due to local error
        let mut num_rej_iter = 0;  // rejected steps due to nonlinear convergence
        let mut dt_lo = f64::INFINITY;  // smallest accepted step size
        let mut dt_hi: f64 = 0.0;  // largest accepted step size

        // iterate until final time
        while t < t_end - t_eps {
            // repeat until the step is accepted
            let mut num_rej = 0;  // consecutive rejections
            loop {
                if num_rej == MAX_REJECT { return Err(VFVChemError::SolverAdaptiveTooManyRejections { max_reject: MAX_REJECT }); }

                // adjust dt to land on the next output time; avoid a small remainder step
                let mut is_land = false;
                if t + dt >= t_write - t_eps {
                    dt = t_write - t;
                    is_land = true;
                } else if t + 1.5 * dt > t_write {
                    dt = 0.5 * (t_write - t);
                }
                let t_try = if is_land { t_write } else { t + dt };

                // attempt a solution with step size dt
                // the first step also solves two half steps; x_full_vec is the single backward euler step
                let mut x_full_vec = Vec::new();
                let mut x_half_vec = Vec::new();
                // time weights; du/dt = tw0 * u + tw1 * u_prv1 + tw2 * u_prv2
                let x_conv = 'attempt: {
                    if num_acc > 0 {
                        let w = dt / dt_prv1;
                        let twgt_bdf2 = ((1.0 + 2.0 * w) / ((1.0 + w) * dt), -(1.0 + w) / dt, w * w / ((1.0 + w) * dt));
                        break 'attempt _solve_step(self, vars, lin, &x_n1_vec, t_try, twgt_bdf2, max_iter, tol_iter, damp, &mut time_asm, &mut time_lin)?;
                    }
                    let twgt_full = (1.0 / dt, -1.0 / dt, 0.0);  // backward euler with dt
                    let twgt_half = (2.0 / dt, -2.0 / dt, 0.0);  // backward euler with dt / 2
                    let Some(x_vec) = _solve_step(self, vars, lin, &x_n1_vec, t_try, twgt_full, max_iter, tol_iter, damp, &mut time_asm, &mut time_lin)? else { break 'attempt None };
                    x_full_vec = x_vec;
                    let Some(x_vec) = _solve_step(self, vars, lin, &x_n1_vec, t + 0.5 * dt, twgt_half, max_iter, tol_iter, damp, &mut time_asm, &mut time_lin)? else { break 'attempt None };
                    x_half_vec = x_vec;
                    _update_prev(vars);
                    _solve_step(self, vars, lin, &x_half_vec, t_try, twgt_half, max_iter, tol_iter, damp, &mut time_asm, &mut time_lin)?
                };

                // predictor x_pred and error factor fac; err ~ dt^ord
                let (x_pred_vec, fac, ord) = if num_acc == 0 {
                    (x_full_vec, 1.0, 2.0)
                } else if num_acc == 1 {
                    let r = dt / dt_prv1;
                    let x_pred_vec: Vec<f64> = (0..mat_size).map(|i| x_n1_vec[i] + r * (x_n1_vec[i] - x_n2_vec[i])).collect();
                    (x_pred_vec, dt / (2.0 * dt + dt_prv1), 2.0)
                } else {
                    let (h, h1, h2) = (dt, dt_prv1, dt_prv2);
                    let l0 = (h + h1) * (h + h1 + h2) / (h1 * (h1 + h2));  // lagrange weights at t + dt
                    let l1 = -h * (h + h1 + h2) / (h1 * h2);
                    let l2 = h * (h + h1) / ((h1 + h2) * h2);
                    let x_pred_vec: Vec<f64> = (0..mat_size).map(|i| l0 * x_n1_vec[i] + l1 * x_n2_vec[i] + l2 * x_n3_vec[i]).collect();
                    let c_cor = h * h * (h + h1).powi(2) / (6.0 * (2.0 * h + h1));
                    let c_prd = h * (h + h1) * (h + h1 + h2) / 6.0;
                    (x_pred_vec, c_cor / (c_cor + c_prd), 3.0)
                };

                // local error; infinite if not converged
                let mut err = f64::INFINITY;
                if let Some(x_vec) = &x_conv {
                    let scale_vec = _solu_scale(vars, &x_n1_vec, x_vec);
                    err = 0.0;
                    for xid in 0..mat_size {
                        let w = 1e-10 + tol_time * scale_vec[xid];
                        err = err.max(fac * (x_vec[xid] - x_pred_vec[xid]).abs() / w);
                    }
                }

                // reject if not converged or if error is too large
                if x_conv.is_none() || err > 1.0 {
                    // raise error at dt_min
                    if dt <= dt_min && x_conv.is_none() { return Err(VFVChemError::SolverAdaptiveDidNotConverge { dt_min }); }
                    if dt <= dt_min { return Err(VFVChemError::SolverAdaptiveTimeErrorTooLarge { dt_min }); }

                    // shrink step size
                    if x_conv.is_none() {
                        dt = (dt * grow_fail).max(dt_min);
                        num_rej_iter += 1;
                        println!("Nonlinear solver did not converge. Reducing step size to {dt:.4e}.");
                    } else {
                        dt = (dt * (SAFETY * err.powf(-1.0 / ord)).max(grow_min)).max(dt_min);
                        num_rej_time += 1;
                        println!("Local error {err:.4e} too large. Reducing step size to {dt:.4e}.");
                    }
                    num_rej += 1;

                    // restore the last accepted state; the first step also restores the initial history
                    _solu_to_unkn(vars, &x_n1_vec);
                    if num_acc == 0 {
                        _update_func(vars, t);
                        _update_grad(vars);
                        _update_prev(vars);
                        _update_prev(vars);
                    }
                    continue;
                }

                // update previous data
                // the half steps of the first step become the history
                if num_acc == 0 {
                    x_n2_vec = x_half_vec;
                    dt_prv1 = 0.5 * dt;
                    dt_prv2 = 0.5 * dt;
                } else {
                    x_n3_vec = x_n2_vec;
                    x_n2_vec = x_n1_vec;
                    dt_prv2 = dt_prv1;
                    dt_prv1 = dt;
                }
                x_n1_vec = x_conv.unwrap();
                t = t_try;
                num_acc += 1;
                dt_lo = dt_lo.min(dt_prv1);
                dt_hi = dt_hi.max(dt_prv1);
                println!("Step accepted. Time: {t:.4e}; Step size: {dt_prv1:.4e}; Local error: {err:.4e}");

                // update step size from the last accepted step
                let mut grow = SAFETY * err.max(1e-10).powf(-0.7 / ord) * err_prv.max(1e-10).powf(0.4 / ord);
                grow = grow.clamp(grow_min, grow_max);
                if num_rej > 0 { grow = grow.min(1.0); }
                dt = (dt_prv1 * grow).min(dt_max);
                err_prv = err;

                // write variables
                if is_land {
                    let time_w0 = Instant::now();
                    write_vars_vtu(vars, num_write)?;
                    let time_w1 = Instant::now();
                    time_wrt += time_w1 - time_w0;
                    num_write += 1;
                    if let Some(ti) = t_out.iter().copied().find(|&ti| ti > t) {
                        t_write = ti;
                    }
                }
                break;
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
        println!("Accepted steps: {num_acc}; Rejected steps: {num_rej_time} (local error), {num_rej_iter} (nonlinear)");
        println!("Step size: {dt_lo:.4e} to {dt_hi:.4e}");

        Ok(())
    }
}

impl<T: TransientBase> TransientAdaptive for T {}

fn _solve_step<T: TransientBase + ?Sized>(
    phys: &T, vars: &mut Variables, lin: &mut Box<dyn LinearBase>, x_init_vec: &Vec<f64>, t: f64, twgt: (f64, f64, f64),
    max_iter: usize, tol: f64, damp: f64, time_asm: &mut Duration, time_lin: &mut Duration,
) -> Result<Option<Vec<f64>>, VFVChemError> {
    // solves one step from x_init to time t with time weights twgt
    // returns None if not converged

    // assembly of A and b at time t
    // used for first iteration
    let time_t0 = Instant::now();
    let mat_size = x_init_vec.len();
    let mut x_iter_vec = x_init_vec.clone();  // damped result
    let mut a_triplet: Vec<(usize, usize, f64)> = Vec::new();
    let mut b_vec: Vec<f64> = vec![0.0; mat_size];
    _solu_to_unkn(vars, &x_iter_vec);
    _update_func(vars, t);
    _update_grad(vars);
    phys.assemble_matrix(vars, &mut a_triplet, &mut b_vec, twgt)?;
    let time_t1 = Instant::now();
    *time_asm += time_t1 - time_t0;

    // iterate to convergence
    for iter in 0..max_iter {
        // solve for undamped solution
        // x_udmp in A_k * x_udmp = b_k
        let time_i0 = Instant::now();
        let x_udmp_vec = lin.solve(&a_triplet, &b_vec, &x_iter_vec)?;
        let time_i1 = Instant::now();
        *time_lin += time_i1 - time_i0;

        // solve for damped solution
        // x_(k+1) in x_damp = (1 - damp) * x_k + damp * x_udmp
        let x_damp_vec: Vec<f64> = (0..mat_size).map(|i| (1.0 - damp) * x_iter_vec[i] + damp * x_udmp_vec[i]).collect();

        // assemble A_(k+1) and b_(k+1) using x_(k+1)
        let time_i2 = Instant::now();
        a_triplet.clear();  // reset A
        b_vec.fill(0.0);  // reset b
        _solu_to_unkn(vars, &x_damp_vec);
        _update_func(vars, t);
        _update_grad(vars);
        phys.assemble_matrix(vars, &mut a_triplet, &mut b_vec, twgt)?;
        let time_i3 = Instant::now();
        *time_asm += time_i3 - time_i2;

        // compute tolerance
        let res = _solu_res(vars, &a_triplet, &b_vec, &x_damp_vec);
        println!("Time: {t:.4e}; Iteration: {iter}; Residual: {res:.4e}");

        // reject if diverged
        if !res.is_finite() {
            return Ok(None);
        }

        // update iteration vector
        x_iter_vec = x_damp_vec;

        // tolerance check
        if res < tol {
            return Ok(Some(x_iter_vec));
        }
    }

    Ok(None)
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
    for vecdom in &mut vars.vecdom {
        vecdom.update_prev();
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
