use crate::base::prelude::*;
use crate::linear::prelude::*;
use crate::transient::transient_base::TransientBase;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::{Duration, Instant};

// Adaptive time stepping for every transient solver (blanket implementation over TransientBase).
//
// Integrator: backward Euler (BDF1), the same time operators as the fixed-step solve. Every time operator and the
// transient continuity term take dt per step and use the previous-step state (value_prev), so the discrete
// conservation identities of the conservative form hold for any sequence of step sizes.
//
// Error estimate (no extra solves): linear extrapolation predictor through the last two accepted states
//   y_p = y_n + (dt / dt_prev) * (y_n - y_{n-1})
// Taylor expansion gives y_BE - y_p = dt * (2 dt + dt_prev) / 2 * y'' and the BE local error is dt^2 / 2 * y'', so
//   LTE ~ dt / (2 dt + dt_prev) * |y_BE - y_p|
// (on the first step, the constant predictor y_p = y_n is used with factor 1/2: a safe overestimate).
// The error norm is over all unknowns
//   max_i |LTE_i| / (atol + rtol * s_i)
// with s_i the field scale of unknown i: max magnitude of its field over y_n and y_BE (|(v_x, v_y)| for vectors),
// so fields of different magnitude (T, c, u, p) are each measured relative to their own size.
// Pressure is a gauge pressure; an absolute pressure (large constant offset) would loosen its error control.
//
// Controller: PI (Gustafsson / Hairer), local order k = 2:
//   fac = safety * err^(-0.7 / k) * err_prev^(0.4 / k), clamped to [fac_min, fac_max], then dt in [dt_min, dt_max]
// a step with err > 1 is rejected and retried with fac = max(fac_min, safety * err^(-1 / k)); a step whose Picard
// iteration does not converge is rejected and retried with dt * fac_fail. No growth is allowed right after a rejection.
// Steps are shortened to land exactly on the requested output times.
//
// Time-dependent functions are evaluated at the end of each step, so a forcing discontinuity (e.g., a source switched
// on at t0) must be listed in out_times; otherwise a pulse that falls inside one step is missed.

//
// Arguments:
//   t_end, dt_init, dt_min, dt_max - end time, initial step, and step size bounds
//   rtol, atol - relative and absolute local error tolerances
//   safety - step size safety factor (0.9)
//   fac_min - largest shrink per step (0.2); fac_max - largest growth per step (2.0)
//   fac_fail - shrink after a Picard failure (0.25)
//   max_reject - consecutive rejections before giving up (30)
//   out_times - output written exactly at these times (plus t = 0 and t_end)
//   max_iter, tol, damp - Picard iteration
//   fixed - if true, dt = dt_init throughout (still hitting the output times); no error control
//   folder_path - folder for adaptive_times.csv and adaptive_stats.csv

pub trait TransientAdaptive: TransientBase {
    // output at num_step + 1 evenly spaced times from 0 to t_end; rtol = 1e-3, atol = 1e-12, safety = 0.9, max_reject = 30
    fn solve_adaptive(
        &mut self, vars: &mut Variables, lin: &mut Box<dyn LinearBase>,
        max_iter: usize, tol: f64, damp: f64,
        t_end: f64, dt_init: f64, dt_min: f64, dt_max: f64, num_step: usize,
        fac_min: f64, fac_max: f64, fac_fail: f64, folder_path: &str
    ) -> Result<(), VFVChemError> {
        if num_step == 0 { return Err(VFVChemError::SolverNonPositiveNumTimeSteps { num_ts: num_step }); }
        let out_times: Vec<f64> = (0..=num_step).map(|i| t_end * (i as f64) / (num_step as f64)).collect();
        self.solve_adaptive_with_options(
            vars, lin, t_end, dt_init, dt_min, dt_max, 1e-3, 1e-12, 0.9, fac_min, fac_max, fac_fail, 30,
            &out_times, max_iter, tol, damp, false, folder_path
        )
    }

    fn solve_adaptive_with_options(
        &mut self, vars: &mut Variables, lin: &mut Box<dyn LinearBase>,
        t_end: f64, dt_init: f64, dt_min: f64, dt_max: f64, rtol: f64, atol: f64, safety: f64,
        fac_min: f64, fac_max: f64, fac_fail: f64, max_reject: usize, out_times: &[f64],
        max_iter: usize, tol: f64, damp: f64, fixed: bool, folder_path: &str
    ) -> Result<(), VFVChemError> {
        self.solve_adaptive_with(
            vars, lin, t_end, dt_init, dt_min, dt_max, rtol, atol, safety, fac_min, fac_max, fac_fail, max_reject,
            out_times, max_iter, tol, damp, fixed, folder_path, &mut |_v: &Variables, _t: f64, _dt: f64| {}
        )
    }

    // hook(vars, t, dt) is called after each accepted step, before the state is copied to value_prev,
    // so value_prev still holds the state at t - dt (useful for balances)
    fn solve_adaptive_with(
        &mut self, vars: &mut Variables, lin: &mut Box<dyn LinearBase>,
        t_end: f64, dt_init: f64, dt_min: f64, dt_max: f64, rtol: f64, atol: f64, safety: f64,
        fac_min: f64, fac_max: f64, fac_fail: f64, max_reject: usize, out_times: &[f64],
        max_iter: usize, tol: f64, damp: f64, fixed: bool, folder_path: &str,
        hook: &mut dyn FnMut(&Variables, f64, f64)
    ) -> Result<(), VFVChemError> {
        // error handling
        if !(t_end > 0.0) { return Err(VFVChemError::SolverNonPositiveTimeStep { dt: t_end }); }
        if !(dt_init > 0.0) || !(dt_min > 0.0) || !(dt_max >= dt_min) { return Err(VFVChemError::SolverNonPositiveTimeStep { dt: dt_init }); }
        if max_iter == 0 { return Err(VFVChemError::SolverNonPositiveMaxIter { max_iter }); }
        if !(tol > 0.0) { return Err(VFVChemError::SolverNonPositiveTolerance { tol }); }
        if !(damp > 0.0 && damp <= 1.0) { return Err(VFVChemError::SolverInvalidDamping { damp }); }
        if !fixed && !(rtol > 0.0 || atol > 0.0) { return Err(VFVChemError::SolverNonPositiveTolerance { tol: rtol }); }

        // start solver
        let time_start = Instant::now();
        println!("Starting adaptive transient solver.");

        // initialize time measurement
        let mut time_dual = Duration::ZERO;
        let mut time_asm = Duration::ZERO;
        let mut time_lin = Duration::ZERO;
        let mut time_wrt = Duration::ZERO;

        // output times; sorted, without duplicates, and ending at t_end
        let t_eps = 1e-12 * t_end;
        let mut outs: Vec<f64> = out_times.iter().cloned().filter(|&t| t > t_eps && t <= t_end).collect();
        outs.sort_by(|a, b| a.total_cmp(b));
        outs.dedup_by(|a, b| (*a - *b).abs() <= t_eps);
        match outs.last_mut() {
            Some(t) if *t >= t_end - t_eps => { *t = t_end; }
            _ => { outs.push(t_end); }
        }

        // file index and time of every output
        if !folder_path.is_empty() {
            create_dir_all(folder_path).map_err(|e| _adaptive_err(folder_path, e))?;
        }
        let times_path = Path::new(folder_path).join("adaptive_times.csv").display().to_string();
        let times_file = File::create(&times_path).map_err(|e| _adaptive_err(&times_path, e))?;
        let mut times_w = BufWriter::new(times_file);
        writeln!(times_w, "file_index,time").map_err(|e| _adaptive_err(&times_path, e))?;

        // initialize dual control volumes and operators
        let time_0 = Instant::now();
        Self::_update_tri(vars)?;
        let time_1 = Instant::now();
        time_dual += time_1 - time_0;
        let mat_size = self.initial_start(vars)?;
        self.initial_operator(vars)?;

        // write initial condition at t = 0; it is the previous state of the first step
        let mut x_n: Vec<f64> = vec![0.0; mat_size];  // last accepted state
        Self::_update_func(vars, 0.0);
        self.initial_dirichlet(vars);  // load Dirichlet boundary values at t = 0 into the initial unknowns
        Self::_unkn_to_solu(vars, &mut x_n);
        Self::_update_func(vars, 0.0);
        Self::_update_grad(vars);
        Self::_update_prev(vars);
        let time_w0 = Instant::now();
        Self::_write(vars, 0)?;
        writeln!(times_w, "0,{:.9e}", 0.0).map_err(|e| _adaptive_err(&times_path, e))?;
        times_w.flush().map_err(|e| _adaptive_err(&times_path, e))?;
        let time_w1 = Instant::now();
        time_wrt += time_w1 - time_w0;

        // statistics
        let mut num_accept = 0;
        let mut num_reject_err = 0;
        let mut num_reject_picard = 0;
        let mut num_picard = 0;
        let mut steps: Vec<(f64, f64, f64)> = Vec::new();  // accepted steps: (t at end of step, dt, error norm)

        // iterate over time steps
        let mut out_idx = 0;
        let mut n_write = 0;
        let mut t = 0.0;
        let mut dt = dt_init.clamp(dt_min, dt_max);
        let mut x_nm1: Option<(Vec<f64>, f64)> = None;  // previous accepted state and its step size
        let mut err_prev: f64 = 1.0;
        let mut n_reject = 0;
        let mut no_grow = false;
        let mut a_triplet: Vec<(usize, usize, f64)> = Vec::new();
        while t < t_end - t_eps {
            // index of the step being attempted; retries after a rejection keep the same index
            let ts = num_accept;

            // land exactly on the next output time; avoid a tiny remainder step
            let t_out = outs[out_idx];
            let mut dt_step = dt;
            let mut hits = false;
            if t + dt_step >= t_out - t_eps {
                dt_step = t_out - t;
                hits = true;
            } else if t + 1.5 * dt_step > t_out {
                dt_step = 0.5 * (t_out - t);
            }
            let t_new = if hits { t_out } else { t + dt_step };

            // initial assembly of A and b at the end of this step
            let time_t0 = Instant::now();
            let mut x_iter = x_n.clone();
            a_triplet.clear();
            let mut b_vec: Vec<f64> = vec![0.0; mat_size];
            Self::_update_func(vars, t_new);
            Self::_update_grad(vars);
            self.update_flow(vars);
            self.assemble_matrix(vars, &mut a_triplet, &mut b_vec, dt_step)?;
            lin.update_sparse(&a_triplet, mat_size)?;
            let time_t1 = Instant::now();
            time_asm += time_t1 - time_t0;

            // iterate to convergence; a failed linear solve or a non-finite residual rejects the step
            let mut is_conv = false;
            let mut res = f64::NAN;
            for iter in 0..max_iter {
                let time_i0 = Instant::now();
                let x_solve = lin.solve(&a_triplet, &b_vec, &x_iter);
                let time_i1 = Instant::now();
                time_lin += time_i1 - time_i0;
                let x_udmp = match x_solve {
                    Ok(x) => x,
                    Err(e) => { println!("Timestep: {ts}; Linear solve failed: {e}"); break; }
                };
                let x_damp: Vec<f64> = (0..mat_size).map(|i| (1.0 - damp) * x_iter[i] + damp * x_udmp[i]).collect();
                a_triplet.clear();
                b_vec = vec![0.0; mat_size];
                Self::_solu_to_unkn(vars, &x_damp);
                Self::_update_func(vars, t_new);
                Self::_update_grad(vars);
                self.update_flow(vars);
                self.assemble_matrix(vars, &mut a_triplet, &mut b_vec, dt_step)?;
                let time_i2 = Instant::now();
                time_asm += time_i2 - time_i1;
                // compute tolerance: rms of the change relative to the field scale of each unknown (1e-12 absolute floor)
                let s_vec = Self::_solu_scale(vars, &x_iter, &x_damp);
                let mut sum = 0.0;
                for i in 0..mat_size {
                    let w = 1e-12 + s_vec[i];
                    sum += ((x_damp[i] - x_iter[i]) / w).powi(2);
                }
                res = (sum / mat_size as f64).sqrt();
                num_picard += 1;
                println!("Timestep: {ts}; Time: {t_new:.4e}; dt: {dt_step:.4e}; Iteration: {iter}; Residual: {res:.4e}");
                x_iter = x_damp;
                if !res.is_finite() { break; }
                if res < tol { is_conv = true; break; }
            }

            // reject on Picard failure; unknowns go back to the last accepted state
            if !is_conv {
                num_reject_picard += 1;
                n_reject += 1;
                if fixed || n_reject > max_reject || dt_step <= dt_min * (1.0 + 1e-12) {
                    return Err(VFVChemError::TransientDidNotConverge { ts: num_accept, max_iter, res });
                }
                Self::_solu_to_unkn(vars, &x_n);
                dt = (dt_step * fac_fail).max(dt_min);
                no_grow = true;
                println!("Timestep: {ts}; Picard did not converge; retry with dt = {dt:.4e}");
                continue;
            }

            // error estimate
            let err = if fixed { 0.0 } else {
                // max_i |LTE_i| / (atol + rtol * s_i), s_i the field scale over y_n and y_BE
                let s_vec = Self::_solu_scale(vars, &x_n, &x_iter);
                let mut e: f64 = 0.0;
                for i in 0..mat_size {
                    let (yp, fac) = match &x_nm1 {
                        Some((xm, dtm)) => (x_n[i] + dt_step / dtm * (x_n[i] - xm[i]), dt_step / (2.0 * dt_step + dtm)),
                        None => (x_n[i], 0.5),
                    };
                    let w = atol + rtol * s_vec[i];
                    e = e.max(fac * (x_iter[i] - yp).abs() / w);
                }
                e
            };

            // reject on error; at dt_min the step is accepted anyway
            if err > 1.0 && dt_step > dt_min * (1.0 + 1e-12) {
                num_reject_err += 1;
                n_reject += 1;
                if n_reject > max_reject { return Err(VFVChemError::TransientDidNotConverge { ts: num_accept, max_iter, res: err }); }
                Self::_solu_to_unkn(vars, &x_n);
                dt = (dt_step * (safety * err.powf(-0.5)).max(fac_min)).max(dt_min);
                no_grow = true;
                println!("Timestep: {ts}; Error: {err:.4e}; rejected; retry with dt = {dt:.4e}");
                continue;
            }

            // accept
            hook(vars, t_new, dt_step);
            Self::_update_prev(vars);
            num_accept += 1;
            println!("Timestep: {ts}; Time: {t_new:.4e}; dt: {dt_step:.4e}; Error: {err:.4e}; accepted");
            steps.push((t_new, dt_step, err));
            x_nm1 = Some((x_n, dt_step));
            x_n = x_iter;
            t = t_new;
            n_reject = 0;

            // write the state at an output time as the next snapshot
            if hits {
                n_write += 1;
                let time_w2 = Instant::now();
                Self::_write(vars, n_write)?;
                writeln!(times_w, "{n_write},{t:.9e}").map_err(|e| _adaptive_err(&times_path, e))?;
                times_w.flush().map_err(|e| _adaptive_err(&times_path, e))?;
                let time_w3 = Instant::now();
                time_wrt += time_w3 - time_w2;
                out_idx += 1;
                if out_idx >= outs.len() { break; }
            }

            // next step size; a step shortened to hit an output time does not limit the next one
            if fixed { dt = dt_init; continue; }
            let e = err.max(1e-10);
            let mut fac = safety * e.powf(-0.35) * err_prev.max(1e-10).powf(0.2);
            fac = fac.clamp(fac_min, fac_max);
            if no_grow { fac = fac.min(1.0); no_grow = false; }
            let base = if hits { dt.max(dt_step) } else { dt_step };
            dt = (base * fac).clamp(dt_min, dt_max);
            err_prev = e;
        }

        // write statistics; summary counts are repeated on every accepted step
        let time_w4 = Instant::now();
        let stats_path = Path::new(folder_path).join("adaptive_stats.csv").display().to_string();
        let stats_file = File::create(&stats_path).map_err(|e| _adaptive_err(&stats_path, e))?;
        let mut stats_w = BufWriter::new(stats_file);
        writeln!(stats_w, "accepted,rejected_error,rejected_picard,picard_iterations,time,dt,error").map_err(|e| _adaptive_err(&stats_path, e))?;
        for &(t_step, dt_step, err) in &steps {
            writeln!(stats_w, "{num_accept},{num_reject_err},{num_reject_picard},{num_picard},{t_step:.9e},{dt_step:.9e},{err:.9e}").map_err(|e| _adaptive_err(&stats_path, e))?;
        }
        stats_w.flush().map_err(|e| _adaptive_err(&stats_path, e))?;
        let time_w5 = Instant::now();
        time_wrt += time_w5 - time_w4;

        // output time measurement and statistics
        let time_end = Instant::now();
        let time_total = time_end - time_start;
        println!("Adaptive solution completed!");
        println!("Total: {:.6} s", time_total.as_secs_f64());
        println!("  Dual Mesh         : {:.6} s", time_dual.as_secs_f64());
        println!("  Matrix Assembly   : {:.6} s", time_asm.as_secs_f64());
        println!("  Matrix Solve      : {:.6} s", time_lin.as_secs_f64());
        println!("  File Output       : {:.6} s", time_wrt.as_secs_f64());
        println!("  Accepted Steps    : {num_accept}");
        println!("  Rejected (Error)  : {num_reject_err}");
        println!("  Rejected (Picard) : {num_reject_picard}");
        println!("  Picard Iterations : {num_picard}");

        Ok(())
    }
}

impl<T: TransientBase> TransientAdaptive for T {}

fn _adaptive_err(file_path: &str, source: std::io::Error) -> VFVChemError {
    VFVChemError::WriteCsv { file_path: file_path.to_string(), source }
}
