use crate::base::prelude::*;
use crate::linear::prelude::*;
use rslab::{gmres, CscLayout, KrylovSettings, LuSolver, LuSymbolic, OrderingMethod, RslabError, SolverSettings};
use std::cell::RefCell;

pub struct RslabGmres {
    settings: SolverSettings,
    layout: Option<CscLayout>,
    symbolic: Option<LuSymbolic>,
    lu: RefCell<Option<LuSolver<f64>>>,
    tol: f64,
    restart: usize,
    max_iter: usize,
}

impl RslabGmres {
    // default: incomplete-LU-preconditioned flexible GMRES with AMD ordering
    // restart = Krylov subspace size per cycle (memory ~ O(n * restart))
    // max_iter = total number of GMRES iterations
    // drop_tol = factor entries below drop_tol relative to their column are dropped
    pub fn new(num_thread: usize) -> Result<Box<dyn LinearBase>, VFVChemError> {
        Self::new_with_options(num_thread, 1e-6, 50, 10_000, 1e-2)
    }

    pub fn new_with_options(num_thread: usize, tol: f64, restart: usize, max_iter: usize, drop_tol: f64) -> Result<Box<dyn LinearBase>, VFVChemError> {
        // error handling
        if num_thread == 0 {
            return Err(VFVChemError::RslabGmresNonPositiveThreads { num_thread });
        }
        if tol.is_nan() || tol <= 0.0 {
            return Err(VFVChemError::RslabGmresNonPositiveTolerance { tol });
        }
        if restart == 0 {
            return Err(VFVChemError::RslabGmresNonPositiveRestart { restart });
        }
        if max_iter == 0 {
            return Err(VFVChemError::RslabGmresNonPositiveMaxIter { max_iter });
        }
        if drop_tol.is_nan() || drop_tol <= 0.0 {
            return Err(VFVChemError::RslabGmresNonPositiveDropTolerance { drop_tol });
        }

        // never-fail preconditioner: pivots below 1e-8 are lifted instead of failing
        let settings = SolverSettings::preconditioner(1e-8).with_drop_tol(drop_tol).with_ordering(OrderingMethod::Amd).with_threads(num_thread);
        let lin = RslabGmres { settings, layout: None, symbolic: None, lu: RefCell::new(None), tol, restart, max_iter };

        // return
        Ok(Box::new(lin))
    }
}

impl LinearBase for RslabGmres {
    fn update_sparse(&mut self, a_triplet: &Vec<(usize, usize, f64)>, mat_size: usize) -> Result<(), VFVChemError> {
        // collect rows, columns, and values
        let a_rows: Vec<usize> = a_triplet.iter().map(|&(row, _, _)| row).collect();
        let a_cols: Vec<usize> = a_triplet.iter().map(|&(_, col, _)| col).collect();
        let a_vals: Vec<f64> = a_triplet.iter().map(|&(_, _, val)| val).collect();

        // layout and symbolic LU factorization
        // reused for every later matrix with this sparsity pattern
        let layout = CscLayout::new(mat_size, &a_rows, &a_cols).map_err(|e| _layout_err(e))?;
        let a_mat = layout.assemble(&a_vals);
        let symbolic = LuSymbolic::analyze(&a_mat, &self.settings).map_err(|e| _symbolic_err(e))?;
        self.layout = Some(layout);
        self.symbolic = Some(symbolic);
        self.lu = RefCell::new(None);

        Ok(())
    }

    fn solve(&self, a_triplet: &Vec<(usize, usize, f64)>, b_vec: &Vec<f64>, x_init: &Vec<f64>) -> Result<Vec<f64>, VFVChemError> {
        // error handling
        let (layout, symbolic) = match (&self.layout, &self.symbolic) {
            (Some(layout), Some(symbolic)) => (layout, symbolic),
            _ => return Err(VFVChemError::RslabGmresMissingSparse {}),
        };

        // assemble matrix using the cached layout
        let a_vals: Vec<f64> = a_triplet.iter().map(|&(_, _, val)| val).collect();
        let a_mat = layout.assemble(&a_vals);

        // incomplete LU preconditioner using the cached symbolic LU
        // refactor reuses the factor buffers of the previous solve
        let lu = match self.lu.take() {
            Some(mut lu) => {
                symbolic.refactor(&a_mat, &self.settings, &mut lu).map_err(|e| _lu_err(e))?;
                lu
            }
            None => symbolic.factor(&a_mat, &self.settings).map_err(|e| _lu_err(e))?,
        };
        let mut lu_cell = self.lu.borrow_mut();
        let lu: &LuSolver<f64> = lu_cell.insert(lu);

        // preconditioned flexible GMRES, starting from x_init
        let krylov = KrylovSettings::default().with_tol(self.tol).with_restart(self.restart).with_max_iter(self.max_iter);
        let res = gmres(&a_mat, b_vec.as_slice(), lu, &krylov, Some(x_init.as_slice())).map_err(|e| _solve_err(e))?;
        if !res.converged {
            return Err(VFVChemError::RslabGmresDidNotConverge { iters: res.iters, res: res.final_res });
        }

        Ok(res.x)
    }
}

fn _layout_err(source: RslabError) -> VFVChemError {
    VFVChemError::RslabGmresLayoutError { message: format!("{source}") }
}

fn _symbolic_err(source: RslabError) -> VFVChemError {
    VFVChemError::RslabGmresSymbolicError { message: format!("{source}") }
}

fn _lu_err(source: RslabError) -> VFVChemError {
    VFVChemError::RslabGmresLuError { message: format!("{source}") }
}

fn _solve_err(source: RslabError) -> VFVChemError {
    VFVChemError::RslabGmresSolveError { message: format!("{source}") }
}
