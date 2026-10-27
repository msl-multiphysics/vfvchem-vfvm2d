use crate::base::prelude::*;
use crate::linear::prelude::*;
use rslab::{CscLayout, LuSolver, LuSymbolic, OrderingMethod, RslabError, SolverSettings};
use std::cell::RefCell;

pub struct RslabLu {
    settings: SolverSettings,
    layout: Option<CscLayout>,
    symbolic: Option<LuSymbolic>,
    lu: RefCell<Option<LuSolver<f64>>>,
}

impl RslabLu {
    // sparse LU with AMD ordering
    // fails on a zero pivot
    pub fn new(num_thread: usize) -> Result<Box<dyn LinearBase>, VFVChemError> {
        // error handling
        if num_thread == 0 {
            return Err(VFVChemError::RslabLuNonPositiveThreads { num_thread });
        }

        // exact factorization; num_thread = 0 would use all cores in rslab
        let settings = SolverSettings::exact().with_ordering(OrderingMethod::Amd).with_threads(num_thread);
        let lin = RslabLu { settings, layout: None, symbolic: None, lu: RefCell::new(None) };

        // return
        Ok(Box::new(lin))
    }
}

impl LinearBase for RslabLu {
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

    fn solve(&self, a_triplet: &Vec<(usize, usize, f64)>, b_vec: &Vec<f64>, _x_init: &Vec<f64>) -> Result<Vec<f64>, VFVChemError> {
        // error handling
        let (layout, symbolic) = match (&self.layout, &self.symbolic) {
            (Some(layout), Some(symbolic)) => (layout, symbolic),
            _ => return Err(VFVChemError::RslabLuMissingSparse {}),
        };

        // assemble matrix using the cached layout
        let a_vals: Vec<f64> = a_triplet.iter().map(|&(_, _, val)| val).collect();
        let a_mat = layout.assemble(&a_vals);

        // numeric factorization using the cached symbolic LU
        // refactor reuses the factor buffers of the previous solve
        let lu = match self.lu.take() {
            Some(mut lu) => {
                symbolic.refactor(&a_mat, &self.settings, &mut lu).map_err(|e| _lu_err(e))?;
                lu
            }
            None => symbolic.factor(&a_mat, &self.settings).map_err(|e| _lu_err(e))?,
        };
        let mut lu_cell = self.lu.borrow_mut();
        let lu = lu_cell.insert(lu);

        // matrix solve
        let x_vec = lu.solve(b_vec).map_err(|e| _solve_err(e))?;

        Ok(x_vec)
    }
}

fn _layout_err(source: RslabError) -> VFVChemError {
    VFVChemError::RslabLuLayoutError { message: format!("{source}") }
}

fn _symbolic_err(source: RslabError) -> VFVChemError {
    VFVChemError::RslabLuSymbolicError { message: format!("{source}") }
}

fn _lu_err(source: RslabError) -> VFVChemError {
    VFVChemError::RslabLuLuError { message: format!("{source}") }
}

fn _solve_err(source: RslabError) -> VFVChemError {
    VFVChemError::RslabLuSolveError { message: format!("{source}") }
}
