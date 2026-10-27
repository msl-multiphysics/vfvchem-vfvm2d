use crate::base::prelude::*;
use crate::linear::prelude::*;
use faer::linalg::solvers::Solve;
use faer::sparse::{Argsort, Pair, SparseColMat, SymbolicSparseColMat};
use faer::sparse::linalg::solvers::{Lu, SymbolicLu};
use faer::{ColMut, Par};
use std::num::NonZeroUsize;

pub struct FaerLu {
    symbolic: SymbolicSparseColMat<usize, usize, usize>,
    argsort: Argsort<usize>,
    symbolic_lu: SymbolicLu<usize>,
}

impl FaerLu {
    pub fn new(num_thread: usize) -> Result<Box<dyn LinearBase>, VFVChemError> {
        // placeholder default
        let (symbolic, argsort) = SymbolicSparseColMat::try_new_from_indices(0, 0, &[]).map_err(|e| _symbolic_err(e))?;
        let symbolic_lu = SymbolicLu::try_new(symbolic.as_ref()).map_err(|e| _symbolic_lu_err(e))?;
        let lin = FaerLu { symbolic, argsort, symbolic_lu };

        // set number of threads
        let num_thread_nz = NonZeroUsize::new(num_thread).unwrap();
        faer::set_global_parallelism(Par::Rayon(num_thread_nz));

        // return
        Ok(Box::new(lin))
    }
}

impl LinearBase for FaerLu {
    fn update_sparse(&mut self, a_triplet: &Vec<(usize, usize, f64)>, mat_size: usize) -> Result<(), VFVChemError> {
        // collect (row, col) pairs
        let a_pair: Vec<Pair<usize, usize>> = a_triplet.iter().map(|&(row, col, _)| Pair::new(row, col)).collect();

        // symbolic sparse matrix and LU factorization
        // reused for every later matrix with this sparsity pattern
        let (symbolic, argsort) = SymbolicSparseColMat::try_new_from_indices(mat_size, mat_size, &a_pair).map_err(|e| _symbolic_err(e))?;
        let symbolic_lu = SymbolicLu::try_new(symbolic.as_ref()).map_err(|e| _symbolic_lu_err(e))?;
        self.symbolic = symbolic;
        self.argsort = argsort;
        self.symbolic_lu = symbolic_lu;

        Ok(())
    }

    fn solve(&self, a_triplet: &Vec<(usize, usize, f64)>, b_vec: &Vec<f64>, _x_init: &Vec<f64>) -> Result<Vec<f64>, VFVChemError> {
        // collect values
        let a_vals: Vec<f64> = a_triplet.iter().map(|&(_, _, val)| val).collect();

        // assemble matrix using symbolic factorization
        let a_mat = SparseColMat::new_from_argsort(self.symbolic.clone(), &self.argsort, &a_vals).map_err(|e| _matrix_err(e))?;

        // numeric factorization using the cached symbolic LU
        // cloning SymbolicLu only copies its reference count
        let lu = Lu::try_new_with_symbolic(self.symbolic_lu.clone(), a_mat.as_ref()).map_err(|e| _lu_err(e))?;

        // matrix solve
        let mut x_vec = b_vec.to_vec();
        lu.solve_in_place(ColMut::from_slice_mut(&mut x_vec));

        Ok(x_vec)
    }
}

fn _symbolic_err(source: faer::sparse::CreationError) -> VFVChemError {
    VFVChemError::FaerLuSymbolicError { message: format!("{source}") }
}

fn _matrix_err(source: faer::sparse::FaerError) -> VFVChemError {
    VFVChemError::FaerLuMatrixError { message: format!("{source}") }
}

fn _symbolic_lu_err(source: faer::sparse::FaerError) -> VFVChemError {
    VFVChemError::FaerLuSymbolicLuError { message: format!("{source}") }
}

fn _lu_err(source: faer::sparse::linalg::LuError) -> VFVChemError {
    VFVChemError::FaerLuLuError { message: format!("{source}") }
}
