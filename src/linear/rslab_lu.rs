use crate::base::prelude::*;
use crate::linear::prelude::*;
use rslab::{CscLayout, LuSolver, LuSymbolic, SolverSettings, RslabError};

pub struct RslabLu {
    settings: SolverSettings,
    layout: CscLayout,
    symbolic: LuSymbolic,
    lu: LuSolver<f64>,
}

impl RslabLu {
    pub fn new(num_thread: usize) -> Result<Box<dyn LinearBase>, VFVChemError> {
        // error handling
        if num_thread == 0 { return Err(VFVChemError::LinearNonPositiveThreads { num_thread }); }
        
        // create settings
        let settings = SolverSettings::exact().with_threads(num_thread);
        
        // placeholder objects
        let layout = CscLayout::new(1, &vec![0], &vec![0]).map_err(|e| _rslab_err(e))?;
        let a_mat = layout.assemble(&vec![1.0]);
        let symbolic = LuSymbolic::analyze(&a_mat, &settings).map_err(|e| _rslab_err(e))?;
        let lu = symbolic.factor(&a_mat, &settings).map_err(|e| _rslab_err(e))?;

        // create struct
        let lin = RslabLu{settings, layout, symbolic, lu};
        Ok(Box::new(lin))

    }
}

impl LinearBase for RslabLu {
    fn set_sparse(&mut self, a_triplet: &Vec<(usize, usize, f64)>, mat_size: usize) -> Result<(), VFVChemError> {
        // collect rows, cols, and values
        let mut rows = Vec::with_capacity(a_triplet.len());
        let mut cols = Vec::with_capacity(a_triplet.len());
        let mut vals = Vec::with_capacity(a_triplet.len());
        for &(r, c, v) in a_triplet {
            rows.push(r);
            cols.push(c);
            vals.push(v);
        }
        
        // set up objects needed for reusing sparsity pattern
        let layout = CscLayout::new(mat_size, &rows, &cols).map_err(|e| _rslab_err(e))?;
        let a_mat = layout.assemble(&vals);
        let symbolic = LuSymbolic::analyze(&a_mat, &self.settings).map_err(|e| _rslab_err(e))?;
        let lu = symbolic.factor(&a_mat, &self.settings).map_err(|e| _rslab_err(e))?;

        // store in struct
        self.layout = layout;
        self.symbolic = symbolic;
        self.lu = lu;

        Ok(())
    }

    fn solve(&mut self, a_triplet: &Vec<(usize, usize, f64)>, b_vec: &Vec<f64>, _x_init: &Vec<f64>) -> Result<Vec<f64>, VFVChemError> {
        // assemble matrix and solve
        let vals: Vec<f64> = a_triplet.iter().map(|&(_, _, v)| v).collect();   
        let a_mat = self.layout.assemble(&vals);
        self.symbolic.refactor(&a_mat, &self.settings, &mut self.lu).map_err(|e| _rslab_err(e))?;
        let x_vec = self.lu.solve(b_vec).map_err(|e| _rslab_err(e))?;
        Ok(x_vec)
    }
}

fn _rslab_err(source: RslabError) -> VFVChemError {
    VFVChemError::LinearRslabError{ source }
}
