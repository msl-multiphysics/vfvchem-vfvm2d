use crate::base::prelude::*;

pub trait LinearBase {
    // to implement in specific solvers
    fn set_sparse(&mut self, a_triplet: &Vec<(usize, usize, f64)>, mat_size: usize) -> Result<(), VFVChemError>;
    fn solve(&mut self, a_triplet: &Vec<(usize, usize, f64)>, b_vec: &Vec<f64>, x_init: &Vec<f64>) -> Result<Vec<f64>, VFVChemError>;
}
