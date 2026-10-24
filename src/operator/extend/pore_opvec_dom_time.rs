use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpVecDomTime {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // superficial velocity (vel); unknown

    // scalars
    pub den_id: usize,   // density (rho)
    pub epsl_id: usize,  // porosity (epsl)
}

impl PoreOpVecDomTime {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, den_id: usize, epsl_id: usize) -> Result<PoreOpVecDomTime, VFVChemError> {
        // adds the time derivative term to the Darcy-Brinkman-Forchheimer equation
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        // ^^^^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpVecDomTime::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_surf -(rho / epsl) * d(vel)/dt dA
        // -(rho_v / epsl_v) * (vel_v - vel_v_prev) * dA / dt

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        let epsl_dom = &vars.scldom[self.epsl_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let den_v = den_dom.value[vid_v] / epsl_dom.value[vid_v];
            let (u_v_prev, v_v_prev) = vel_dom.value_prev[vid_v];

            // load into matrix and vector
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v, -den_v * area / dt)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v, -den_v * area / dt)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -den_v * u_v_prev * area / dt)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -den_v * v_v_prev * area / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for PoreOpVecDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
