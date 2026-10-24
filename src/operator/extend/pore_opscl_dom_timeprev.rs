use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclDomTimePrevious {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (c_i)
    pub den_id: usize,   // density (rho)
    pub epsl_id: usize,  // porosity (epsl)
}

impl PoreOpSclDomTimePrevious {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, den_id: usize, epsl_id: usize) -> Result<PoreOpSclDomTimePrevious, VFVChemError> {
        // adds the time derivative term of a porous medium with the fluid mass of the previous time step
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( ... ) + epsl * R_v,i + a_v * R_s,i
        // ^^^^^^^^^^^^^^^^^^^^
        // m_prev * (c - c_prev) / dt with m = epsl * rho; the remainder is c_v times the time derivative term of PoreOpSclDomDivergence

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclDomTimePrevious::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.den_id = den_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_surf -epsl * rho * dc/dt dA
        // -epsl_v_prev * rho_v_prev * (c_v - c_v_prev) * dA / dt

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let den_dom = &vars.scldom[self.den_id];
        let epsl_dom = &vars.scldom[self.epsl_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let wgt_v_prev = epsl_dom.value_prev[vid_v] * den_dom.value_prev[vid_v];
            let unk_v_prev = unk_dom.value_prev[vid_v];

            // load into matrix and vector
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -wgt_v_prev * area / dt)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -wgt_v_prev * unk_v_prev * area / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for PoreOpSclDomTimePrevious {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
