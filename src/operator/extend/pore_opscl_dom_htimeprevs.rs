use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclDomSolidTimePrevious {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (T_s)
    pub wgt_id: usize,   // solid volumetric heat capacity (rho * c)_s
    pub epsl_id: usize,  // porosity (epsl); weight is 1 - epsl
}

impl PoreOpSclDomSolidTimePrevious {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, wgt_id: usize, epsl_id: usize) -> Result<PoreOpSclDomSolidTimePrevious, VFVChemError> {
        // adds the solid time derivative term of a porous medium
        // (1 - epsl) * (rho * c)_s * dT_s/dt = ...
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // same as OpSclDomTimePrevious with the weight (1 - epsl) * (rho * c)_s at the previous time

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, wgt_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclDomSolidTimePrevious::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.wgt_id = wgt_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_surf -(1 - epsl) * m * dT/dt dA
        // -((1 - epsl) * m)_v_prev * (T_v - T_v_prev) * dA / dt

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let wgt_dom = &vars.scldom[self.wgt_id];
        let epsl_dom = &vars.scldom[self.epsl_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let wgt_v_prev = (1.0 - epsl_dom.value_prev[vid_v]) * wgt_dom.value_prev[vid_v];
            let unk_v_prev = unk_dom.value_prev[vid_v];

            // load into matrix and vector
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -wgt_v_prev * area / dt)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -wgt_v_prev * unk_v_prev * area / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for PoreOpSclDomSolidTimePrevious {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
