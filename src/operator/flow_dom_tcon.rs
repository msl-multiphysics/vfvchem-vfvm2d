use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct FlowDomTimeContinuity {
    // domain
    pub dom_id: usize,

    // unknowns
    pub pres_id: usize,  // pressure (p)

    // others
    pub den_id: usize,  // density (rho)
}

impl FlowDomTimeContinuity {
    pub fn new(vars: &Variables, dom_id: usize, pres_id: usize, den_id: usize) -> Result<FlowDomTimeContinuity, VFVChemError> {
        // adds the time derivative term to the continuity equation
        // d(rho)/dt = -div(rho * vel)
        // ^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, pres_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = FlowDomTimeContinuity::default();
        oper.dom_id = dom_id;
        oper.pres_id = pres_id;
        oper.den_id = den_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, dt: f64) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the LHS of Ax = b
        // -d(rho)/dt
        //
        // get surface integral
        // int_surf -d(rho)/dt dA
        //
        // discretize into centroid values
        // -d(rho_c)/dt * dA
        //
        // discretize using backward euler
        // -(rho_c - rho_prev)/dt * dA

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let offc = tri.vert_offc[vid_v];

            // get centroid properties
            let den_c = den_dom.interp_vert(offc, vid_v);
            let den_c_prev = den_dom.interp_vert_prev(offc, vid_v);

            // load into matrix
            self.add_scl_b(vars, b_vec, self.pres_id, vid_v, (den_c - den_c_prev) * area / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for FlowDomTimeContinuity {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
