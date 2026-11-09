use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct HeatDomTime {
    // domain
    pub dom_id: usize,

    // unknowns
    pub unk_id: usize,  // unknown scalar (c_i)

    // others
    pub den_id: usize,  // density (rho)
    pub wgt_id: usize,  // weight (m)
}

impl HeatDomTime {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, den_id: usize, wgt_id: usize) -> Result<HeatDomTime, VFVChemError> {
        // adds the weighted time derivative term to the scalar transport equations
        // rho * m * d(c_i)/dt + rho * m * vel . grad(c_i) = -div( -D_ij * grad(c_j) ) + L_i c_i + R_i
        // ^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, wgt_id)?;

        // create struct
        let mut oper = HeatDomTime::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.den_id = den_id;
        oper.wgt_id = wgt_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, dt: f64) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the LHS of Ax = b
        // -rho * m * dc/dt
        //
        // get surface integral
        // int_surf -rho * m * dc/dt dA
        //
        // discretize into centroid values
        // -rho_c * m_c * d(c_c)/dt * dA
        //
        // discretize using backward euler
        // -rho_c * m_c * (c_c - c_prev)/dt * dA

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let den_dom = &vars.scldom[self.den_id];
        let wgt_dom = &vars.scldom[self.wgt_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let offc = tri.vert_offc[vid_v];

            // get centroid properties
            let unk_o = unk_dom.correct_vert(offc, vid_v);
            let unk_c_prev = unk_dom.interp_vert_prev(offc, vid_v);
            let den_c = den_dom.interp_vert(offc, vid_v);
            let wgt_c = wgt_dom.interp_vert(offc, vid_v);

            // load into matrix
            let coeff_v = -den_c * wgt_c * area / dt;
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, coeff_v)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, coeff_v * (-unk_o + unk_c_prev))?;
        }

        Ok(())
    }
}

impl OperatorBase for HeatDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
