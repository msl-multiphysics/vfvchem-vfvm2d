use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct MoleDomTime {
    // domain
    pub dom_id: usize,

    // unknowns
    pub conc_id: usize,   // concentration (c_i); mole solute per mass solvent

    // others
    pub den_id: usize,  // density (rho)
}

impl MoleDomTime {
    pub fn new(vars: &Variables, dom_id: usize, conc_id: usize, den_id: usize) -> Result<MoleDomTime, VFVChemError> {
        // adds the time derivative term to the mass transfer equation
        // rho * d(c_i)/dt + rho * vel . grad(c_i) = -div( -rho * D_ij * grad(c_j) ) + L_i * c_i + R_i
        // ^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, conc_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = MoleDomTime::default();
        oper.dom_id = dom_id;
        oper.conc_id = conc_id;
        oper.den_id = den_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, dt: f64) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the LHS of Ax = b
        // -rho * dc/dt
        //
        // get surface integral
        // int_surf -rho * dc/dt dA
        //
        // discretize into centroid values
        // -rho_c * d(c_c)/dt * dA
        //
        // discretize using backward euler
        // -rho_c * (c_c - c_prev)/dt * dA

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let conc_dom = &vars.scldom[self.conc_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let offc = tri.vert_offc[vid_v];

            // get centroid properties
            let conc_o = conc_dom.correct_vert(offc, vid_v);
            let conc_c_prev = conc_dom.interp_vert_prev(offc, vid_v);
            let den_c = den_dom.interp_vert(offc, vid_v);

            // load into matrix
            let coeff_v = -den_c * area / dt;
            self.add_sclscl_a(vars, a_triplet, self.conc_id, vid_v, self.conc_id, vid_v, coeff_v)?;
            self.add_scl_b(vars, b_vec, self.conc_id, vid_v, coeff_v * (-conc_o + conc_c_prev))?;
        }

        Ok(())
    }
}

impl OperatorBase for MoleDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
