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
        // rho * d(c_i)/dt + rho * vel . grad(c_i) = -div( -rho * D_ij * grad(c_j) ) + L_i * c_i + S_i
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

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, twgt: (f64, f64, f64)) -> Result<(), VFVChemError> {
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
        // discretize using time weights
        // d(c_c)/dt = tw0 * c_c + tw1 * c_prv1 + tw2 * c_prv2

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let conc_dom = &vars.scldom[self.conc_id];
        let den_dom = &vars.scldom[self.den_id];
        let (tw0, tw1, tw2) = twgt;

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let offc = tri.vert_offc[vid_v];

            // get centroid properties
            let conc_o = conc_dom.correct_vert(offc, vid_v);
            let conc_c_prv1 = conc_dom.interp_vert_prv1(offc, vid_v);
            let conc_c_prv2 = conc_dom.interp_vert_prv2(offc, vid_v);
            let den_c = den_dom.interp_vert(offc, vid_v);

            // load into matrix
            let coeff_v = -den_c * area;
            self.add_sclscl_a(vars, a_triplet, self.conc_id, vid_v, self.conc_id, vid_v, coeff_v * tw0)?;
            self.add_scl_b(vars, b_vec, self.conc_id, vid_v, coeff_v * (-tw0 * conc_o - tw1 * conc_c_prv1 - tw2 * conc_c_prv2))?;
        }

        Ok(())
    }
}

impl OperatorBase for MoleDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
