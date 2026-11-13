use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct HeatDomTime {
    // domain
    pub dom_id: usize,

    // unknowns
    pub temp_id: usize,   // temperature (T)

    // others
    pub den_id: usize,  // density (rho)
    pub hcap_id: usize,  // heat capacity (cp)
}

impl HeatDomTime {
    pub fn new(vars: &Variables, dom_id: usize, temp_id: usize, den_id: usize, hcap_id: usize) -> Result<HeatDomTime, VFVChemError> {
        // adds the time derivative term to the heat transfer equation
        // rho * cp * d(T)/dt + rho * cp * vel . grad(T) = -div( -k * grad(T) ) + L * T + R
        // ^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, temp_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_scldom(vars, dom_id, hcap_id)?;

        // create struct
        let mut oper = HeatDomTime::default();
        oper.dom_id = dom_id;
        oper.temp_id = temp_id;
        oper.den_id = den_id;
        oper.hcap_id = hcap_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, twgt: (f64, f64, f64)) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the LHS of Ax = b
        // -rho * cp * d(T)/dt
        //
        // get surface integral
        // int_surf -rho * cp * d(T)/dt dA
        //
        // discretize into centroid values
        // -rho_c * cp_c * d(T_c)/dt * dA
        //
        // discretize using time weights
        // d(T_c)/dt = tw0 * T_c + tw1 * T_prv1 + tw2 * T_prv2

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let temp_dom = &vars.scldom[self.temp_id];
        let den_dom = &vars.scldom[self.den_id];
        let hcap_dom = &vars.scldom[self.hcap_id];
        let (tw0, tw1, tw2) = twgt;

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let offc = tri.vert_offc[vid_v];

            // get centroid properties
            let temp_o = temp_dom.correct_vert(offc, vid_v);
            let temp_c_prv1 = temp_dom.interp_vert_prv1(offc, vid_v);
            let temp_c_prv2 = temp_dom.interp_vert_prv2(offc, vid_v);
            let den_c = den_dom.interp_vert(offc, vid_v);
            let hcap_c = hcap_dom.interp_vert(offc, vid_v);

            // load into matrix
            let coeff_v = -den_c * hcap_c * area;
            self.add_sclscl_a(vars, a_triplet, self.temp_id, vid_v, self.temp_id, vid_v, coeff_v * tw0)?;
            self.add_scl_b(vars, b_vec, self.temp_id, vid_v, coeff_v * (-tw0 * temp_o - tw1 * temp_c_prv1 - tw2 * temp_c_prv2))?;
        }

        Ok(())
    }
}

impl OperatorBase for HeatDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
