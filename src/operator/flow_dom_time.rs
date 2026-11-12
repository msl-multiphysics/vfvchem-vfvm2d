use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct FlowDomTime {
    // domain
    pub dom_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub den_id: usize,  // density (rho)
}

impl FlowDomTime {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, den_id: usize) -> Result<FlowDomTime, VFVChemError> {
        // adds the advection term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        // ^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = FlowDomTime::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, tc: (f64, f64, f64)) -> Result<(), VFVChemError> {
        // the LHS of the PDE is for the time derivative
        // the RHS of the PDE is discretized into the LHS of Ax = b
        // invert the sign to transfer the term to the LHS of Ax = b
        // -rho * d(vel)/dt
        //
        // get surface integral
        // int_surf -rho * d(vel)/dt dA
        //
        // discretize into centroid values
        // -rho_c * d(vel_c)/dt * dA
        //
        // discretize using time coefficients tc = (tc0, tc1, tc2)
        // bdf1: tc = (1, -1, 0) / dt; bdf2: tc = (3/2, -2, 1/2) / dt
        // -rho_c * (tc0 * vel_c + tc1 * vel_prv1 + tc2 * vel_prv2) * dA

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        let (tc0, tc1, tc2) = tc;

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let offc = tri.vert_offc[vid_v];

            // get centroid properties
            let (velx_o, vely_o) = vel_dom.correct_vert(offc, vid_v);
            let (velx_c_prv1, vely_c_prv1) = vel_dom.interp_vert_prv1(offc, vid_v);
            let (velx_c_prv2, vely_c_prv2) = vel_dom.interp_vert_prv2(offc, vid_v);
            let den_c = den_dom.interp_vert(offc, vid_v);

            // load into matrix
            let coeff_v = -den_c * area;
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v, coeff_v * tc0)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v, coeff_v * tc0)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, coeff_v * (-tc0 * velx_o - tc1 * velx_c_prv1 - tc2 * velx_c_prv2))?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, coeff_v * (-tc0 * vely_o - tc1 * vely_c_prv1 - tc2 * vely_c_prv2))?;
        }

        Ok(())
    }
}

impl OperatorBase for FlowDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
