use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpVecItfTime {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,   // superficial velocity (vel^a); unknown
    pub den_a_id: usize,   // density (rho^a)
    pub epsl_a_id: usize,  // porosity (epsl^a)

    // domain b
    pub vel_b_id: usize,   // superficial velocity (vel^b); unknown
    pub den_b_id: usize,   // density (rho^b)
    pub epsl_b_id: usize,  // porosity (epsl^b)
}

impl PoreOpVecItfTime {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, den_a_id: usize, epsl_a_id: usize, vel_b_id: usize, den_b_id: usize, epsl_b_id: usize) -> Result<PoreOpVecItfTime, VFVChemError> {
        // adds the time derivative term to the combined control volumes of interface nodes
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        // ^^^^^^^^^^^^^^^^^^^^^^^^
        // each domain contributes its half of the control volume; loads into the rows of domain b

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, den_a_id, den_b_id)?;
        check_itf_scldom(vars, itf_id, epsl_a_id, epsl_b_id)?;

        // create struct
        let mut oper = PoreOpVecItfTime::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.den_a_id = den_a_id;
        oper.epsl_a_id = epsl_a_id;
        oper.vel_b_id = vel_b_id;
        oper.den_b_id = den_b_id;
        oper.epsl_b_id = epsl_b_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_surf -(rho / epsl) * d(vel)/dt dA
        // -(rho_v / epsl_v) * (vel_v - vel_v_prev) * dA / dt

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, den_id, epsl_id) in [(self.vel_a_id, self.den_a_id, self.epsl_a_id), (self.vel_b_id, self.den_b_id, self.epsl_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];
            let vel_dom = &vars.vecdom[vel_id];
            let den_dom = &vars.scldom[den_id];
            let epsl_dom = &vars.scldom[epsl_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let area = tri.cvol_area[vid_v];

                // get vertex properties
                let den_v = den_dom.value[vid_v] / epsl_dom.value[vid_v];
                let (u_v_prev, v_v_prev) = vel_dom.value_prev[vid_v];

                // load into matrix and vector
                self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 0, vid_b_v, vel_id, 0, vid_v, -den_v * area / dt)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 1, vid_b_v, vel_id, 1, vid_v, -den_v * area / dt)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 0, vid_b_v, -den_v * u_v_prev * area / dt)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 1, vid_b_v, -den_v * v_v_prev * area / dt)?;
            }
        }

        Ok(())
    }
}

impl OperatorBase for PoreOpVecItfTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
