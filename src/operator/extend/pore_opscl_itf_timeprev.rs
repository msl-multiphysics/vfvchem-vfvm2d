use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclItfTimePrevious {
    // interface
    pub itf_id: usize,

    // scalars - domain a
    pub unk_a_id: usize,   // unknown scalar (c_i^a)
    pub den_a_id: usize,   // density (rho^a)
    pub epsl_a_id: usize,  // porosity (epsl^a)

    // scalars - domain b
    pub unk_b_id: usize,   // unknown scalar (c_i^b)
    pub den_b_id: usize,   // density (rho^b)
    pub epsl_b_id: usize,  // porosity (epsl^b)
}

impl PoreOpSclItfTimePrevious {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, den_a_id: usize, epsl_a_id: usize, unk_b_id: usize, den_b_id: usize, epsl_b_id: usize) -> Result<PoreOpSclItfTimePrevious, VFVChemError> {
        // adds the time derivative term of a porous medium with the fluid mass of the previous time step to the combined control volumes of interface nodes
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( ... ) + epsl * R_v,i + a_v * R_s,i
        // ^^^^^^^^^^^^^^^^^^^^
        // m_prev * (c - c_prev) / dt with m = epsl * rho; the remainder is c_v times the time derivative term of PoreOpSclItfDivergence
        // each domain contributes its half of the control volume; loads into the rows of domain b

        // input checking
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_itf_scldom(vars, itf_id, den_a_id, den_b_id)?;
        check_itf_scldom(vars, itf_id, epsl_a_id, epsl_b_id)?;

        // create struct
        let mut oper = PoreOpSclItfTimePrevious::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.den_a_id = den_a_id;
        oper.epsl_a_id = epsl_a_id;
        oper.unk_b_id = unk_b_id;
        oper.den_b_id = den_b_id;
        oper.epsl_b_id = epsl_b_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_surf -epsl * rho * dc/dt dA
        // -epsl_v_prev * rho_v_prev * (c_v - c_v_prev) * dA / dt

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];

        // iterate through both sides of the interface
        for (unk_id, den_id, epsl_id) in [(self.unk_a_id, self.den_a_id, self.epsl_a_id), (self.unk_b_id, self.den_b_id, self.epsl_b_id)] {
            // get objects on this side
            let dom_id = vars.scldom[unk_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];
            let unk_dom = &vars.scldom[unk_id];
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
                let wgt_v_prev = epsl_dom.value_prev[vid_v] * den_dom.value_prev[vid_v];
                let unk_v_prev = unk_dom.value_prev[vid_v];

                // load into matrix and vector
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, unk_id, vid_v, -wgt_v_prev * area / dt)?;
                self.add_scl_b(vars, b_vec, self.unk_b_id, vid_b_v, -wgt_v_prev * unk_v_prev * area / dt)?;
            }
        }

        Ok(())
    }
}

impl OperatorBase for PoreOpSclItfTimePrevious {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
