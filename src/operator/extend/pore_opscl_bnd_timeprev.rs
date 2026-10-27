use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclBndTimePrevious {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (c_i)
    pub den_id: usize,   // density (rho)
    pub epsl_id: usize,  // porosity (epsl)
}

impl PoreOpSclBndTimePrevious {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, den_id: usize, epsl_id: usize) -> Result<PoreOpSclBndTimePrevious, VFVChemError> {
        // adds the time derivative term of a porous medium with the fluid mass of the previous time step to the control volumes of boundary nodes
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( ... ) + epsl * R_v,i + a_v * R_s,i
        // ^^^^^^^^^^^^^^^^^^^^
        // m_prev * (c - c_prev) / dt with m = epsl * rho; the remainder is c_v times the time derivative term of PoreOpSclDomDivergence

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;
        check_bnd_scldom(vars, bnd_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclBndTimePrevious::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.den_id = den_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_surf -epsl * rho * dc/dt dA
        // -epsl_v_prev * rho_v_prev * (c_v - c_v_prev) * dA / dt

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let den_dom = &vars.scldom[self.den_id];
        let epsl_dom = &vars.scldom[self.epsl_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

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

impl OperatorBase for PoreOpSclBndTimePrevious {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
