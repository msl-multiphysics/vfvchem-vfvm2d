use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclItfSolidTimePrevious {
    // interface
    pub itf_id: usize,

    // scalars - rows
    pub row_id: usize,  // row scalar (T_s^b with continuity; T_s on this side without)

    // scalars - this side
    pub unk_id: usize,   // unknown scalar (T_s)
    pub wgt_id: usize,   // solid volumetric heat capacity (rho * c)_s
    pub epsl_id: usize,  // porosity (epsl); weight is 1 - epsl
}

impl PoreOpSclItfSolidTimePrevious {
    pub fn new(vars: &Variables, itf_id: usize, unk_id: usize, row_id: usize, wgt_id: usize, epsl_id: usize) -> Result<PoreOpSclItfSolidTimePrevious, VFVChemError> {
        // adds the solid time derivative term of a porous medium to the interface nodes
        // (1 - epsl) * (rho * c)_s * dT_s/dt = ...
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // same as PoreOpSclDomSolidTimePrevious; covers the half of the control volume on the side of the unknown; loads into the rows of row_id
        // row_id is T_s^b for a solid temperature with continuity (use once for each side), or unk_id for a solid temperature without

        // input checking
        let dom_id = vars.scldom[unk_id].dom_id;
        let itf = &vars.itf[itf_id];
        if row_id != unk_id { check_itf_scldom(vars, itf_id, unk_id, row_id)?; }
        else if dom_id != itf.dom_a_id && dom_id != itf.dom_b_id { return Err(VFVChemError::OperatorScalarDomainNotOnInterface { scldom_a_id: unk_id, scldom_b_id: row_id, itf_id }) }
        check_dom_scldom(vars, dom_id, wgt_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclItfSolidTimePrevious::default();
        oper.itf_id = itf_id;
        oper.row_id = row_id;
        oper.unk_id = unk_id;
        oper.wgt_id = wgt_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_surf -(1 - epsl) * m * dT/dt dA
        // -((1 - epsl) * m)_v_prev * (T_v - T_v_prev) * dA / dt

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_id = vars.scldom[self.unk_id].dom_id;
        let dom = &vars.dom[dom_id];
        let dom_r = &vars.dom[vars.scldom[self.row_id].dom_id];
        let tri = &vars.tri[dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let wgt_dom = &vars.scldom[self.wgt_id];
        let epsl_dom = &vars.scldom[self.epsl_id];

        // iterate through nodes
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
            let vid_r_v = dom_r.num_pts() + dom_r.nds_nid[&gnd_v];

            // get geometry
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let wgt_v_prev = (1.0 - epsl_dom.value_prev[vid_v]) * wgt_dom.value_prev[vid_v];
            let unk_v_prev = unk_dom.value_prev[vid_v];

            // load into matrix and vector
            self.add_sclscl_a(vars, a_triplet, self.row_id, vid_r_v, self.unk_id, vid_v, -wgt_v_prev * area / dt)?;
            self.add_scl_b(vars, b_vec, self.row_id, vid_r_v, -wgt_v_prev * unk_v_prev * area / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for PoreOpSclItfSolidTimePrevious {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
