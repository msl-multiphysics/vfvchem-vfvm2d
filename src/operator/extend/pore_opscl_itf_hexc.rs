use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclItfHeatExchange {
    // interface
    pub itf_id: usize,

    // scalars - rows
    pub row_id: usize,  // row scalar (T^b with continuity; T on this side without)

    // scalars - this side
    pub unk_id: usize,   // unknown scalar (T_f or T_s)
    pub drv_id: usize,   // other phase temperature (T_s or T_f)
    pub htc_id: usize,   // interphase heat transfer coefficient (h_sf)
    pub area_id: usize,  // specific surface area (a_v)
    pub epsl_id: usize,  // porosity (epsl); for the solid mask
}

impl PoreOpSclItfHeatExchange {
    pub fn new(vars: &Variables, itf_id: usize, unk_id: usize, row_id: usize, drv_id: usize, htc_id: usize, area_id: usize, epsl_id: usize) -> Result<PoreOpSclItfHeatExchange, VFVChemError> {
        // adds the interphase heat exchange term of a porous medium to the interface nodes
        // fluid: ... = ... + h_sf * a_v * (T_s - T_f) + ...
        // solid: ... = ... - h_sf * a_v * (T_s - T_f) + ...
        //                  ^^^^^^^^^^^^^^^^^^^^^^^^^^
        // same as PoreOpSclDomHeatExchange; covers the half of the control volume on the side of the unknown; loads into the rows of row_id
        // row_id is the unknown on domain b for a temperature with continuity (use once for each side), or unk_id for a temperature without
        // skipped on vertices where the solid is masked (T_s = T_f)

        // input checking
        let dom_id = vars.scldom[unk_id].dom_id;
        let itf = &vars.itf[itf_id];
        if row_id != unk_id { check_itf_scldom(vars, itf_id, unk_id, row_id)?; }
        else if dom_id != itf.dom_a_id && dom_id != itf.dom_b_id { return Err(VFVChemError::OperatorScalarDomainNotOnInterface { scldom_a_id: unk_id, scldom_b_id: row_id, itf_id }) }
        check_dom_scldom(vars, dom_id, drv_id)?;
        check_dom_scldom(vars, dom_id, htc_id)?;
        check_dom_scldom(vars, dom_id, area_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclItfHeatExchange::default();
        oper.itf_id = itf_id;
        oper.row_id = row_id;
        oper.unk_id = unk_id;
        oper.drv_id = drv_id;
        oper.htc_id = htc_id;
        oper.area_id = area_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclItfHeatExchange {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf h * a_v * (T_drv - T_unk) dA
        // (h * a_v)_v * (T_drv,v - T_unk,v) * dA

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_id = vars.scldom[self.unk_id].dom_id;
        let dom = &vars.dom[dom_id];
        let dom_r = &vars.dom[vars.scldom[self.row_id].dom_id];
        let tri = &vars.tri[dom_id];
        let htc_dom = &vars.scldom[self.htc_id];
        let area_dom = &vars.scldom[self.area_id];

        // iterate through nodes
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
            let vid_r_v = dom_r.num_pts() + dom_r.nds_nid[&gnd_v];

            // skip masked vertices
            if pore_helper::solid_masked(vars, dom_id, self.epsl_id, vid_v) { continue; }

            // get geometry
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let coeff = htc_dom.value[vid_v] * area_dom.value[vid_v] * area;

            // load into matrix
            self.add_sclscl_a(vars, a_triplet, self.row_id, vid_r_v, self.drv_id, vid_v,  coeff)?;
            self.add_sclscl_a(vars, a_triplet, self.row_id, vid_r_v, self.unk_id, vid_v, -coeff)?;
        }

        Ok(())
    }
}
