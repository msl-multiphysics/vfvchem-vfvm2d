use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclBndHeatExchange {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (T_f or T_s)
    pub drv_id: usize,   // other phase temperature (T_s or T_f)
    pub htc_id: usize,   // interphase heat transfer coefficient (h_sf)
    pub area_id: usize,  // specific surface area (a_v)
    pub epsl_id: usize,  // porosity (epsl); for the solid mask
}

impl PoreOpSclBndHeatExchange {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, drv_id: usize, htc_id: usize, area_id: usize, epsl_id: usize) -> Result<PoreOpSclBndHeatExchange, VFVChemError> {
        // adds the interphase heat exchange term of a porous medium
        // fluid: ... = ... + h_sf * a_v * (T_s - T_f) + ...
        // solid: ... = ... - h_sf * a_v * (T_s - T_f) + ...
        //                  ^^^^^^^^^^^^^^^^^^^^^^^^^^
        // both are h_sf * a_v * (T_drv - T_unk); implicit in both temperatures
        // skipped on vertices where the solid is masked (T_s = T_f)

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, drv_id)?;
        check_bnd_scldom(vars, bnd_id, htc_id)?;
        check_bnd_scldom(vars, bnd_id, area_id)?;
        check_bnd_scldom(vars, bnd_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclBndHeatExchange::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.drv_id = drv_id;
        oper.htc_id = htc_id;
        oper.area_id = area_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclBndHeatExchange {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf h * a_v * (T_drv - T_unk) dA
        // (h * a_v)_v * (T_drv,v - T_unk,v) * dA

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let htc_dom = &vars.scldom[self.htc_id];
        let area_dom = &vars.scldom[self.area_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // skip masked vertices
            if pore_helper::solid_masked(vars, bnd.dom_id, self.epsl_id, vid_v) { continue; }

            // get geometry
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let coeff = htc_dom.value[vid_v] * area_dom.value[vid_v] * area;

            // load into matrix
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.drv_id, vid_v,  coeff)?;
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -coeff)?;
        }

        Ok(())
    }
}
