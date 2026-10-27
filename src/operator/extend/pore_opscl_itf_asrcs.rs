use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclItfSurfaceSource {
    // interface
    pub itf_id: usize,

    // scalars - rows
    pub row_id: usize,  // row scalar (T_s^b with continuity; T_s on this side without)

    // scalars - this side
    pub unk_id: usize,   // unknown scalar (T_s)
    pub src_id: usize,   // source per solid surface area (q_a)
    pub area_id: usize,  // specific surface area (a_v)
}

impl PoreOpSclItfSurfaceSource {
    pub fn new(vars: &Variables, itf_id: usize, unk_id: usize, row_id: usize, src_id: usize, area_id: usize) -> Result<PoreOpSclItfSurfaceSource, VFVChemError> {
        // adds the solid surface source term of a porous medium to the interface nodes
        // (1 - epsl) * (rho * c)_s * dT_s/dt = ... + a_v * (s_a * T_s + q_a)
        //                                                               ^^^
        // same as PoreOpSclDomSource with the weight a_v; covers the half of the control volume on the side of the unknown; loads into the rows of row_id
        // row_id is T_s^b for a solid temperature with continuity (use once for each side), or unk_id for a solid temperature without

        // input checking
        let dom_id = vars.scldom[unk_id].dom_id;
        let itf = &vars.itf[itf_id];
        if row_id != unk_id { check_itf_scldom(vars, itf_id, unk_id, row_id)?; }
        else if dom_id != itf.dom_a_id && dom_id != itf.dom_b_id { return Err(VFVChemError::OperatorScalarDomainNotOnInterface { scldom_a_id: unk_id, scldom_b_id: row_id, itf_id }) }
        check_dom_scldom(vars, dom_id, src_id)?;
        check_dom_scldom(vars, dom_id, area_id)?;

        // create struct
        let mut oper = PoreOpSclItfSurfaceSource::default();
        oper.itf_id = itf_id;
        oper.row_id = row_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;
        oper.area_id = area_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclItfSurfaceSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf a_v * q dA
        // (a_v * q)_c * dA; centroid skewness correction as in OpSclDomSource

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_id = vars.scldom[self.unk_id].dom_id;
        let dom = &vars.dom[dom_id];
        let dom_r = &vars.dom[vars.scldom[self.row_id].dom_id];
        let tri = &vars.tri[dom_id];

        // iterate through nodes
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
            let vid_r_v = dom_r.num_pts() + dom_r.nds_nid[&gnd_v];

            // get geometry
            let (sx, sy) = tri.cvol_skew[vid_v];
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let (src_v, (src_x_v, src_y_v)) = pore_helper::weighted(vars, self.src_id, self.area_id, vid_v);

            // get source at centroid
            let src_c = src_v + src_x_v * sx + src_y_v * sy;

            // load into vector
            self.add_scl_b(vars, b_vec, self.row_id, vid_r_v, -src_c * area)?;
        }

        Ok(())
    }
}
