use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclBndSource {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub src_id: usize,  // source (R_i)
    pub wgt_id: usize,  // weight (epsl or a_v)
}

impl PoreOpSclBndSource {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, src_id: usize, wgt_id: usize) -> Result<PoreOpSclBndSource, VFVChemError> {
        // adds the weighted source term of a porous medium to the control volumes of boundary nodes
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( ... ) + epsl * (s_v,i * c_i + R_v,i) + a_v * (s_s,i * c_i + R_s,i)
        //                                                                                    ^^^^^                         ^^^^^
        // the weight is epsl for homogeneous (per fluid volume) and a_v for heterogeneous (per solid surface area) reactions

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, src_id)?;
        check_bnd_scldom(vars, bnd_id, wgt_id)?;

        // create struct
        let mut oper = PoreOpSclBndSource::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;
        oper.wgt_id = wgt_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclBndSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf w * R dA
        // (w * R)_c * dA; centroid skewness correction as in OpSclDomSource with w * R in place of R

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let (sx, sy) = tri.cvol_skew[vid_v];
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let (src_v, (src_x_v, src_y_v)) = pore_helper::weighted(vars, self.src_id, self.wgt_id, vid_v);

            // get source at centroid
            let src_c = src_v + src_x_v * sx + src_y_v * sy;

            // load into vector
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -src_c * area)?;
        }

        Ok(())
    }
}
