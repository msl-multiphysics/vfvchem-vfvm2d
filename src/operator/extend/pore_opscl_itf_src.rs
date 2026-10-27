use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclItfSource {
    // interface
    pub itf_id: usize,

    // scalars - domain a
    pub unk_a_id: usize,  // unknown scalar (c_i^a)
    pub src_a_id: usize,  // source (R_i^a)
    pub wgt_a_id: usize,  // weight (epsl^a or a_v^a)

    // scalars - domain b
    pub unk_b_id: usize,  // unknown scalar (c_i^b)
    pub src_b_id: usize,  // source (R_i^b)
    pub wgt_b_id: usize,  // weight (epsl^b or a_v^b)
}

impl PoreOpSclItfSource {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, src_a_id: usize, wgt_a_id: usize, unk_b_id: usize, src_b_id: usize, wgt_b_id: usize) -> Result<PoreOpSclItfSource, VFVChemError> {
        // adds the weighted source term of a porous medium to the combined control volumes of interface nodes
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( ... ) + epsl * (s_v,i * c_i + R_v,i) + a_v * (s_s,i * c_i + R_s,i)
        //                                                                                    ^^^^^                         ^^^^^
        // the weight is epsl for homogeneous (per fluid volume) and a_v for heterogeneous (per solid surface area) reactions
        // each domain contributes its half of the control volume; loads into the rows of domain b

        // input checking
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_itf_scldom(vars, itf_id, src_a_id, src_b_id)?;
        check_itf_scldom(vars, itf_id, wgt_a_id, wgt_b_id)?;

        // create struct
        let mut oper = PoreOpSclItfSource::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.src_a_id = src_a_id;
        oper.wgt_a_id = wgt_a_id;
        oper.unk_b_id = unk_b_id;
        oper.src_b_id = src_b_id;
        oper.wgt_b_id = wgt_b_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclItfSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf w * R dA
        // (w * R)_c * dA; centroid skewness correction as in OpSclDomSource with w * R in place of R

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];

        // iterate through both sides of the interface
        for (src_id, wgt_id) in [(self.src_a_id, self.wgt_a_id), (self.src_b_id, self.wgt_b_id)] {
            // get objects on this side
            let dom_id = vars.scldom[src_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let (sx, sy) = tri.cvol_skew[vid_v];
                let area = tri.cvol_area[vid_v];

                // get vertex properties
                let (src_v, (src_x_v, src_y_v)) = pore_helper::weighted(vars, src_id, wgt_id, vid_v);

                // get source at centroid
                let src_c = src_v + src_x_v * sx + src_y_v * sy;

                // load into vector
                self.add_scl_b(vars, b_vec, self.unk_b_id, vid_b_v, -src_c * area)?;
            }
        }

        Ok(())
    }
}
