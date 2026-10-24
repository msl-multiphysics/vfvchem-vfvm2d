use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclItfProduction {
    // interface
    pub itf_id: usize,

    // scalars - domain a
    pub unk_a_id: usize,  // unknown scalar (c_i^a)
    pub prd_a_id: usize,  // production coefficient (s_i^a)
    pub wgt_a_id: usize,  // weight (epsl^a or a_v^a)

    // scalars - domain b
    pub unk_b_id: usize,  // unknown scalar (c_i^b)
    pub prd_b_id: usize,  // production coefficient (s_i^b)
    pub wgt_b_id: usize,  // weight (epsl^b or a_v^b)
}

impl PoreOpSclItfProduction {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, prd_a_id: usize, wgt_a_id: usize, unk_b_id: usize, prd_b_id: usize, wgt_b_id: usize) -> Result<PoreOpSclItfProduction, VFVChemError> {
        // adds the weighted production term of a porous medium to the combined control volumes of interface nodes
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( ... ) + epsl * (s_v,i * c_i + R_v,i) + a_v * (s_s,i * c_i + R_s,i)
        //                                                                      ^^^^^^^^^^^                   ^^^^^^^^^^^
        // the weight is epsl for homogeneous (per fluid volume) and a_v for heterogeneous (per solid surface area) reactions
        // s_i is lagged (evaluated with the previous iterate); s_i * c_i is implicit on the diagonal
        // each domain contributes its half of the control volume; loads into the rows of domain b
        // c_i^a = c_i^b by continuity, so both halves act on c_i^b

        // input checking
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_itf_scldom(vars, itf_id, prd_a_id, prd_b_id)?;
        check_itf_scldom(vars, itf_id, wgt_a_id, wgt_b_id)?;

        // create struct
        let mut oper = PoreOpSclItfProduction::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.prd_a_id = prd_a_id;
        oper.wgt_a_id = wgt_a_id;
        oper.unk_b_id = unk_b_id;
        oper.prd_b_id = prd_b_id;
        oper.wgt_b_id = wgt_b_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclItfProduction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf w * s * c dA
        // (w * s)_c * c_c * dA; centroid skewness correction as in OpSclDomProduction with w * s in place of s

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];

        // iterate through both sides of the interface
        for (unk_id, prd_id, wgt_id) in [(self.unk_a_id, self.prd_a_id, self.wgt_a_id), (self.unk_b_id, self.prd_b_id, self.wgt_b_id)] {
            // get objects on this side
            let dom_id = vars.scldom[prd_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];
            let unk_dom = &vars.scldom[unk_id];

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
                let (prd_v, (prd_x_v, prd_y_v)) = pore_helper::weighted(vars, prd_id, wgt_id, vid_v);
                let (unk_x_v, unk_y_v) = unk_dom.grad[vid_v];

                // get production coefficient at centroid and skewness correction of the unknown
                let prd_c = prd_v + prd_x_v * sx + prd_y_v * sy;
                let unk_skew = unk_x_v * sx + unk_y_v * sy;

                // load into matrix and vector
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_b_id, vid_b_v, prd_c * area)?;
                self.add_scl_b(vars, b_vec, self.unk_b_id, vid_b_v, -prd_v * unk_skew * area)?;
            }
        }

        Ok(())
    }
}
