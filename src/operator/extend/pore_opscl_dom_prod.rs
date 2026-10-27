use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclDomProduction {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub prd_id: usize,  // production coefficient (s_i)
    pub wgt_id: usize,  // weight (epsl or a_v)
}

impl PoreOpSclDomProduction {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, prd_id: usize, wgt_id: usize) -> Result<PoreOpSclDomProduction, VFVChemError> {
        // adds the weighted production term of a porous medium
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( ... ) + epsl * (s_v,i * c_i + R_v,i) + a_v * (s_s,i * c_i + R_s,i)
        //                                                                      ^^^^^^^^^^^                   ^^^^^^^^^^^
        // the weight is epsl for homogeneous (per fluid volume) and a_v for heterogeneous (per solid surface area) reactions
        // s_i is lagged (evaluated with the previous iterate); s_i * c_i is implicit on the diagonal

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, prd_id)?;
        check_dom_scldom(vars, dom_id, wgt_id)?;

        // create struct
        let mut oper = PoreOpSclDomProduction::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.prd_id = prd_id;
        oper.wgt_id = wgt_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclDomProduction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf w * s * c dA
        // (w * s)_c * c_c * dA; centroid skewness correction as in OpSclDomProduction with w * s in place of s

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let (sx, sy) = tri.cvol_skew[vid_v];
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let (prd_v, (prd_x_v, prd_y_v)) = pore_helper::weighted(vars, self.prd_id, self.wgt_id, vid_v);
            let (unk_x_v, unk_y_v) = unk_dom.grad[vid_v];

            // get production coefficient at centroid and skewness correction of the unknown
            let prd_c = prd_v + prd_x_v * sx + prd_y_v * sy;
            let unk_skew = unk_x_v * sx + unk_y_v * sy;

            // load into matrix and vector
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, prd_c * area)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -prd_v * unk_skew * area)?;
        }

        Ok(())
    }
}
