use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclBndProduction {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub prd_id: usize,  // production coefficient (s_i)
    pub wgt_id: usize,  // weight (epsl or a_v)
}

impl PoreOpSclBndProduction {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, prd_id: usize, wgt_id: usize) -> Result<PoreOpSclBndProduction, VFVChemError> {
        // adds the weighted production term of a porous medium to the control volumes of boundary nodes
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( ... ) + epsl * (s_v,i * c_i + R_v,i) + a_v * (s_s,i * c_i + R_s,i)
        //                                                                      ^^^^^^^^^^^                   ^^^^^^^^^^^
        // the weight is epsl for homogeneous (per fluid volume) and a_v for heterogeneous (per solid surface area) reactions
        // s_i is lagged (evaluated with the previous iterate); s_i * c_i is implicit on the diagonal

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, prd_id)?;
        check_bnd_scldom(vars, bnd_id, wgt_id)?;

        // create struct
        let mut oper = PoreOpSclBndProduction::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.prd_id = prd_id;
        oper.wgt_id = wgt_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclBndProduction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf w * s * c dA
        // (w * s)_c * c_c * dA; centroid skewness correction as in OpSclDomProduction with w * s in place of s

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];

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
