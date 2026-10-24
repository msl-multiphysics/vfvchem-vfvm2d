use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclDomSource {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub src_id: usize,  // source (R_i)
    pub wgt_id: usize,  // weight (epsl or a_v)
}

impl PoreOpSclDomSource {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, src_id: usize, wgt_id: usize) -> Result<PoreOpSclDomSource, VFVChemError> {
        // adds the weighted source term of a porous medium
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( ... ) + epsl * (s_v,i * c_i + R_v,i) + a_v * (s_s,i * c_i + R_s,i)
        //                                                                                    ^^^^^                         ^^^^^
        // the weight is epsl for homogeneous (per fluid volume) and a_v for heterogeneous (per solid surface area) reactions

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, src_id)?;
        check_dom_scldom(vars, dom_id, wgt_id)?;

        // create struct
        let mut oper = PoreOpSclDomSource::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;
        oper.wgt_id = wgt_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclDomSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf w * R dA
        // (w * R)_c * dA; centroid skewness correction as in OpSclDomSource with w * R in place of R

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
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
