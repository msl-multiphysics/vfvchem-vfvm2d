use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclDomSolidProduction {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (T_s)
    pub prd_id: usize,   // production coefficient per solid volume (s_sv)
    pub epsl_id: usize,  // porosity (epsl); weight is 1 - epsl
}

impl PoreOpSclDomSolidProduction {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, prd_id: usize, epsl_id: usize) -> Result<PoreOpSclDomSolidProduction, VFVChemError> {
        // adds the solid volumetric production term of a porous medium
        // (1 - epsl) * (rho * c)_s * dT_s/dt = ... + (1 - epsl) * (s_sv * T_s + q_sv) + ...
        //                                            ^^^^^^^^^^^^^^^^^^^^^^^
        // same as PoreOpSclDomProduction with the weight 1 - epsl

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, prd_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclDomSolidProduction::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.prd_id = prd_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclDomSolidProduction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf (1 - epsl) * s * T dA
        // ((1 - epsl) * s)_c * T_c * dA; centroid skewness correction as in OpSclDomProduction

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
            let (prd_v, (prd_x_v, prd_y_v)) = pore_helper::solid_weighted(vars, self.prd_id, self.epsl_id, vid_v);
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
