use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct PoreOpSclDomSolidSource {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (T_s)
    pub src_id: usize,   // source per solid volume (q_sv)
    pub epsl_id: usize,  // porosity (epsl); weight is 1 - epsl
}

impl PoreOpSclDomSolidSource {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, src_id: usize, epsl_id: usize) -> Result<PoreOpSclDomSolidSource, VFVChemError> {
        // adds the solid volumetric source term of a porous medium
        // (1 - epsl) * (rho * c)_s * dT_s/dt = ... + (1 - epsl) * (s_sv * T_s + q_sv) + ...
        //                                                                    ^^^^
        // same as PoreOpSclDomSource with the weight 1 - epsl

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, src_id)?;
        check_dom_scldom(vars, dom_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclDomSolidSource::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclDomSolidSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf (1 - epsl) * q dA
        // ((1 - epsl) * q)_c * dA; centroid skewness correction as in OpSclDomSource

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
            let (src_v, (src_x_v, src_y_v)) = pore_helper::solid_weighted(vars, self.src_id, self.epsl_id, vid_v);

            // get source at centroid
            let src_c = src_v + src_x_v * sx + src_y_v * sy;

            // load into vector
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -src_c * area)?;
        }

        Ok(())
    }
}
