use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclBndSolidSource {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (T_s)
    pub src_id: usize,   // source per solid volume (q_sv)
    pub epsl_id: usize,  // porosity (epsl); weight is 1 - epsl
}

impl PoreOpSclBndSolidSource {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, src_id: usize, epsl_id: usize) -> Result<PoreOpSclBndSolidSource, VFVChemError> {
        // adds the solid volumetric source term of a porous medium
        // (1 - epsl) * (rho * c)_s * dT_s/dt = ... + (1 - epsl) * (s_sv * T_s + q_sv) + ...
        //                                                                    ^^^^
        // same as PoreOpSclBndSource with the weight 1 - epsl

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, src_id)?;
        check_bnd_scldom(vars, bnd_id, epsl_id)?;

        // create struct
        let mut oper = PoreOpSclBndSolidSource::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;
        oper.epsl_id = epsl_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclBndSolidSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf (1 - epsl) * q dA
        // ((1 - epsl) * q)_c * dA; centroid skewness correction as in OpSclDomSource

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
            let (src_v, (src_x_v, src_y_v)) = pore_helper::solid_weighted(vars, self.src_id, self.epsl_id, vid_v);

            // get source at centroid
            let src_c = src_v + src_x_v * sx + src_y_v * sy;

            // load into vector
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -src_c * area)?;
        }

        Ok(())
    }
}
