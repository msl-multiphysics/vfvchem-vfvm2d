use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpVecBndSource {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown
    pub fce_id: usize,  // body force (fce)
}

impl OpVecBndSource {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, fce_id: usize) -> Result<OpVecBndSource, VFVChemError> {
        // adds the body force term to the control volumes of boundary nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                                                                                    ^^^

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_vecdom(vars, bnd_id, fce_id)?;

        // create struct
        let mut oper = OpVecBndSource::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.fce_id = fce_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecBndSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf fce dA
        // fce_c * dA
        //
        // apply centroid skewness correction
        // fce_c = fce_v + grad(fce_v) . skew

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let fce_dom = &vars.vecdom[self.fce_id];

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
            let (fce_x_v, fce_y_v) = fce_dom.value[vid_v];
            let (fce_x_x_v, fce_x_y_v, fce_y_x_v, fce_y_y_v) = fce_dom.grad[vid_v];

            // get body force at centroid
            let fce_x_c = fce_x_v + fce_x_x_v * sx + fce_x_y_v * sy;
            let fce_y_c = fce_y_v + fce_y_x_v * sx + fce_y_y_v * sy;

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -fce_x_c * area)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -fce_y_c * area)?;
        }

        Ok(())
    }
}
