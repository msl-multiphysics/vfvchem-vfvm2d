use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpVecDomSource {
    // domain
    pub dom_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub fce_id: usize,  // body force (fce)
}

impl OpVecDomSource {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, fce_id: usize) -> Result<OpVecDomSource, VFVChemError> {
        // adds the advection term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                                                                                     ^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_vecdom(vars, dom_id, fce_id)?;

        // create struct
        let mut oper = OpVecDomSource::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.fce_id = fce_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecDomSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf fce dA
        //
        // discretize into centroid values
        // fce_c * dA

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let fce_dom = &vars.vecdom[self.fce_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let offc = tri.vert_offc[vid_v];
            
            // get centroid properties
            let (fcex_c, fcey_c) = fce_dom.interp_vert(offc, vid_v);

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -fcex_c * area)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -fcey_c * area)?;
        }

        Ok(())
    }
}
