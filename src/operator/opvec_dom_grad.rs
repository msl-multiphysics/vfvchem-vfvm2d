use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpVecDomGradient {
    // domain
    pub dom_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub pres_id: usize,  // pressure (p)
}

impl OpVecDomGradient {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, pres_id: usize) -> Result<OpVecDomGradient, VFVChemError> {
        // adds the advection term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                           ^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, pres_id)?;

        // create struct
        let mut oper = OpVecDomGradient::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecDomGradient {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {  
        // get surface integral
        // int_surf -grad(p) dA
        //
        // discretize into centroid values
        // -grad(p) * dA
        //
        // notes:
        // gradients are discretized implicitly using gradient weights
        // divergence theorem is not used to avoid explicit corrections

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let area = tri.vert_area[vid_v];
            let list_n = &tri.ngh_vid[vid_v];
            let wgt_n = &tri.grad_wgt[vid_v];

            // iterate through neighbors
            for (loc, &vid_n) in list_n.iter().enumerate() {
                let (w_x, w_y) = wgt_n[loc];
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_n, -w_x * area)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_v,  w_x * area)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_n, -w_y * area)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_v,  w_y * area)?;
            }
        }

        Ok(())
    }
}
