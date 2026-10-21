use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecBndPressureTraction {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub pres_id: usize,  // pressure (p); unknown
    pub val_id: usize,  // boundary pressure (p_set)
}

impl OpVecBndPressureTraction {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, pres_id: usize, val_id: usize) -> Result<OpVecBndPressureTraction, VFVChemError> {
        // prescribes the pressure on the boundary half-edges of the control volumes of boundary nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                          ^^^^^^^^
        // p_b = p_set; must be used with OpVecBndGradient

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, pres_id)?;
        check_bnd_sclbnd(vars, bnd_id, val_id)?;

        // create struct
        let mut oper = OpVecBndPressureTraction::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;
        oper.val_id = val_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecBndPressureTraction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_bnd -p n dl
        // sum_bnd -p_b * n * dl
        //
        // remove the boundary half-edges of OpVecBndGradient
        // sum_bnd (0.75 * p_v + 0.25 * p_n) * n * dl
        //
        // replace with the prescribed pressure
        // sum_bnd -(0.75 * p_set_v + 0.25 * p_set_n) * n * dl

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let val_bnd = &vars.sclbnd[self.val_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_n = &tri.bnd_vid[vid_v];
            let half_n = &tri.bnd_half[vid_v];
            let norm_n = &tri.bnd_norm[vid_v];

            // get vertex properties
            let val_v = val_bnd.value[idx_v];

            // iterate through boundary half-edges
            for (blc, &vid_n) in list_n.iter().enumerate() {
                // get indices
                let nid_n = vid_n - dom.num_pts();
                let gnd_n = dom.nds_gnd[nid_n];

                // get geometry
                let half = half_n[blc];
                let (nx, ny) = norm_n[blc];

                // get neighbor properties
                // use vertex pressure if neighbor is on another boundary
                let val_n = match bnd.nds_idx.get(&gnd_n) {
                    Some(&idx_n) => val_bnd.value[idx_n],
                    None => val_v,
                };

                // load into matrix
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_v, 0.75 * half * nx)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_v, 0.75 * half * ny)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_n, 0.25 * half * nx)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_n, 0.25 * half * ny)?;

                // load into vector
                let val_b = 0.75 * val_v + 0.25 * val_n;
                self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, val_b * half * nx)?;
                self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, val_b * half * ny)?;
            }
        }

        Ok(())
    }
}
