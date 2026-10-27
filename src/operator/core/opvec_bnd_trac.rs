use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpVecBndTraction {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub pres_id: usize,  // pressure (p); unknown
    pub val_id: usize,   // outflow pressure (p_set)
}

impl OpVecBndTraction {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, pres_id: usize, val_id: usize) -> Result<OpVecBndTraction, VFVChemError> {
        // adds the pressure gradient term to the control volumes of nodes on a traction outflow boundary
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                          ^^^^^^^^
        // sets the traction on the boundary half-edges
        // (-p * I + mu * grad(vel)) . n = -p_set * n
        // the viscous boundary operators omit the boundary half-edges, so the viscous flux there is zero
        // replaces OpSclBndDirichlet on pressure and OpVecBndGradient; pressure stays unknown on the boundary
        // and the continuity equation of the boundary nodes is kept

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, pres_id)?;
        check_bnd_sclbnd(vars, bnd_id, val_id)?;

        // create struct
        let mut oper = OpVecBndTraction::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;
        oper.val_id = val_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecBndTraction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // int_surf -grad(p) dA
        // int_lin -p n dl
        // sum_ngh -p_f * n * dl + sum_bnd -p_b * n * dl
        // sum_ngh -p_f * n * dl + sum_bnd -(0.75 * p_v + 0.25 * p_n) * n * dl
        //
        // expand the least-squares gradients; see Triangle
        // sum_m -(c_x, c_y) * p_m
        //
        // on half-edges along this boundary, replace the interpolated pressure with the outflow pressure
        // sum_bnd (0.75 * p_v + 0.25 * p_n) * n * dl - (0.75 * p_set_v + 0.25 * p_set_n) * n * dl

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let val_bnd = &vars.sclbnd[self.val_id];

        // get pressure gradient coefficients
        let grad_vec = &tri.grad_coeff;

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

            // load gradient into matrix
            for &(vid_m, c_x, c_y) in &grad_vec[vid_v] {
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_m, -c_x)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_m, -c_y)?;
            }

            // iterate through boundary half-edges
            for (blc, &vid_n) in list_n.iter().enumerate() {
                // get indices
                let nid_n = vid_n - dom.num_pts();
                let gnd_n = dom.nds_gnd[nid_n];

                // skip if neighbor is on another boundary
                let Some(&idx_n) = bnd.nds_idx.get(&gnd_n) else { continue };

                // get geometry
                let half = half_n[blc];
                let (nx, ny) = norm_n[blc];

                // get neighbor properties
                let val_n = val_bnd.value[idx_n];

                // load into matrix and vector
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_v, 0.75 * half * nx)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_v, 0.75 * half * ny)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_n, 0.25 * half * nx)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_n, 0.25 * half * ny)?;
                self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, (0.75 * val_v + 0.25 * val_n) * half * nx)?;
                self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, (0.75 * val_v + 0.25 * val_n) * half * ny)?;
            }
        }

        Ok(())
    }
}
