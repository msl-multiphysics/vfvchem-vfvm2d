use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct FlowBndTranspose {
    // boundary
    pub bnd_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub visc_id: usize,  // viscosity (mu)
}

impl FlowBndTranspose {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, visc_id: usize) -> Result<FlowBndTranspose, VFVChemError> {
        // adds the boundary term of the transpose term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                      ^^^^^^^^^               ^^^^^^^^^^^
        //
        // used on pressure boundaries where the momentum equation is solved
        // together with the natural condition mu * grad(vel) . n = 0, this gives a do-nothing outlet

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, visc_id)?;

        // create struct
        let mut oper = FlowBndTranspose::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for FlowBndTranspose {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf div( mu * grad(vel)^T ) dA
        //
        // apply divergence theorem
        // int_lin mu * grad(vel)^T . n dl
        //
        // subset the boundary term
        // int_bnd mu * grad(vel)^T . n dl
        //
        // discretize using linear interpolation of vertex values
        // mu * grad(vel)^T = (mu * grad(vel)^T)_v * (1 - r/l) + (mu * grad(vel)^T)_n * r/l
        //
        // integrate over the half-edge owned by the vertex
        // int_0^(l/2) mu * grad(vel)^T . n dl
        // (0.75 * (mu * grad(vel)^T)_v + 0.25 * (mu * grad(vel)^T)_n) . n * half

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let visc_dom = &vars.scldom[self.visc_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_b = &tri.bnd_vid[vid_v];
            let half_b = &tri.bnd_half[vid_v];
            let norm_b = &tri.bnd_norm[vid_v];

            // get vertex properties
            let (velx_x_v, velx_y_v, vely_x_v, vely_y_v) = vel_dom.grad[vid_v];
            let visc_v = visc_dom.value[vid_v];

            // iterate through boundary half-edges
            let mut sumx = 0.0;
            let mut sumy = 0.0;
            for (blc, &vid_n) in list_b.iter().enumerate() {
                // get geometry
                let half = half_b[blc];
                let (nx, ny) = norm_b[blc];

                // get neighbor properties
                let (velx_x_n, velx_y_n, vely_x_n, vely_y_n) = vel_dom.grad[vid_n];
                let visc_n = visc_dom.value[vid_n];

                // add to sum
                let trnx_v = visc_v * (velx_x_v * nx + vely_x_v * ny);
                let trny_v = visc_v * (velx_y_v * nx + vely_y_v * ny);
                let trnx_n = visc_n * (velx_x_n * nx + vely_x_n * ny);
                let trny_n = visc_n * (velx_y_n * nx + vely_y_n * ny);
                sumx += (0.75 * trnx_v + 0.25 * trnx_n) * half;
                sumy += (0.75 * trny_v + 0.25 * trny_n) * half;
            }

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -sumx)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -sumy)?;
        }

        Ok(())
    }
}
