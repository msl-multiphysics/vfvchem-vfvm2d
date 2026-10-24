use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpVecBndDrag {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // superficial velocity (vel); unknown

    // scalars
    pub den_id: usize,   // density (rho)
    pub visc_id: usize,  // fluid viscosity (mu)
    pub epsl_id: usize,  // porosity (epsl)
    pub diam_id: usize,  // particle diameter (diam)
}

impl PoreOpVecBndDrag {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, den_id: usize, visc_id: usize, epsl_id: usize, diam_id: usize) -> Result<PoreOpVecBndDrag, VFVChemError> {
        // adds the Darcy and Forchheimer drag to the control volumes of boundary nodes
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        //                                                                                                                                               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // K and F are from the Ergun relation; see PoreOpVecDomDrag

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;
        check_bnd_scldom(vars, bnd_id, visc_id)?;
        check_bnd_scldom(vars, bnd_id, epsl_id)?;
        check_bnd_scldom(vars, bnd_id, diam_id)?;

        // create struct
        let mut oper = PoreOpVecBndDrag::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.visc_id = visc_id;
        oper.epsl_id = epsl_id;
        oper.diam_id = diam_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpVecBndDrag {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -(mu / K + rho * F / sqrt(K) * |vel|) * vel dA
        // -(mu_v / K_v + rho_v * F_v / sqrt(K_v) * |vel_v|) * vel_v * dA
        //
        // linearize the Forchheimer term about the previous iterate (vel_k)
        // |vel| * vel = 2 * |vel_k| * vel - |vel_k| * vel_k
        // lagging |vel| (Picard) stalls where the Forchheimer drag dominates; triplets are always added to keep the sparsity pattern fixed

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        let visc_dom = &vars.scldom[self.visc_id];
        let epsl_dom = &vars.scldom[self.epsl_id];
        let diam_dom = &vars.scldom[self.diam_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let area = tri.cvol_area[vid_v];

            // get vertex properties
            let (u_v, v_v) = vel_dom.value[vid_v];
            let speed = (u_v * u_v + v_v * v_v).sqrt();
            let (kinv, fsk) = PoreOpVecDomDrag::closure(epsl_dom.value[vid_v], diam_dom.value[vid_v]);
            let drag_lin = visc_dom.value[vid_v] * kinv;
            let drag_quad = den_dom.value[vid_v] * fsk * speed;

            // load into matrix and vector
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v, -(drag_lin + 2.0 * drag_quad) * area)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v, -(drag_lin + 2.0 * drag_quad) * area)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -drag_quad * u_v * area)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -drag_quad * v_v * area)?;
        }

        Ok(())
    }
}
