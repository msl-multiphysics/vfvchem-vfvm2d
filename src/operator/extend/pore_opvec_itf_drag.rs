use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpVecItfDrag {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,    // superficial velocity (vel^a); unknown
    pub den_a_id: usize,    // density (rho^a)
    pub visc_a_id: usize,   // fluid viscosity (mu^a)
    pub perm_a_id: usize,   // permeability (K^a)
    pub forch_a_id: usize,  // Forchheimer coefficient (F^a)

    // domain b
    pub vel_b_id: usize,    // superficial velocity (vel^b); unknown
    pub den_b_id: usize,    // density (rho^b)
    pub visc_b_id: usize,   // fluid viscosity (mu^b)
    pub perm_b_id: usize,   // permeability (K^b)
    pub forch_b_id: usize,  // Forchheimer coefficient (F^b)
}

impl PoreOpVecItfDrag {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, den_a_id: usize, visc_a_id: usize, perm_a_id: usize, forch_a_id: usize, vel_b_id: usize, den_b_id: usize, visc_b_id: usize, perm_b_id: usize, forch_b_id: usize) -> Result<PoreOpVecItfDrag, VFVChemError> {
        // adds the Darcy and Forchheimer drag to the combined control volumes of interface nodes
        // (rho / epsl) * d(vel)/dt + (rho / epsl) * vel . grad(vel / epsl) = -grad(p) + div( mu_eff * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) - (mu / K) * vel - (rho * F / sqrt(K)) * |vel| * vel + fce
        //                                                                                                                                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // K and F are user inputs; see PoreOpVecDomDrag
        // each domain contributes its half of the control volume; loads into the rows of domain b

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, den_a_id, den_b_id)?;
        check_itf_scldom(vars, itf_id, visc_a_id, visc_b_id)?;
        check_itf_scldom(vars, itf_id, perm_a_id, perm_b_id)?;
        check_itf_scldom(vars, itf_id, forch_a_id, forch_b_id)?;

        // create struct
        let mut oper = PoreOpVecItfDrag::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.den_a_id = den_a_id;
        oper.visc_a_id = visc_a_id;
        oper.perm_a_id = perm_a_id;
        oper.forch_a_id = forch_a_id;
        oper.vel_b_id = vel_b_id;
        oper.den_b_id = den_b_id;
        oper.visc_b_id = visc_b_id;
        oper.perm_b_id = perm_b_id;
        oper.forch_b_id = forch_b_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpVecItfDrag {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -(mu / K + rho * F / sqrt(K) * |vel|) * vel dA
        // -(mu_v / K_v + rho_v * F_v / sqrt(K_v) * |vel_v|) * vel_v * dA
        //
        // linearize the Forchheimer term about the previous iterate (vel_k)
        // |vel| * vel = 2 * |vel_k| * vel - |vel_k| * vel_k
        // triplets are always added to keep the sparsity pattern fixed

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, den_id, visc_id, perm_id, forch_id) in [(self.vel_a_id, self.den_a_id, self.visc_a_id, self.perm_a_id, self.forch_a_id), (self.vel_b_id, self.den_b_id, self.visc_b_id, self.perm_b_id, self.forch_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tri = &vars.tri[dom_id];
            let vel_dom = &vars.vecdom[vel_id];
            let den_dom = &vars.scldom[den_id];
            let visc_dom = &vars.scldom[visc_id];
            let perm_dom = &vars.scldom[perm_id];
            let forch_dom = &vars.scldom[forch_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let area = tri.cvol_area[vid_v];

                // get vertex properties
                let (u_v, v_v) = vel_dom.value[vid_v];
                let speed = (u_v * u_v + v_v * v_v).sqrt();
                if !(perm_dom.value[vid_v] > 0.0) { return Err(VFVChemError::PoreNonPositivePermeability { dom_id }) }
                let (kinv, fsk) = PoreOpVecDomDrag::coeff(perm_dom.value[vid_v], forch_dom.value[vid_v]);
                let drag_lin = visc_dom.value[vid_v] * kinv;
                let drag_quad = den_dom.value[vid_v] * fsk * speed;

                // load into matrix and vector
                self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 0, vid_b_v, vel_id, 0, vid_v, -(drag_lin + 2.0 * drag_quad) * area)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 1, vid_b_v, vel_id, 1, vid_v, -(drag_lin + 2.0 * drag_quad) * area)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 0, vid_b_v, -drag_quad * u_v * area)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 1, vid_b_v, -drag_quad * v_v * area)?;
            }
        }

        Ok(())
    }
}
