use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpSclBndDivergence {
    // boundary
    pub bnd_id: usize,

    // unknowns
    pub vel_id: usize,  // velocity (vel)

    // others
    pub den_id: usize,  // density (rho)
    pub pres_id: usize,  // pressure (p)
}

impl OpSclBndDivergence {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, den_id: usize, pres_id: usize) -> Result<OpSclBndDivergence, VFVChemError> {
        // adds the divergence term to the continuity equation
        // d(rho)/dt = -div(rho * vel)
        //             ^^^^^^^^^^^^^^^
        
        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;
        check_bnd_scldom(vars, bnd_id, pres_id)?;

        // create struct
        let mut oper = OpSclBndDivergence::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.pres_id = pres_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_surf -div(rho * vel) dA
        //
        // apply divergence theorem
        // int_lin -rho * vel . n dl
        //
        // subset the boundary term
        // int_bnd -rho * vel . n dl
        //
        // discretize flux using linear interpolation
        // rho * vel = rho_v * vel_v * (1 - r/l) + rho_n * vel_n * r/l
        //
        // integrate over the half-edge owned by the vertex
        // int_0^(l/2) -rho * vel . n dl
        // -(0.75 * rho_v * vel_v + 0.25 * rho_n * vel_n) . n * half

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let den_dom = &vars.scldom[self.den_id];

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
            let den_v = den_dom.value[vid_v];

            // iterate through boundary half-edges
            for (blc, &vid_n) in list_b.iter().enumerate() {
                // get geometry
                let half = half_b[blc];
                let (nx, ny) = norm_b[blc];

                // get neighbor properties
                let den_n = den_dom.value[vid_n];

                // load into matrix
                let coeff_v = -0.75 * den_v * half;
                let coeff_n = -0.25 * den_n * half;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_v, coeff_v * nx)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 0, vid_n, coeff_n * nx)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_v, coeff_v * ny)?;
                self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, 1, vid_n, coeff_n * ny)?;
            }

        }

        Ok(())
    }
}
