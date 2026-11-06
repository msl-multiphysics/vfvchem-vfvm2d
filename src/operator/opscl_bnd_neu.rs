use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct OpSclBndNeumann {
    // boundary
    pub bnd_id: usize,

    // unknowns
    pub unk_id: usize,   // unknown scalar (c_i)

    // others
    pub flux_id: usize,  // flux (N_set)
}

impl OpSclBndNeumann {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, flux_id: usize) -> Result<OpSclBndNeumann, VFVChemError> {
        // sets the outward boundary flux
        // N_i . n = N_set

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_sclbnd(vars, bnd_id, flux_id)?;
        
        // create struct
        let mut oper = OpSclBndNeumann::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.flux_id = flux_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndNeumann {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // get surface integral
        // int_bnd N . n dl
        //
        // discretize flux using linear interpolation
        // N = N_v * (1 - r/l) + N_n * r/l
        //
        // integrate over the half-edge owned by the vertex
        // int_0^(l/2) N dl
        // (0.75 * N_v + 0.25 * N_n) . n * half
        
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let flux_bnd = &vars.sclbnd[self.flux_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let half_b = &tri.bnd_half[vid_v];
            let list_b = &tri.bnd_vid[vid_v];
            
            // get vertex properties
            let flux_v = flux_bnd.value[idx_v];

            // iterate through boundary half-edges
            for (blc, &vid_n) in list_b.iter().enumerate() {
                // get indices
                // idx_n not found -> neighbor on different boundary
                let nid_n = vid_n - dom.num_pts();
                let gnd_n = dom.nds_gnd[nid_n];
                let idx_n = match bnd.nds_idx.get(&gnd_n) {
                    Some(&idx_n) => idx_n,
                    None => idx_v,  // use idx_v properties if idx_n is a different boundary
                };

                // get geometry
                let half = half_b[blc];

                // get neighbor properties
                let flux_n = flux_bnd.value[idx_n];

                // load into matrix
                let coeff_v = (0.75 * flux_v + 0.25 * flux_n) * half;
                self.add_scl_b(vars, b_vec, self.unk_id, vid_v, coeff_v)?;
            }
        }

        Ok(())
    }
}
