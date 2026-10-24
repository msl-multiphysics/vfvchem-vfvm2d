use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclItfAdvection {
    // interface
    pub itf_id: usize,

    // scalars - rows
    pub unk_b_id: usize,  // unknown scalar (c_i^b)

    // scalars - this side
    pub unk_id: usize,           // unknown scalar (c_i); unknown
    pub diff_id: Option<usize>,  // diffusion coefficient (G_ii); for the Peclet number
    pub epsl_id: usize,          // porosity (epsl)
    pub tort_id: usize,          // tortuosity (tau)
    pub den_id: usize,           // density (rho)
    pub displ_id: usize,         // longitudinal dispersivity (alpha_L)
    pub dispt_id: usize,         // transverse dispersivity (alpha_T)

    // vectors - this side
    pub vel_id: usize,  // superficial velocity (vel)

    // flow - this side
    pub flw_id: usize,  // flow data; for the Rhie-Chow mass flux
}

impl PoreOpSclItfAdvection {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, unk_b_id: usize, diff_id: Option<usize>, epsl_id: usize, tort_id: usize, den_id: usize, vel_id: usize, displ_id: usize, dispt_id: usize, flw_id: usize) -> Result<PoreOpSclItfAdvection, VFVChemError> {
        // adds the advection term of a porous medium in terms of the superficial mass flux to the combined control volumes of interface nodes
        // epsl * rho * dc_i/dt + rho * vel . grad(c_i) = -div( -G_eff,ij * grad(c_j) - rho * D_disp * grad(c_i) ) + epsl * R_v,i + a_v * R_s,i
        //                        ^^^^^^^^^^^^^^^^^^^^^
        // covers the half of the control volume on the side of the flow data; loads into the rows of domain b
        // use once for each side, together with PoreOpSclItfDiffusion and PoreOpSclItfDispersion; the power law scheme corrects their central difference
        // the Peclet number uses G_eff,ii + rho * D_nn; a missing diffusion coefficient drops G_eff,ii
        // the boundary half-edges are upwinded from the node and contribute nothing

        // input checking
        // the unknown on this side is either unk_a_id or unk_b_id
        let dom_id = vars.sclflw[flw_id].dom_id;
        let unk_id = if vars.scldom[unk_a_id].dom_id == dom_id { unk_a_id } else { unk_b_id };
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_dom_scldom(vars, dom_id, unk_id)?;
        if let Some(diff_id) = diff_id { check_dom_scldom(vars, dom_id, diff_id)?; }
        check_dom_scldom(vars, dom_id, epsl_id)?;
        check_dom_scldom(vars, dom_id, tort_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, displ_id)?;
        check_dom_scldom(vars, dom_id, dispt_id)?;

        // create struct
        let mut oper = PoreOpSclItfAdvection::default();
        oper.itf_id = itf_id;
        oper.unk_b_id = unk_b_id;
        oper.unk_id = unk_id;
        oper.diff_id = diff_id;
        oper.epsl_id = epsl_id;
        oper.tort_id = tort_id;
        oper.den_id = den_id;
        oper.vel_id = vel_id;
        oper.displ_id = displ_id;
        oper.dispt_id = dispt_id;
        oper.flw_id = flw_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclItfAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -rho * vel . grad(c) dA
        // sum_ngh -F_f * (c_f - c_v) with the superficial mass flux F = rho * vel * dl
        //
        // apply power law scheme
        // sum_ngh (G_f * (A(|Pe_f|) - 1) + max(-F_f, 0)) * (c_n - c_v), G_f = G_eff,f + rho_f * D_nn,f
        //
        // conservative form minus c_v times the discrete continuity residual of the same F_f; see OpSclDomAdvectionUnity

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_id = vars.scldom[self.unk_id].dom_id;
        let dom = &vars.dom[dom_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];
        let tri = &vars.tri[dom_id];
        let den_dom = &vars.scldom[self.den_id];
        let flow_dom = &vars.sclflw[self.flw_id];

        // iterate through nodes
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
            let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];
            let flux_n = &flow_dom.flux[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (nx, ny) = norm_n[loc];
                let (sx, sy) = skew_n[loc];

                // get face effective diffusion coefficient; zero if absent
                let diff_f = match self.diff_id {
                    Some(diff_id) => {
                        let (diff_v, (diff_x_v, diff_y_v)) = pore_helper::diff_eff(vars, diff_id, self.epsl_id, self.tort_id, vid_v);
                        let (diff_n, (diff_x_n, diff_y_n)) = pore_helper::diff_eff(vars, diff_id, self.epsl_id, self.tort_id, vid_n);
                        let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                        let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n));  // skew correction
                        if diff_f > 0.0 { diff_f } else { diff_m }  // guard against negative diffusion
                    }
                    None => 0.0,
                };

                // get face dispersion coefficient; same as PoreOpSclItfDispersion
                let (disp_nn, _) = pore_helper::disp_face(vars, self.vel_id, self.epsl_id, self.displ_id, self.dispt_id, vid_v, vid_n, (nx, ny));
                let den_f = 0.5 * (den_dom.value[vid_v] + den_dom.value[vid_n]);

                // get power law term
                let cond = (diff_f + den_f * disp_nn) * face / dist;
                let flux = flux_n[loc];
                let peclet = if cond != 0.0 { flux / cond } else { 0.0 };
                let pow_law = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff = cond * (pow_law - 1.0) + (-flux).max(0.0);
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_id, vid_v, -coeff)?;
            }

        }

        Ok(())
    }
}
