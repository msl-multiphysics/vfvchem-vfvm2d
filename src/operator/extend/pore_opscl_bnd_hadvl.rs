use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct PoreOpSclBndHeatAdvection {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,    // unknown scalar (T_f)
    pub wgt_id: usize,    // fluid volumetric heat capacity (rho * cp)_f
    pub cond_id: usize,   // fluid thermal conductivity (k_f); for the Peclet number
    pub epsl_id: usize,   // porosity (epsl)
    pub den_id: usize,    // density (rho)
    pub displ_id: usize,  // longitudinal thermal dispersivity (alpha_L)
    pub dispt_id: usize,  // transverse thermal dispersivity (alpha_T)

    // vectors
    pub vel_id: usize,  // superficial velocity (vel)

    // flow
    pub flw_id: usize,  // flow data; for the Rhie-Chow mass flux
}

impl PoreOpSclBndHeatAdvection {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, wgt_id: usize, cond_id: usize, epsl_id: usize, den_id: usize, vel_id: usize, displ_id: usize, dispt_id: usize, flw_id: usize) -> Result<PoreOpSclBndHeatAdvection, VFVChemError> {
        // adds the fluid heat advection term of a porous medium in terms of the superficial mass flux
        // epsl * (rho * cp)_f * dT_f/dt + (rho * cp)_f * vel . grad(T_f) = -div( -epsl * k_f * grad(T_f) - (rho * cp)_f * D_disp * grad(T_f) ) + ...
        //                                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // same as OpSclDomAdvection with the weight (rho * cp)_f; must be used with PoreOpSclBndHeatDiffusion and PoreOpSclBndDispersion
        // the Peclet number uses epsl * k_f + (rho * cp)_f * D_nn; reduces to OpSclDomAdvection for epsl = 1 and alpha_L = alpha_T = 0

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, wgt_id)?;
        check_bnd_scldom(vars, bnd_id, cond_id)?;
        check_bnd_scldom(vars, bnd_id, epsl_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, displ_id)?;
        check_bnd_scldom(vars, bnd_id, dispt_id)?;
        check_bnd_sclflw(vars, bnd_id, flw_id)?;

        // create struct
        let mut oper = PoreOpSclBndHeatAdvection::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.wgt_id = wgt_id;
        oper.cond_id = cond_id;
        oper.epsl_id = epsl_id;
        oper.den_id = den_id;
        oper.vel_id = vel_id;
        oper.displ_id = displ_id;
        oper.dispt_id = dispt_id;
        oper.flw_id = flw_id;

        Ok(oper)
    }
}

impl OperatorBase for PoreOpSclBndHeatAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_surf -(rho * cp)_f * vel . grad(T_f) dA
        // sum_ngh -F_f * (T_f - T_v) with the weight flux F = ((rho * cp)_f / rho) * rho * vel * dl
        //
        // apply power law scheme
        // sum_ngh (G_f * (A(|Pe_f|) - 1) + max(-F_f, 0)) * (T_n - T_v), G_f = epsl_f * k_f + (rho * cp)_f * D_nn,f

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tri = &vars.tri[bnd.dom_id];
        let wgt_dom = &vars.scldom[self.wgt_id];
        let den_dom = &vars.scldom[self.den_id];
        let flow_dom = &vars.sclflw[self.flw_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_n = &tri.ngh_vid[vid_v];
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];
            let flux_n = &flow_dom.flux[vid_v];

            // get vertex properties
            let wgt_v = wgt_dom.value[vid_v];
            let (wgt_x_v, wgt_y_v) = wgt_dom.grad[vid_v];
            let (diff_v, (diff_x_v, diff_y_v)) = pore_helper::weighted(vars, self.cond_id, self.epsl_id, vid_v);
            let den_v = den_dom.value[vid_v];
            let (den_x_v, den_y_v) = den_dom.grad[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (nx, ny) = norm_n[loc];
                let (sx, sy) = skew_n[loc];

                // get neighbor properties
                let wgt_n = wgt_dom.value[vid_n];
                let (wgt_x_n, wgt_y_n) = wgt_dom.grad[vid_n];
                let (diff_n, (diff_x_n, diff_y_n)) = pore_helper::weighted(vars, self.cond_id, self.epsl_id, vid_n);
                let den_n = den_dom.value[vid_n];
                let (den_x_n, den_y_n) = den_dom.grad[vid_n];

                // get face weight, density, and effective conductivity
                let wgt_f = 0.5 * (wgt_v + wgt_n) + 0.5 * (sx * (wgt_x_v + wgt_x_n) + sy * (wgt_y_v + wgt_y_n));  // skew-corrected arithmetic mean
                let den_f = 0.5 * (den_v + den_n) + 0.5 * (sx * (den_x_v + den_x_n) + sy * (den_y_v + den_y_n));  // skew-corrected arithmetic mean
                let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n));  // skew correction
                let diff_f = if diff_f > 0.0 { diff_f } else { diff_m };  // guard against negative diffusion

                // get face dispersion coefficient; same as PoreOpSclBndDispersion with (rho * cp)_f in place of rho
                let (disp_nn, _) = pore_helper::disp_face(vars, self.vel_id, self.epsl_id, self.displ_id, self.dispt_id, vid_v, vid_n, (nx, ny));
                let wgt_m = 0.5 * (wgt_v + wgt_n);

                // get power law term
                let cond = (diff_f + wgt_m * disp_nn) * face / dist;
                let flux = if den_f != 0.0 { wgt_f / den_f * flux_n[loc] } else { 0.0 };
                let peclet = if cond != 0.0 { flux / cond } else { 0.0 };
                let pow_law = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff = cond * (pow_law - 1.0) + (-flux).max(0.0);
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -coeff)?;
            }
        }

        Ok(())
    }
}
