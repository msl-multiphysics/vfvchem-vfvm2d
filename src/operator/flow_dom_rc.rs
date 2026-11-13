use crate::base::prelude::*;
use crate::operator::prelude::*;

#[derive(Default)]
pub struct FlowDomRhieChow {
    // domain
    pub dom_id: usize,

    // unknowns
    pub pres_id: usize,  // pressure (p)
    
    // others
    pub vel_id: usize,  // velocity (vel)
    pub den_id: usize,   // density (rho)
}

impl FlowDomRhieChow {
    pub fn new(vars: &Variables, dom_id: usize, pres_id: usize, vel_id: usize, den_id: usize) -> Result<FlowDomRhieChow, VFVChemError> {
        // adds the rhie-chow correction to the continuity equation face velocities
        // must be applied before dirichlet boundaries and pressure point

        // input checking
        check_dom_scldom(vars, dom_id, pres_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = FlowDomRhieChow::default();
        oper.dom_id = dom_id;
        oper.pres_id = pres_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>, twgt: (f64, f64, f64)) -> Result<(), VFVChemError> {
        // adds the rhie-chow correction for transient equations; used instead of apply
        // must be applied before dirichlet boundaries and pressure point
        //
        // start with discretized equations
        // expand a into steady and transient components
        // a vel = mom - b grad(p)
        // a = a_s + a_t
        // a_t = rho_c * dA * tw0
        //
        // the `apply` function below uses a (with a_s + a_t)
        // this uses the combined d
        // d = b / a = b / (a_s + a_t) -> (eq. 1)
        // 
        // define d_s that uses only a_s and removes a_t
        // d_s = b / a_s -> (eq. 2)
        // 
        // express d_s in terms of d and a_t
        // a_s = a - a_t = b / d - a_t
        // d_s = b / (b / d - a_t)
        // 
        // use d_s instead of d in rhie-chow
        // rho_f * (d_s grad(p) . n)_f * dl - rho_f * d_sf (p_n - p_v) / dn * dl
        //
        // notes:
        // this is the same code as in `apply`, but using d_s = b / (b / d - a_t)
        // this is timestep-independent, but not the complete rhie-chow interpolation.

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        let (tw0, _, _) = twgt;
        
        // get steady rhie-chow coefficients
        let rc_vec = _get_rc(dom, tri, vel_dom, a_triplet)?;
        let mut rcs_vec = Vec::with_capacity(dom.num_vert());
        for vid in 0..dom.num_vert() {
            let area = tri.vert_area[vid];
            let den_c = den_dom.interp_vert(tri.vert_offc[vid], vid);
            let coeff_t = -den_c * area * tw0;  // a_t
            let rc = rc_vec[vid];
            let rcs = area / (area / rc - coeff_t);
            rcs_vec.push(rcs);
        }

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let list_n = &tri.ngh_vid[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // get vertex properties
            let rcs_v = rcs_vec[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let dist = dist_n[loc];
                let face = face_n[loc];
                let (nx, ny) = norm_n[loc];
                let skew = skew_n[loc];

                // get neighbor properties
                let rcs_n = rcs_vec[vid_n];

                // get face properties
                let den_f = den_dom.interp_face(skew, vid_v, vid_n);
                let rcs_f = 0.5 * (rcs_v + rcs_n);

                // load into matrix
                let coeff_v = den_f * rcs_f * face / dist;  // -rho_f * d_f (p_n - p_v) / dn * dl
                self.add_sclscl_a(vars, a_triplet, self.pres_id, vid_v, self.pres_id, vid_n, -coeff_v)?;
                self.add_sclscl_a(vars, a_triplet, self.pres_id, vid_v, self.pres_id, vid_v,  coeff_v)?;

                // load gradient term into matrix; implicit to avoid lagging behind the difference term
                // (d grad(p))_f = 0.5 * (d_v grad(p)_v + d_n grad(p)_n); grad(p)_w = sum_m wgt_m (p_m - p_w)
                for vid_w in [vid_v, vid_n] {
                    let rcs_w = rcs_vec[vid_w];
                    for (loc_w, &vid_m) in tri.ngh_vid[vid_w].iter().enumerate() {
                        let (wgt_x, wgt_y) = tri.grad_wgt[vid_w][loc_w];
                        let coeff_g = 0.5 * den_f * rcs_w * (wgt_x * nx + wgt_y * ny) * face;  // rho_f * (d grad(p) . n)_f * dl
                        self.add_sclscl_a(vars, a_triplet, self.pres_id, vid_v, self.pres_id, vid_m,  coeff_g)?;
                        self.add_sclscl_a(vars, a_triplet, self.pres_id, vid_v, self.pres_id, vid_w, -coeff_g)?;
                    }
                }
            }
        }

        Ok(())
    }
}

impl OperatorBase for FlowDomRhieChow {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, _b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // start with discretized equations
        // a and b are coefficients of vel_v and grad(p); mom are the other terms
        // a vel = mom - b grad(p)
        //
        // divide by a_v
        // vel = (mom / a) - (b / a) grad(p)
        // vel = vel' - d grad(p)
        //
        // interpolate velocity at face
        // vel_f = vel_f' - (d grad(p))_f -> (eq. 1)
        //
        // define rhie-chow face velocity (vel_f")
        // similar to eq. 1 but grad(p) is replaced with difference
        // (vel" . n)_f = (vel' . n)_f - d_f (p_n - p_v) / dn -> (eq. 2)
        //
        // eliminatate vel_f' from equations
        // (vel" . n)_f = (vel_f + (d grad(p))_f) . n - d_f (p_n - p_v) / dn
        // (vel" . n)_f = (vel . n)_f + (d grad(p) . n)_f - d_f (p_n - p_v) / dn
        //
        // multiply by rho_f * dl for consistency with divergence operator
        // rho_f * (vel" . n)_f * dl = rho_f * (vel . n)_f * dl + rho_f * (d grad(p) . n)_f * dl - rho_f * d_f (p_n - p_v) / dn * dl
        //
        // rhie-chow replaces face velocity in divergence with vel_f"
        // divergence term already has vel_f; add the pressure terms only
        // rho_f * (d grad(p) . n)_f * dl - rho_f * d_f (p_n - p_v) / dn * dl
        //
        // notes:
        // d is scalar because the x and y momentum equations have the same coefficients
        // grad(p) is implicit; a lagged grad(p) does not cancel the implicit difference term, and the iteration diverges
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let tri = &vars.tri[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];
        
        // get rhie-chow coefficients
        let rc_vec = _get_rc(dom, tri, vel_dom, a_triplet)?;

        // iterate through points
        for vid_v in 0..dom.num_vert() {
            // get geometry
            let dist_n = &tri.ngh_dist[vid_v];
            let face_n = &tri.face_len[vid_v];
            let list_n = &tri.ngh_vid[vid_v];
            let norm_n = &tri.face_norm[vid_v];
            let skew_n = &tri.face_skew[vid_v];

            // get vertex properties
            let rc_v = rc_vec[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let dist = dist_n[loc];
                let face = face_n[loc];
                let (nx, ny) = norm_n[loc];
                let skew = skew_n[loc];

                // get neighbor properties
                let rc_n = rc_vec[vid_n];

                // get face properties
                let den_f = den_dom.interp_face(skew, vid_v, vid_n);
                let rc_f = 0.5 * (rc_v + rc_n);

                // load into matrix
                let coeff_v = den_f * rc_f * face / dist;  // -rho_f * d_f (p_n - p_v) / dn * dl
                self.add_sclscl_a(vars, a_triplet, self.pres_id, vid_v, self.pres_id, vid_n, -coeff_v)?;
                self.add_sclscl_a(vars, a_triplet, self.pres_id, vid_v, self.pres_id, vid_v,  coeff_v)?;

                // load gradient term into matrix; implicit to avoid lagging behind the difference term
                // (d grad(p))_f = 0.5 * (d_v grad(p)_v + d_n grad(p)_n); grad(p)_w = sum_m wgt_m (p_m - p_w)
                for vid_w in [vid_v, vid_n] {
                    let rc_w = rc_vec[vid_w];
                    for (loc_w, &vid_m) in tri.ngh_vid[vid_w].iter().enumerate() {
                        let (wgt_x, wgt_y) = tri.grad_wgt[vid_w][loc_w];
                        let coeff_g = 0.5 * den_f * rc_w * (wgt_x * nx + wgt_y * ny) * face;  // rho_f * (d grad(p) . n)_f * dl
                        self.add_sclscl_a(vars, a_triplet, self.pres_id, vid_v, self.pres_id, vid_m,  coeff_g)?;
                        self.add_sclscl_a(vars, a_triplet, self.pres_id, vid_v, self.pres_id, vid_w, -coeff_g)?;
                    }
                }
            }
        }

        Ok(())
    }
}

fn _get_rc(dom: &Domain, tri: &Triangle, vel_dom: &VectorDomain, a_triplet: &Vec<(usize, usize, f64)>) -> Result<Vec<f64>, VFVChemError> {
    // get start indices
    let num_vert = dom.num_vert();
    let velx_start = match vel_dom.vecdom_type {
        VecDomType::Unknown { start } => start,
        _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
    };
    let vely_start = velx_start + num_vert;

    // initialize coefficient vectors
    let mut a_vec = vec![0.0; num_vert];

    // collect coefficients from triplets
    // only collect coefficients of velx; vely has the same coefficients
    for &(r, c, v) in a_triplet.iter() {
        if r == c && velx_start <= r && r < vely_start {  // vecx on diagonal
            a_vec[r - velx_start] += v;
        }
    }

    // get rhie-chow coefficients
    let mut rc_vec = Vec::with_capacity(num_vert);
    for vid in 0..num_vert {
        let area = tri.vert_area[vid];
        let rc = area / a_vec[vid];
        rc_vec.push(rc);
    }

    // return
    Ok(rc_vec)
}
