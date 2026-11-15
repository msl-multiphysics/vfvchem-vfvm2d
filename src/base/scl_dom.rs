use crate::base::prelude::*;

// input: time, (x, y), slice of unknown values at a point
// output: function value at the point
pub type SclDomFunc = Box<dyn Fn(f64, (f64, f64), &[f64]) -> f64 + Send + Sync>;

#[derive(Default)]
pub enum SclDomType {
    #[default] Constant,
    Function {
        func: SclDomFunc,  
        unkn_id: Vec<usize>,
    },
    Unknown {
        start: usize,
    },
}

#[derive(Default)]
pub struct ScalarDomain {
    // ids
    pub name: String,
    pub scldom_id: usize,
    pub dom_id: usize,  // domain this is attached to
    
    // value data
    pub scldom_type: SclDomType,
    pub value: Vec<f64>,  // vid -> value (v)
    pub value_prv1: Vec<f64>,  // 1st previous timestep
    pub value_prv2: Vec<f64>,  // 2nd previous timestep

    // gradient data
    pub grad: Vec<(f64, f64)>,  // vid -> (dv/dx, dv/dy)
    pub grad_prv1: Vec<(f64, f64)>,
    pub grad_prv2: Vec<(f64, f64)>,
}

impl ScalarDomain {
    pub fn new_cnst(scldom_id: usize, dom: &Domain, value: f64, name: &str) -> Result<ScalarDomain, VFVChemError> {
        // default class
        let mut scldom = ScalarDomain::default();
        scldom.name = name.to_string();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // point and node data
        scldom.scldom_type = SclDomType::Constant;
        scldom.value = vec![value; dom.num_vert()];
        scldom.value_prv1 = scldom.value.clone();
        scldom.value_prv2 = scldom.value.clone();

        // gradient data
        scldom.grad = vec![(0.0, 0.0); dom.num_vert()];
        scldom.grad_prv1 = scldom.grad.clone();
        scldom.grad_prv2 = scldom.grad.clone();

        Ok(scldom)
    }

    pub fn new_func(scldom_id: usize, dom: &Domain, scldom_all: &Vec<ScalarDomain>, func: SclDomFunc, unkn_id: Vec<usize>, name: &str) -> Result<ScalarDomain, VFVChemError> {
        // check that unkn_ids refer to unknown scalars
        for &i in unkn_id.iter() {
            match scldom_all[i].scldom_type {
                SclDomType::Unknown { start: _ } => { continue; }
                _ => { return Err(VFVChemError::ScalarDomainFunctionDependsOnNonUnknown { scldom_id: i }) }
            }
        }
        
        // default class
        let mut scldom = ScalarDomain::default();
        scldom.name = name.to_string();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // point data
        scldom.scldom_type = SclDomType::Function { func, unkn_id };
        scldom.value = vec![f64::NAN; dom.num_vert()];  // placeholder
        scldom.value_prv1 = scldom.value.clone();
        scldom.value_prv2 = scldom.value.clone();

        // gradient data
        scldom.grad = vec![(f64::NAN, f64::NAN); dom.num_vert()];  // placeholder
        scldom.grad_prv1 = scldom.grad.clone();
        scldom.grad_prv2 = scldom.grad.clone();

        Ok(scldom)
    }

    pub fn new_unkn(scldom_id: usize, dom: &Domain, value_init: f64, name: &str) -> Result<ScalarDomain, VFVChemError> {
        // default class
        let mut scldom = ScalarDomain::default();
        scldom.name = name.to_string();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // point data
        scldom.scldom_type = SclDomType::Unknown { start: usize::MAX };  // placeholder
        scldom.value = vec![value_init; dom.num_vert()];
        scldom.value_prv1 = scldom.value.clone();
        scldom.value_prv2 = scldom.value.clone();

        // gradient data
        scldom.grad = vec![(0.0, 0.0); dom.num_vert()];
        scldom.grad_prv1 = scldom.grad.clone();
        scldom.grad_prv2 = scldom.grad.clone();

        Ok(scldom)
    }

    pub fn update_prev(&mut self) {
        self.value_prv2 = self.value_prv1.clone();
        self.grad_prv2 = self.grad_prv1.clone();
        self.value_prv1 = self.value.clone();
        self.grad_prv1 = self.grad.clone();
    }

    pub fn update_grad(&mut self, tri: &Triangle) {
        // initialize
        self.grad.clear();

        // iterate through neighbors
        for (vid_c, list_n) in tri.ngh_vid.iter().enumerate() {
            let mut grad_x = 0.0;
            let mut grad_y = 0.0;
            for (loc, &vid_n) in list_n.iter().enumerate() {
                let (w_x, w_y) = tri.grad_wgt[vid_c][loc];
                let delta_value = self.value[vid_n] - self.value[vid_c];
                grad_x += w_x * delta_value;
                grad_y += w_y * delta_value;
            }
            self.grad.push((grad_x, grad_y));
        }
    }

    pub fn update_func(vars: &mut Variables, scldom_id: usize, t: f64) {
        // does not use &mut self to prevent borrowing problems
        // assume that scldoms in unkn_id have the same dom as scldom_id

        // outline:
        // collect values from unknown-type scldom (scldom_u)
        // update the selected function-type scldom (scldom_f)

        // compute updated values
        let scldom_f = &vars.scldom[scldom_id];
        match &scldom_f.scldom_type {
            SclDomType::Function { func, unkn_id } => {
                // get number of points
                let num_vert = scldom_f.value.len();
                let dom = &vars.dom[scldom_f.dom_id];

                // collect positions
                let mut pos = Vec::with_capacity(num_vert);  // vid -> (x, y)
                for &gpt in &dom.pts_gpt {
                    pos.push(vars.pts.pts_pos[gpt]);
                }
                for &gnd in &dom.nds_gnd {
                    pos.push(vars.nds.nds_pos[gnd]);
                }
                
                // collect values from unknown scalars
                let mut value_u = Vec::with_capacity(num_vert);  // vid -> values of unknowns
                for vid in 0..num_vert {
                    let mut value_sub = Vec::with_capacity(unkn_id.len());  // values of unknowns
                    for &i in unkn_id {
                        let scldom_u = &vars.scldom[i];
                        value_sub.push(scldom_u.value[vid]);
                    }
                    value_u.push(value_sub);
                }

                // compute function value
                let mut value_f = Vec::with_capacity(num_vert);
                for vid in 0..num_vert {
                    let value_sub = func(t, pos[vid], value_u[vid].as_slice());
                    value_f.push(value_sub)
                }

                // modify scldom values
                let scldom_f = &mut vars.scldom[scldom_id];
                scldom_f.value = value_f;

            }
            _ => { return; }  // skip if not function-type
        };

    }

    pub fn interp_face(&self, skew: (f64, f64), vid_v: usize, vid_n: usize) -> f64 {
        // get value and gradient
        let value_v = self.value[vid_v];
        let value_n = self.value[vid_n];
        let (gradx_v, grady_v) = self.grad[vid_v];
        let (gradx_n, grady_n) = self.grad[vid_n];
        let (sx, sy) = skew;

        // compute value
        let value_m = 0.5 * (value_v + value_n);
        let value_s = 0.5 * ((gradx_v + gradx_n) * sx + (grady_v + grady_n) * sy);
        let value_f = value_m + value_s;
        
        value_f
    }

    pub fn interp_face_harmonic(&self, skew: (f64, f64), vid_v: usize, vid_n: usize) -> f64 {
        // get value and gradient
        let value_v = self.value[vid_v];
        let value_n = self.value[vid_n];
        let (gradx_v, grady_v) = self.grad[vid_v];
        let (gradx_n, grady_n) = self.grad[vid_n];
        let (sx, sy) = skew;

        // compute value
        let value_m = if value_v * value_n > 0.0 { 2.0 * value_v * value_n / (value_v + value_n) } else { 0.0 };  // guard against zero or opposite signs
        let value_s = 0.5 * ((gradx_v + gradx_n) * sx + (grady_v + grady_n) * sy);
        let value_f = value_m + value_s;

        value_f
    }

    pub fn correct_face(&self, skew: (f64, f64), vid_v: usize, vid_n: usize) -> f64 {
        // get value and gradient
        let (gradx_v, grady_v) = self.grad[vid_v];
        let (gradx_n, grady_n) = self.grad[vid_n];
        let (sx, sy) = skew;

        // compute value
        let value_s = 0.5 * ((gradx_v + gradx_n) * sx + (grady_v + grady_n) * sy);

        value_s
    }

    pub fn interp_vert(&self, offc: (f64, f64), vid_v: usize) -> f64 {
        // get value and gradient
        let value_v = self.value[vid_v];
        let (gradx_v, grady_v) = self.grad[vid_v];
        let (ox, oy) = offc;

        // compute value
        let value_o = ox * gradx_v + oy * grady_v;
        let value_c = value_v + value_o;

        value_c
    }

    pub fn interp_vert_prv1(&self, offc: (f64, f64), vid_v: usize) -> f64 {
        // get value and gradient
        let value_v = self.value_prv1[vid_v];
        let (gradx_v, grady_v) = self.grad_prv1[vid_v];
        let (ox, oy) = offc;

        // compute value
        let value_o = ox * gradx_v + oy * grady_v;
        let value_c = value_v + value_o;

        value_c
    }

    pub fn interp_vert_prv2(&self, offc: (f64, f64), vid_v: usize) -> f64 {
        // get value and gradient
        let value_v = self.value_prv2[vid_v];
        let (gradx_v, grady_v) = self.grad_prv2[vid_v];
        let (ox, oy) = offc;

        // compute value
        let value_o = ox * gradx_v + oy * grady_v;
        let value_c = value_v + value_o;

        value_c
    }

    pub fn correct_vert(&self, offc: (f64, f64), vid_v: usize) -> f64 {
        // get value and gradient
        let (gradx_v, grady_v) = self.grad[vid_v];
        let (ox, oy) = offc;

        // compute value
        let value_o = ox * gradx_v + oy * grady_v;

        value_o
    }

}
