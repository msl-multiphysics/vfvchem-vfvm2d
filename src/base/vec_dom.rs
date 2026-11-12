use crate::base::prelude::*;

// input: time, (x, y), slice of unknown values at a point
// output: function value at the point
pub type VecDomFunc = Box<dyn Fn(f64, (f64, f64), &[f64]) -> (f64, f64) + Send + Sync>;

#[derive(Default)]
pub enum VecDomType {
    #[default] Constant,
    Function {
        func: VecDomFunc,  
        unkn_id: Vec<usize>,
    },
    Unknown {
        start: usize,
    },
}

#[derive(Default)]
pub struct VectorDomain {
    // ids
    pub name: String,
    pub vecdom_id: usize,
    pub dom_id: usize,  // domain this is attached to
    
    // value data
    pub vecdom_type: VecDomType,
    pub value: Vec<(f64, f64)>,  // vid -> (u, v)
    pub value_prv1: Vec<(f64, f64)>,  // 1st previous timestep
    pub value_prv2: Vec<(f64, f64)>,  // 2nd previous timestep

    // gradient data
    pub grad: Vec<(f64, f64, f64, f64)>,  // vid -> (du/dx, du/dy, dv/dx, dv/dy)
    pub grad_prv1: Vec<(f64, f64, f64, f64)>,
    pub grad_prv2: Vec<(f64, f64, f64, f64)>,
}

impl VectorDomain {
    pub fn new_cnst(vecdom_id: usize, dom: &Domain, value: (f64, f64), name: &str) -> Result<VectorDomain, VFVChemError> {
        // default class
        let mut vecdom = VectorDomain::default();
        vecdom.name = name.to_string();
        vecdom.vecdom_id = vecdom_id;
        vecdom.dom_id = dom.dom_id;

        // point and node data
        vecdom.vecdom_type = VecDomType::Constant;
        vecdom.value = vec![value; dom.num_vert()];
        vecdom.value_prv1 = vecdom.value.clone();
        vecdom.value_prv2 = vecdom.value.clone();

        // gradient data
        vecdom.grad = vec![(0.0, 0.0, 0.0, 0.0); dom.num_vert()];
        vecdom.grad_prv1 = vecdom.grad.clone();
        vecdom.grad_prv2 = vecdom.grad.clone();

        Ok(vecdom)
    }

    pub fn new_func(vecdom_id: usize, dom: &Domain, scldom_all: &Vec<ScalarDomain>, func: VecDomFunc, unkn_id: Vec<usize>, name: &str) -> Result<VectorDomain, VFVChemError> {
        // check that unkn_ids refer to unknown scalars
        for &i in unkn_id.iter() {
            match scldom_all[i].scldom_type {
                SclDomType::Unknown { start: _ } => { continue; }
                _ => { return Err(VFVChemError::VectorDomainFunctionDependsOnNonUnknown { scldom_id: i }) }
            }
        }
        
        // default class
        let mut vecdom = VectorDomain::default();
        vecdom.name = name.to_string();
        vecdom.vecdom_id = vecdom_id;
        vecdom.dom_id = dom.dom_id;

        // point data
        vecdom.vecdom_type = VecDomType::Function { func, unkn_id };
        vecdom.value = vec![(f64::NAN, f64::NAN); dom.num_vert()];  // placeholder
        vecdom.value_prv1 = vecdom.value.clone();
        vecdom.value_prv2 = vecdom.value.clone();

        // gradient data
        vecdom.grad = vec![(f64::NAN, f64::NAN, f64::NAN, f64::NAN); dom.num_vert()];  // placeholder
        vecdom.grad_prv1 = vecdom.grad.clone();
        vecdom.grad_prv2 = vecdom.grad.clone();

        Ok(vecdom)
    }

    pub fn new_unkn(vecdom_id: usize, dom: &Domain, value_init: (f64, f64), name: &str) -> Result<VectorDomain, VFVChemError> {
        // default class
        let mut vecdom = VectorDomain::default();
        vecdom.name = name.to_string();
        vecdom.vecdom_id = vecdom_id;
        vecdom.dom_id = dom.dom_id;

        // point data
        vecdom.vecdom_type = VecDomType::Unknown { start: usize::MAX };  // placeholder
        vecdom.value = vec![value_init; dom.num_vert()];
        vecdom.value_prv1 = vecdom.value.clone();
        vecdom.value_prv2 = vecdom.value.clone();

        // gradient data
        vecdom.grad = vec![(0.0, 0.0, 0.0, 0.0); dom.num_vert()];
        vecdom.grad_prv1 = vecdom.grad.clone();
        vecdom.grad_prv2 = vecdom.grad.clone();

        Ok(vecdom)
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
            let mut grad_ux = 0.0;
            let mut grad_uy = 0.0;
            let mut grad_vx = 0.0;
            let mut grad_vy = 0.0;
            for (loc, &vid_n) in list_n.iter().enumerate() {
                let (w_x, w_y) = tri.grad_wgt[vid_c][loc];
                let delta_u = self.value[vid_n].0 - self.value[vid_c].0;
                let delta_v = self.value[vid_n].1 - self.value[vid_c].1;
                grad_ux += w_x * delta_u;
                grad_uy += w_y * delta_u;
                grad_vx += w_x * delta_v;
                grad_vy += w_y * delta_v;
            }
            self.grad.push((grad_ux, grad_uy, grad_vx, grad_vy));
        }
    }

    pub fn update_func(vars: &mut Variables, vecdom_id: usize, t: f64) {
        // does not use &mut self to prevent borrowing problems
        // assume that scldoms in unkn_id have the same dom as vecdom_id

        // outline:
        // collect values from unknown-type scldom (scldom_u)
        // update the selected function-type vecdom (vecdom_f)

        // compute updated values
        let vecdom_f = &vars.vecdom[vecdom_id];
        match &vecdom_f.vecdom_type {
            VecDomType::Function { func, unkn_id } => {
                // get number of points
                let num_vert = vecdom_f.value.len();
                let dom = &vars.dom[vecdom_f.dom_id];

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

                // modify vecdom values
                let vecdom_f = &mut vars.vecdom[vecdom_id];
                vecdom_f.value = value_f;

            }
            _ => { return; }  // skip if not function-type
        };

    }

    pub fn interp_face(&self, skew: (f64, f64), vid_v: usize, vid_n: usize) -> (f64, f64) {
        // get value and gradient
        let (valu_v, valv_v) = self.value[vid_v];
        let (valu_n, valv_n) = self.value[vid_n];
        let (gradux_v, graduy_v, gradvx_v, gradvy_v) = self.grad[vid_v];
        let (gradux_n, graduy_n, gradvx_n, gradvy_n) = self.grad[vid_n];
        let (sx, sy) = skew;

        // compute value
        let valu_m = 0.5 * (valu_v + valu_n);
        let valv_m = 0.5 * (valv_v + valv_n);
        let valu_s = 0.5 * ((gradux_v + gradux_n) * sx + (graduy_v + graduy_n) * sy);
        let valv_s = 0.5 * ((gradvx_v + gradvx_n) * sx + (gradvy_v + gradvy_n) * sy);
        let valu_f = valu_m + valu_s;
        let valv_f = valv_m + valv_s;
        
        (valu_f, valv_f)
    }

    pub fn correct_face(&self, skew: (f64, f64), vid_v: usize, vid_n: usize) -> (f64, f64) {
        // get value and gradient
        let (gradux_v, graduy_v, gradvx_v, gradvy_v) = self.grad[vid_v];
        let (gradux_n, graduy_n, gradvx_n, gradvy_n) = self.grad[vid_n];
        let (sx, sy) = skew;

        // compute value
        let valu_s = 0.5 * ((gradux_v + gradux_n) * sx + (graduy_v + graduy_n) * sy);
        let valv_s = 0.5 * ((gradvx_v + gradvx_n) * sx + (gradvy_v + gradvy_n) * sy);

        (valu_s, valv_s)
    }

    pub fn interp_vert(&self, offc: (f64, f64), vid_v: usize) -> (f64, f64){
        // get value and gradient
        let (valu_v, valv_v) = self.value[vid_v];
        let (gradux_v, graduy_v, gradvx_v, gradvy_v) = self.grad[vid_v];
        let (ox, oy) = offc;

        // compute value
        let valu_o = ox * gradux_v + oy * graduy_v;
        let valv_o = ox * gradvx_v + oy * gradvy_v;
        let valu_c = valu_v + valu_o;
        let valv_c = valv_v + valv_o;

        (valu_c, valv_c)
    }

    pub fn interp_vert_prv1(&self, offc: (f64, f64), vid_v: usize) -> (f64, f64){
        // get value and gradient
        let (valu_v, valv_v) = self.value_prv1[vid_v];
        let (gradux_v, graduy_v, gradvx_v, gradvy_v) = self.grad_prv1[vid_v];
        let (ox, oy) = offc;

        // compute value
        let valu_o = ox * gradux_v + oy * graduy_v;
        let valv_o = ox * gradvx_v + oy * gradvy_v;
        let valu_c = valu_v + valu_o;
        let valv_c = valv_v + valv_o;

        (valu_c, valv_c)
    }

    pub fn interp_vert_prv2(&self, offc: (f64, f64), vid_v: usize) -> (f64, f64){
        // get value and gradient
        let (valu_v, valv_v) = self.value_prv2[vid_v];
        let (gradux_v, graduy_v, gradvx_v, gradvy_v) = self.grad_prv2[vid_v];
        let (ox, oy) = offc;

        // compute value
        let valu_o = ox * gradux_v + oy * graduy_v;
        let valv_o = ox * gradvx_v + oy * gradvy_v;
        let valu_c = valu_v + valu_o;
        let valv_c = valv_v + valv_o;

        (valu_c, valv_c)
    }

    pub fn correct_vert(&self, offc: (f64, f64), vid_v: usize) -> (f64, f64){
        // get value and gradient
        let (gradux_v, graduy_v, gradvx_v, gradvy_v) = self.grad[vid_v];
        let (ox, oy) = offc;

        // compute value
        let valu_o = ox * gradux_v + oy * graduy_v;
        let valv_o = ox * gradvx_v + oy * gradvy_v;

        (valu_o, valv_o)
    }

}
