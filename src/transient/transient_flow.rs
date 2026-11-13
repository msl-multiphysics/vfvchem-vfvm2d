use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::transient::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct TransientFlow {
    // domain data
    pub flow_dom: Vec<usize>,              // domain
    pub flow_vel: HashMap<usize, usize>,   // velocity
    pub flow_pres: HashMap<usize, usize>,  // pressure
    pub flow_den: HashMap<usize, usize>,   // density
    pub flow_visc: HashMap<usize, usize>,  // viscosity
    pub flow_fce: HashMap<usize, usize>,   // body force

    // boundary data
    pub vel_bnd: Vec<usize>,              // boundary with velocity BC
    pub vel_vel: HashMap<usize, usize>,   // velocity
    pub pres_bnd: Vec<usize>,             // boundary with pressure BC
    pub pres_pres: HashMap<usize, usize>, // pressure

    // interface data
    pub fcnt_itf: Vec<usize>,             // interface with continuity

    // reference pressure
    pub pref_use: bool,  // true -> set reference pressure
    pub pref_dom: usize,  // domain with reference pressure
    pub pref_nid: usize,  // node with reference pressure
    pub pref_pres: f64,   // reference pressure

    // operators
    pub oper_dom_flow: Vec<(FlowDomTime, FlowDomAdvection, FlowDomGradient, FlowDomDiffusion, FlowDomTranspose, FlowDomDilatatory, VecDomSource, FlowDomTimeContinuity, FlowDomDivergence)>,
    pub oper_dom_rc: Vec<FlowDomRhieChow>,
    pub oper_bnd_vel: Vec<(VecBndDirichlet, FlowBndDivergence)>,
    pub oper_bnd_pres: Vec<SclBndDirichlet>,
    pub oper_itf_fcnt: Vec<(VecItfContinuity, SclItfContinuity)>,
    pub oper_pts_ref: FlowNdsReference,
}

impl TransientFlow {
    pub fn new() -> Result<TransientFlow, VFVChemError> {
        let phys = TransientFlow::default();
        Ok(phys)
    }

    pub fn add_flow_dom(&mut self, dom_id: usize, vel_id: usize, pres_id: usize, den_id: usize, visc_id: usize, fce_id: usize) -> Result<(), VFVChemError> {
        self.flow_dom.push(dom_id);
        self.flow_vel.insert(dom_id, vel_id);
        self.flow_pres.insert(dom_id, pres_id);
        self.flow_den.insert(dom_id, den_id);
        self.flow_visc.insert(dom_id, visc_id);
        self.flow_fce.insert(dom_id, fce_id);
        Ok(())
    }

    pub fn add_vel_bnd(&mut self, bnd_id: usize, vel_id: usize) -> Result<(), VFVChemError> {
        self.vel_bnd.push(bnd_id);
        self.vel_vel.insert(bnd_id, vel_id);
        Ok(())
    }

    pub fn add_pres_bnd(&mut self, bnd_id: usize, pres_id: usize) -> Result<(), VFVChemError> {
        self.pres_bnd.push(bnd_id);
        self.pres_pres.insert(bnd_id, pres_id);
        Ok(())
    }

    pub fn add_fcnt_itf(&mut self, itf_id: usize) -> Result<(), VFVChemError> {
        self.fcnt_itf.push(itf_id);
        Ok(())
    }

    pub fn set_pres_ref(&mut self, dom_id: usize, nid: usize, pres: f64) -> Result<(), VFVChemError> {
        self.pref_use = true;
        self.pref_dom = dom_id;
        self.pref_nid = nid;
        self.pref_pres = pres;
        Ok(())
    }

}

impl TransientBase for TransientFlow {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknowns
        let mut xid = 0;  // index in solution vector
        for (&dom_id, &vel_id) in self.flow_vel.iter() {
            vars.vecdom[vel_id].vecdom_type = VecDomType::Unknown { start: xid };
            xid += 2 * vars.dom[dom_id].num_vert();
        }
        for (&dom_id, &pres_id) in self.flow_pres.iter() {
            vars.scldom[pres_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_vert();
        }

        // last xid is the matrix size
        Ok(xid)
    }

    fn initial_dirichlet(&mut self, vars: &mut Variables, mat_size: usize) -> Result<(), VFVChemError> {
        // iterate over dirichlet boundaries
        for (oper_dir, _) in self.oper_bnd_vel.iter_mut() {
            oper_dir.apply_initial(vars, mat_size)?;
        }
        for oper_dir in self.oper_bnd_pres.iter_mut() {
            oper_dir.apply_initial(vars, mat_size)?;
        }

        // iterate over continuity interfaces
        for (oper_vcnt, oper_pcnt) in self.oper_itf_fcnt.iter_mut() {
            oper_vcnt.apply_initial(vars, mat_size)?;
            oper_pcnt.apply_initial(vars, mat_size)?;
        }

        Ok(())
    }

    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError> {
        // clear operators from a previous solve
        self.oper_dom_flow.clear();
        self.oper_dom_rc.clear();
        self.oper_bnd_vel.clear();
        self.oper_bnd_pres.clear();
        self.oper_itf_fcnt.clear();
        
        // domain
        for &dom_id in self.flow_dom.iter() {
            let dom_vel_id = self.flow_vel[&dom_id];
            let dom_pres_id = self.flow_pres[&dom_id];
            let dom_den_id = self.flow_den[&dom_id];
            let dom_visc_id = self.flow_visc[&dom_id];
            let dom_fce_id = self.flow_fce[&dom_id];
            let oper_time = FlowDomTime::new(vars, dom_id, dom_vel_id, dom_den_id)?;
            let oper_adv = FlowDomAdvection::new(vars, dom_id, dom_vel_id, dom_den_id, dom_visc_id)?;
            let oper_grad = FlowDomGradient::new(vars, dom_id, dom_vel_id, dom_pres_id)?;
            let oper_diff = FlowDomDiffusion::new(vars, dom_id, dom_vel_id, dom_visc_id)?;
            let oper_dtrn = FlowDomTranspose::new(vars, dom_id, dom_vel_id, dom_visc_id)?;
            let oper_ddil = FlowDomDilatatory::new(vars, dom_id, dom_vel_id, dom_visc_id)?;
            let oper_src = VecDomSource::new(vars, dom_id, dom_vel_id, dom_fce_id)?;
            let oper_tden = FlowDomTimeContinuity::new(vars, dom_id, dom_pres_id, dom_den_id)?;
            let oper_div = FlowDomDivergence::new(vars, dom_id, dom_pres_id, dom_vel_id, dom_den_id)?;
            let oper_rc = FlowDomRhieChow::new(vars, dom_id, dom_pres_id, dom_vel_id, dom_den_id)?;
            self.oper_dom_flow.push((oper_time, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_tden, oper_div));
            self.oper_dom_rc.push(oper_rc);
        }

        // interface
        for &itf_id in self.fcnt_itf.iter() {
            let dom_a_id = vars.itf[itf_id].dom_a_id;
            let dom_b_id = vars.itf[itf_id].dom_b_id;
            let dom_vel_a_id = self.flow_vel[&dom_a_id];
            let dom_vel_b_id = self.flow_vel[&dom_b_id];
            let dom_pres_a_id = self.flow_pres[&dom_a_id];
            let dom_pres_b_id = self.flow_pres[&dom_b_id];
            let oper_vcnt = VecItfContinuity::new(vars, itf_id, dom_vel_a_id, dom_vel_b_id)?;
            let oper_pcnt = SclItfContinuity::new(vars, itf_id, dom_pres_a_id, dom_pres_b_id)?;
            self.oper_itf_fcnt.push((oper_vcnt, oper_pcnt));
        }

        // dirichlet
        for &bnd_id in self.vel_bnd.iter() {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.flow_vel[&dom_id];
            let dom_pres_id = self.flow_pres[&dom_id];
            let dom_den_id = self.flow_den[&dom_id];
            let bnd_vel_id = self.vel_vel[&bnd_id];
            let oper_dir = VecBndDirichlet::new(vars, bnd_id, dom_vel_id, bnd_vel_id)?;
            let oper_div = FlowBndDivergence::new(vars, bnd_id, dom_vel_id, dom_den_id, dom_pres_id)?;
            self.oper_bnd_vel.push((oper_dir, oper_div));
        }
        for &bnd_id in self.pres_bnd.iter() {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_pres_id = self.flow_pres[&dom_id];
            let bnd_pres_id = self.pres_pres[&bnd_id];
            let oper_dir = SclBndDirichlet::new(vars, bnd_id, dom_pres_id, bnd_pres_id)?;
            self.oper_bnd_pres.push(oper_dir);
        }

        // reference pressure
        if self.pref_use {
            let dom_id = self.pref_dom;
            let unk_id = self.flow_pres[&dom_id];
            self.oper_pts_ref = FlowNdsReference::new(vars, dom_id, unk_id, self.pref_nid, self.pref_pres)?;
        }

        Ok(())
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, twgt: (f64, f64, f64)) -> Result<(), VFVChemError> {
        // domain
        for (oper_time, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_tden, oper_div) in self.oper_dom_flow.iter() {
            oper_time.apply_time(vars, a_triplet, b_vec, twgt)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_tden.apply_time(vars, a_triplet, b_vec, twgt)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
        }

        // boundary
        for (_, oper_div) in self.oper_bnd_vel.iter() {
            oper_div.apply(vars, a_triplet, b_vec)?;
        }

        // rhie-chow stabilization
        for oper_rc in self.oper_dom_rc.iter() {
            oper_rc.apply_time(vars, a_triplet, b_vec, twgt)?;
        }

        // continuity
        for (oper_vcnt, oper_pcnt) in self.oper_itf_fcnt.iter() {
            oper_vcnt.apply(vars, a_triplet, b_vec)?;
            oper_pcnt.apply(vars, a_triplet, b_vec)?;
        }

        // dirichlet
        for (oper_dir, _) in self.oper_bnd_vel.iter() {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }
        for oper_dir in self.oper_bnd_pres.iter() {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }

        // reference pressure
        if self.pref_use {
            self.oper_pts_ref.apply(vars, a_triplet, b_vec)?;
        }

        Ok(())
    }
}
