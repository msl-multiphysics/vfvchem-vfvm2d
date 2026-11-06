use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::transient::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct TransientHeat {
    // domain data
    pub itr_dom: Vec<usize>,              // domain
    pub itr_temp: HashMap<usize, usize>,  // temperature
    pub itr_den: HashMap<usize, usize>,   // density
    pub itr_hcap: HashMap<usize, usize>,  // heat capacity
    pub itr_cond: HashMap<usize, usize>,  // thermal conductivity
    pub itr_hlin: HashMap<usize, usize>,  // linear heat source
    pub itr_hsrc: HashMap<usize, usize>,  // heat source

    // boundary data
    pub temp_bnd: Vec<usize>,              // boundary with temperature BC
    pub temp_temp: HashMap<usize, usize>,  // temperature
    pub hflx_bnd: Vec<usize>,              // boundary with heat flux BC
    pub hflx_hflx: HashMap<usize, usize>,  // heat flux

    // interface data
    pub hcnt_itf: Vec<usize>,              // interface with continuity

    // operators
    pub oper_dom_itr: Vec<(OpSclDomTimeWeight, OpSclDomDiffusion, OpSclDomLinear, OpSclDomSource)>,
    pub oper_bnd_temp: Vec<OpSclBndDirichlet>,
    pub oper_bnd_hflx: Vec<OpSclBndNeumann>,
    pub oper_itf_hcnt: Vec<OpSclItfContinuity>,
}

impl TransientHeat {
    pub fn new() -> Result<TransientHeat, VFVChemError> {
        let phys = TransientHeat::default();
        Ok(phys)
    }

    pub fn add_heat_dom(&mut self, dom_id: usize, temp_id: usize, den_id: usize, hcap_id: usize, cond_id: usize, hlin_id: usize, hsrc_id: usize) -> Result<(), VFVChemError> {
        self.itr_dom.push(dom_id);
        self.itr_temp.insert(dom_id, temp_id);
        self.itr_den.insert(dom_id, den_id);
        self.itr_hcap.insert(dom_id, hcap_id);
        self.itr_cond.insert(dom_id, cond_id);
        self.itr_hlin.insert(dom_id, hlin_id);
        self.itr_hsrc.insert(dom_id, hsrc_id);
        Ok(())
    }

    pub fn add_temp_bnd(&mut self, bnd_id: usize, temp_id: usize) -> Result<(), VFVChemError> {
        self.temp_bnd.push(bnd_id);
        self.temp_temp.insert(bnd_id, temp_id);
        Ok(())
    }

    pub fn add_hflx_bnd(&mut self, bnd_id: usize, hflx_id: usize) -> Result<(), VFVChemError> {
        self.hflx_bnd.push(bnd_id);
        self.hflx_hflx.insert(bnd_id, hflx_id);
        Ok(())
    }

    pub fn add_hcnt_itf(&mut self, itf_id: usize) -> Result<(), VFVChemError> {
        self.hcnt_itf.push(itf_id);
        Ok(())
    }
}

impl TransientBase for TransientHeat {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknowns
        let mut xid = 0;  // index in solution vector
        for (&dom_id, &temp_id) in self.itr_temp.iter() {
            vars.scldom[temp_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_vert();
        }

        // last xid is the matrix size
        Ok(xid)
    }

    fn initial_dirichlet(&mut self, vars: &mut Variables, mat_size: usize) -> Result<(), VFVChemError> {
        // iterate over dirichlet boundaries
        for oper_dir in self.oper_bnd_temp.iter_mut() {
            oper_dir.apply_initial(vars, mat_size)?;
        }

        // iterate over continuity interfaces
        for oper_cont in self.oper_itf_hcnt.iter_mut() {
            oper_cont.apply_initial(vars, mat_size)?;
        }

        Ok(())
    }

    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError> {
        // clear operators from a previous solve
        self.oper_dom_itr.clear();
        self.oper_bnd_temp.clear();
        self.oper_bnd_hflx.clear();
        self.oper_itf_hcnt.clear();

        // domain
        for &dom_id in self.itr_dom.iter() {
            let dom_temp_id = self.itr_temp[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let dom_hcap_id = self.itr_hcap[&dom_id];
            let dom_cond_id = self.itr_cond[&dom_id];
            let dom_hlin_id = self.itr_hlin[&dom_id];
            let dom_hsrc_id = self.itr_hsrc[&dom_id];
            let oper_twgt = OpSclDomTimeWeight::new(vars, dom_id, dom_temp_id, dom_den_id, dom_hcap_id)?;
            let oper_diff = OpSclDomDiffusion::new(vars, dom_id, dom_temp_id, dom_cond_id, dom_temp_id)?;
            let oper_lin = OpSclDomLinear::new(vars, dom_id, dom_temp_id, dom_hlin_id)?;
            let oper_src = OpSclDomSource::new(vars, dom_id, dom_temp_id, dom_hsrc_id)?;
            self.oper_dom_itr.push((oper_twgt, oper_diff, oper_lin, oper_src));
        }

        // boundary
        for &bnd_id in self.hflx_bnd.iter() {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_hflx_id = self.hflx_hflx[&bnd_id];
            let oper_neu = OpSclBndNeumann::new(vars, bnd_id, dom_temp_id, bnd_hflx_id)?;
            self.oper_bnd_hflx.push(oper_neu);
        }

        // interface
        for &itf_id in self.hcnt_itf.iter() {
            let dom_a_id = vars.itf[itf_id].dom_a_id;
            let dom_b_id = vars.itf[itf_id].dom_b_id;
            let dom_temp_a_id = self.itr_temp[&dom_a_id];
            let dom_temp_b_id = self.itr_temp[&dom_b_id];
            let oper_cont = OpSclItfContinuity::new(vars, itf_id, dom_temp_a_id, dom_temp_b_id)?;
            self.oper_itf_hcnt.push(oper_cont);
        }

        // dirichlet
        for &bnd_id in self.temp_bnd.iter() {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_temp_id = self.temp_temp[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(vars, bnd_id, dom_temp_id, bnd_temp_id)?;
            self.oper_bnd_temp.push(oper_dir);
        }

        Ok(())
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>, dt: f64) -> Result<(), VFVChemError> {
        // domain
        for (oper_twgt, oper_diff, oper_lin, oper_src) in self.oper_dom_itr.iter() {
            oper_twgt.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_lin.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }

        // boundary
        for oper_neu in self.oper_bnd_hflx.iter() {
            oper_neu.apply(vars, a_triplet, b_vec)?;
        }

        // continuity
        for oper_cont in self.oper_itf_hcnt.iter() {
            oper_cont.apply(vars, a_triplet, b_vec)?;
        }

        // dirichlet
        for oper_dir in self.oper_bnd_temp.iter() {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }

        Ok(())
    }
}
