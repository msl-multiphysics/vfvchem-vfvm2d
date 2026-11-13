use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::steady::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct SteadyMole {
    // domain data
    pub num_comp: usize,                                        // number of components
    pub mole_dom: Vec<Vec<usize>>,                              // comp_id -> domain
    pub mole_conc: Vec<HashMap<usize, usize>>,                  // comp_id -> concentration
    pub mole_den: Vec<HashMap<usize, usize>>,                   // comp_id -> density
    pub mole_diff: Vec<HashMap<usize, HashMap<usize, usize>>>,  // comp_i, comp_j -> diffusion coefficient
    pub mole_nlin: Vec<HashMap<usize, usize>>,                  // comp_id -> linear molar source
    pub mole_nsrc: Vec<HashMap<usize, usize>>,                  // comp_id -> molar source

    // boundary data
    pub conc_bnd: Vec<Vec<usize>>,              // comp_id -> boundary with concentration BC
    pub conc_conc: Vec<HashMap<usize, usize>>,  // comp_id -> concentration
    pub nflx_bnd: Vec<Vec<usize>>,              // comp_id -> boundary with molar flux BC
    pub nflx_nflx: Vec<HashMap<usize, usize>>,  // comp_id -> molar flux

    // interface data
    pub ncnt_itf: Vec<Vec<usize>>,              // comp_id -> interface with continuity

    // operators
    pub oper_dom_mole: Vec<(Vec<MoleDomDiffusion>, SclDomLinear, SclDomSource)>,
    pub oper_bnd_conc: Vec<SclBndDirichlet>,
    pub oper_bnd_nflx: Vec<SclBndNeumann>,
    pub oper_itf_ncnt: Vec<SclItfContinuity>,
}

impl SteadyMole {
    pub fn new(num_comp: usize) -> Result<SteadyMole, VFVChemError> {
        // input checking
        if num_comp == 0 { return Err(VFVChemError::MoleNonPositiveComponents {}) }

        // default struct
        let mut phys = SteadyMole::default();
        phys.num_comp = num_comp;

        // initialize fields
        phys.mole_dom = vec![Vec::new(); num_comp];
        phys.mole_conc = vec![HashMap::new(); num_comp];
        phys.mole_den = vec![HashMap::new(); num_comp];
        phys.mole_diff = vec![HashMap::new(); num_comp];
        phys.mole_nlin = vec![HashMap::new(); num_comp];
        phys.mole_nsrc = vec![HashMap::new(); num_comp];
        phys.conc_bnd = vec![Vec::new(); num_comp];
        phys.conc_conc = vec![HashMap::new(); num_comp];
        phys.nflx_bnd = vec![Vec::new(); num_comp];
        phys.nflx_nflx = vec![HashMap::new(); num_comp];
        phys.ncnt_itf = vec![Vec::new(); num_comp];

        Ok(phys)
    }

    pub fn add_mole_dom(&mut self, comp_id: usize, dom_id: usize, conc_id: usize, den_id: usize, diff_id: HashMap<usize, usize>, nlin_id: usize, nsrc_id: usize) -> Result<(), VFVChemError> {
        self.mole_dom[comp_id].push(dom_id);
        self.mole_conc[comp_id].insert(dom_id, conc_id);
        self.mole_den[comp_id].insert(dom_id, den_id);
        self.mole_diff[comp_id].insert(dom_id, diff_id);
        self.mole_nlin[comp_id].insert(dom_id, nlin_id);
        self.mole_nsrc[comp_id].insert(dom_id, nsrc_id);
        Ok(())
    }

    pub fn add_conc_bnd(&mut self, comp_id: usize, bnd_id: usize, conc_id: usize) -> Result<(), VFVChemError> {
        self.conc_bnd[comp_id].push(bnd_id);
        self.conc_conc[comp_id].insert(bnd_id, conc_id);
        Ok(())
    }

    pub fn add_nflx_bnd(&mut self, comp_id: usize, bnd_id: usize, nflx_id: usize) -> Result<(), VFVChemError> {
        self.nflx_bnd[comp_id].push(bnd_id);
        self.nflx_nflx[comp_id].insert(bnd_id, nflx_id);
        Ok(())
    }

    pub fn add_ncnt_itf(&mut self, comp_id: usize, itf_id: usize) -> Result<(), VFVChemError> {
        self.ncnt_itf[comp_id].push(itf_id);
        Ok(())
    }
}

impl SteadyBase for SteadyMole {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknowns
        let mut xid = 0;  // index in solution vector
        for comp_id in 0..self.num_comp {
        for (&dom_id, &conc_id) in self.mole_conc[comp_id].iter() {
            vars.scldom[conc_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_vert();
        }}
        
        // last xid is the matrix size
        Ok(xid)
    }

    fn initial_dirichlet(&mut self, vars: &mut Variables, mat_size: usize) -> Result<(), VFVChemError> {
        // iterate over dirichlet boundaries
        for oper_dir in self.oper_bnd_conc.iter_mut() {
            oper_dir.apply_initial(vars, mat_size)?;
        }

        // iterate over continuity interfaces
        for oper_cont in self.oper_itf_ncnt.iter_mut() {
            oper_cont.apply_initial(vars, mat_size)?;
        }

        Ok(())
    }

    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError> {
        // clear operators from a previous solve
        self.oper_dom_mole.clear();
        self.oper_bnd_conc.clear();
        self.oper_bnd_nflx.clear();
        self.oper_itf_ncnt.clear();

        // domain
        for comp_id in 0..self.num_comp {
        for &dom_id in self.mole_dom[comp_id].iter() {
            // properties for comp_i
            let dom_conc_id = self.mole_conc[comp_id][&dom_id];
            let dom_den_id = self.mole_den[comp_id][&dom_id];
            let dom_nlin_id = self.mole_nlin[comp_id][&dom_id];
            let dom_nsrc_id = self.mole_nsrc[comp_id][&dom_id];
            let oper_lin = SclDomLinear::new(vars, dom_id, dom_conc_id, dom_nlin_id)?;
            let oper_src = SclDomSource::new(vars, dom_id, dom_conc_id, dom_nsrc_id)?;

            // cross-properties for comp_i, comp_j
            let mut oper_diff_vec = Vec::new();
            for (&comp_j, &dom_diff_id) in self.mole_diff[comp_id][&dom_id].iter() {
                let dom_drv_id = self.mole_conc[comp_j][&dom_id];
                let oper_diff = MoleDomDiffusion::new(vars, dom_id, dom_conc_id, dom_den_id, dom_diff_id, dom_drv_id)?;
                oper_diff_vec.push(oper_diff);
            }

            // store operators
            self.oper_dom_mole.push((oper_diff_vec, oper_lin, oper_src));
        }}

        // boundary
        for comp_id in 0..self.num_comp {
        for &bnd_id in self.nflx_bnd[comp_id].iter() {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_conc_id = self.mole_conc[comp_id][&dom_id];
            let bnd_nflx_id = self.nflx_nflx[comp_id][&bnd_id];
            let oper_neu = SclBndNeumann::new(vars, bnd_id, dom_conc_id, bnd_nflx_id)?;
            self.oper_bnd_nflx.push(oper_neu);
        }}

        // interface
        for comp_id in 0..self.num_comp {
        for &itf_id in self.ncnt_itf[comp_id].iter() {
            let dom_a_id = vars.itf[itf_id].dom_a_id;
            let dom_b_id = vars.itf[itf_id].dom_b_id;
            let dom_conc_a_id = self.mole_conc[comp_id][&dom_a_id];
            let dom_conc_b_id = self.mole_conc[comp_id][&dom_b_id];
            let oper_cont = SclItfContinuity::new(vars, itf_id, dom_conc_a_id, dom_conc_b_id)?;
            self.oper_itf_ncnt.push(oper_cont);
        }}

        // dirichlet
        for comp_id in 0..self.num_comp {
        for &bnd_id in self.conc_bnd[comp_id].iter() {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_conc_id = self.mole_conc[comp_id][&dom_id];
            let bnd_conc_id = self.conc_conc[comp_id][&bnd_id];
            let oper_dir = SclBndDirichlet::new(vars, bnd_id, dom_conc_id, bnd_conc_id)?;
            self.oper_bnd_conc.push(oper_dir);
        }}

        Ok(())
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<(usize, usize, f64)>, b_vec: &mut Vec<f64>) -> Result<(), VFVChemError> {
        // domain
        for (oper_diff_vec, oper_lin, oper_src) in self.oper_dom_mole.iter() {
            // properties for comp_i
            oper_lin.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;

            // cross-properties for comp_i, comp_j
            for oper_diff in oper_diff_vec.iter() {
                oper_diff.apply(vars, a_triplet, b_vec)?;
            }
        }

        // boundary
        for oper_neu in self.oper_bnd_nflx.iter() {
            oper_neu.apply(vars, a_triplet, b_vec)?;
        }

        // continuity
        for oper_cont in self.oper_itf_ncnt.iter() {
            oper_cont.apply(vars, a_triplet, b_vec)?;
        }

        // dirichlet
        for oper_dir in self.oper_bnd_conc.iter() {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }

        Ok(())
    }
}
