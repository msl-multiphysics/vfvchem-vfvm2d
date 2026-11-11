use crate::base::prelude::*;
use std::fs::create_dir_all;

#[derive(Default)]
pub struct Variables {
    // geometry
    pub pts: Points,
    pub nds: Nodes,
    pub tri: Vec<Triangle>,

    // grouping
    pub dom: Vec<Domain>,
    pub bnd: Vec<Boundary>,
    pub itf: Vec<Interface>,

    // scalar
    pub scldom: Vec<ScalarDomain>,
    pub sclbnd: Vec<ScalarBoundary>,
    pub sclitf: Vec<ScalarInterface>,

    // vector
    pub vecdom: Vec<VectorDomain>,
    pub vecbnd: Vec<VectorBoundary>,
    pub vecitf: Vec<VectorInterface>,

    // output file
    pub folder: String,
}

impl Variables {
    // geometry
    pub fn new_obj(obj_path: &str, nds_space: f64, pts_space: f64, folder: &str) -> Result<Variables, VFVChemError> {
        create_dir_all(folder).unwrap();
        let mut vars = Variables::default();
        vars.folder = folder.to_string();
        vars.nds = Nodes::new_obj(obj_path, nds_space)?;
        vars.pts = Points::new(&vars.nds, pts_space)?;
        Ok(vars)
    }

    pub fn new_obj_options(obj_path: &str, nds_space: f64, off_space: f64, crn_space: f64, crn_ratio: f64, pts_space: f64, bnd_ratio: f64, max_try: usize, seed: u64, folder: &str) -> Result<Variables, VFVChemError> {
        create_dir_all(folder).unwrap();
        let mut vars = Variables::default();
        vars.folder = folder.to_string();
        vars.nds = Nodes::new_obj_options(obj_path, nds_space, off_space, crn_space, crn_ratio)?;
        vars.pts = Points::new_options(&vars.nds, pts_space, bnd_ratio, max_try, seed)?;
        Ok(vars)
    }

    pub fn new_bounds(x_min: f64, y_min: f64, x_max: f64, y_max: f64, nds_space: f64, pts_space: f64, folder: &str) -> Result<Variables, VFVChemError> {
        create_dir_all(folder).unwrap();
        let mut vars = Variables::default();
        vars.folder = folder.to_string();
        vars.nds = Nodes::new_bounds(x_min, y_min, x_max, y_max, nds_space)?;
        vars.pts = Points::new(&vars.nds, pts_space)?;
        Ok(vars)
    }

    pub fn new_bounds_options(x_min: f64, y_min: f64, x_max: f64, y_max: f64, nds_space: f64, off_space: f64, crn_space: f64, crn_ratio: f64, pts_space: f64, bnd_ratio: f64, max_try: usize, seed: u64, folder: &str) -> Result<Variables, VFVChemError> {
        create_dir_all(folder).unwrap();
        let mut vars = Variables::default();
        vars.folder = folder.to_string();
        vars.nds = Nodes::new_bounds_options(x_min, y_min, x_max, y_max, nds_space, off_space, crn_space, crn_ratio)?;
        vars.pts = Points::new_options(&vars.nds, pts_space, bnd_ratio, max_try, seed)?;
        Ok(vars)
    }
    
    // groupings
    pub fn add_dom(&mut self, grp2d_id: &str) -> Result<usize, VFVChemError> {
        let dom_id = self.dom.len();
        let dom = Domain::new(dom_id, &self.pts, &self.nds, grp2d_id)?;
        let tri = Triangle::new(dom_id, &self.pts, &self.nds, &dom)?;  // tri_id = dom_id;
        self.dom.push(dom);
        self.tri.push(tri);
        Ok(dom_id)
    }

    pub fn add_bnd(&mut self, dom_id: usize, grp1d_id: &str) -> Result<usize, VFVChemError> {
        let bnd_id = self.bnd.len();
        let dom = &self.dom[dom_id];
        let bnd = Boundary::new(bnd_id, &self.nds, dom, grp1d_id)?;
        self.bnd.push(bnd);
        Ok(bnd_id)
    }

    pub fn add_itf(&mut self, dom_a_id: usize, dom_b_id: usize, grp1d_id: &str) -> Result<usize, VFVChemError> {
        let itf_id = self.itf.len();
        let dom_a = &self.dom[dom_a_id];
        let dom_b = &self.dom[dom_b_id];
        let itf = Interface::new(itf_id, &self.nds, dom_a, dom_b, grp1d_id)?;
        self.itf.push(itf);
        Ok(itf_id)
    }

    // scalar
    pub fn add_scldom_cnst(&mut self, dom_id: usize, value: f64, name: &str) -> Result<usize, VFVChemError> {
        let scldom_id = self.scldom.len();
        let dom = &self.dom[dom_id];
        let scldom = ScalarDomain::new_cnst(scldom_id, dom, value, name)?;
        self.scldom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_scldom_func(&mut self, dom_id: usize, func: SclDomFunc, unkn_id: Vec<usize>, name: &str) -> Result<usize, VFVChemError> {
        let scldom_id = self.scldom.len();
        let dom = &self.dom[dom_id];
        let scldom = ScalarDomain::new_func(scldom_id, dom, &self.scldom, func, unkn_id, name)?;
        self.scldom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_scldom_unkn(&mut self, dom_id: usize, value_init: f64, name: &str) -> Result<usize, VFVChemError> {
        let scldom_id = self.scldom.len();
        let dom = &self.dom[dom_id];
        let scldom = ScalarDomain::new_unkn(scldom_id, dom, value_init, name)?;
        self.scldom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_sclbnd_cnst(&mut self, bnd_id: usize, value: f64) -> Result<usize, VFVChemError> {
        let sclbnd_id = self.sclbnd.len();
        let bnd = &self.bnd[bnd_id];
        let sclbnd = ScalarBoundary::new_cnst(sclbnd_id, bnd, value)?;
        self.sclbnd.push(sclbnd);
        Ok(sclbnd_id)
    }

    pub fn add_sclbnd_func(&mut self, bnd_id: usize, func: SclBndFunc, unkn_id: Vec<usize>) -> Result<usize, VFVChemError> {
        let sclbnd_id = self.sclbnd.len();
        let bnd = &self.bnd[bnd_id];
        let sclbnd = ScalarBoundary::new_func(sclbnd_id, bnd, &self.scldom, func, unkn_id)?;
        self.sclbnd.push(sclbnd);
        Ok(sclbnd_id)
    }

    pub fn add_sclitf_cnst(&mut self, itf_id: usize, value: f64) -> Result<usize, VFVChemError> {
        let sclitf_id = self.sclitf.len();
        let itf = &self.itf[itf_id];
        let sclitf = ScalarInterface::new_cnst(sclitf_id, itf, value)?;
        self.sclitf.push(sclitf);
        Ok(sclitf_id)
    }

    pub fn add_sclitf_func(&mut self, itf_id: usize, func: SclItfFunc, unkn_a_id: Vec<usize>, unkn_b_id: Vec<usize>) -> Result<usize, VFVChemError> {
        let sclitf_id = self.sclitf.len();
        let itf = &self.itf[itf_id];
        let sclitf = ScalarInterface::new_func(sclitf_id, itf, &self.scldom, func, unkn_a_id, unkn_b_id)?;
        self.sclitf.push(sclitf);
        Ok(sclitf_id)
    }

    // vector
    pub fn add_vecdom_cnst(&mut self, dom_id: usize, value: (f64, f64), name: &str) -> Result<usize, VFVChemError> {
        let vecdom_id = self.vecdom.len();
        let dom = &self.dom[dom_id];
        let vecdom = VectorDomain::new_cnst(vecdom_id, dom, value, name)?;
        self.vecdom.push(vecdom);
        Ok(vecdom_id)
    }

    pub fn add_vecdom_func(&mut self, dom_id: usize, func: VecDomFunc, unkn_id: Vec<usize>, name: &str) -> Result<usize, VFVChemError> {
        let vecdom_id = self.vecdom.len();
        let dom = &self.dom[dom_id];
        let vecdom = VectorDomain::new_func(vecdom_id, dom, &self.scldom, func, unkn_id, name)?;
        self.vecdom.push(vecdom);
        Ok(vecdom_id)
    }

    pub fn add_vecdom_unkn(&mut self, dom_id: usize, value_init: (f64, f64), name: &str) -> Result<usize, VFVChemError> {
        let vecdom_id = self.vecdom.len();
        let dom = &self.dom[dom_id];
        let vecdom = VectorDomain::new_unkn(vecdom_id, dom, value_init, name)?;
        self.vecdom.push(vecdom);
        Ok(vecdom_id)
    }

    pub fn add_vecbnd_cnst(&mut self, bnd_id: usize, value: (f64, f64)) -> Result<usize, VFVChemError> {
        let vecbnd_id = self.vecbnd.len();
        let bnd = &self.bnd[bnd_id];
        let vecbnd = VectorBoundary::new_cnst(vecbnd_id, bnd, value)?;
        self.vecbnd.push(vecbnd);
        Ok(vecbnd_id)
    }

    pub fn add_vecbnd_func(&mut self, bnd_id: usize, func: VecBndFunc, unkn_id: Vec<usize>) -> Result<usize, VFVChemError> {
        let vecbnd_id = self.vecbnd.len();
        let bnd = &self.bnd[bnd_id];
        let vecbnd = VectorBoundary::new_func(vecbnd_id, bnd, &self.scldom, func, unkn_id)?;
        self.vecbnd.push(vecbnd);
        Ok(vecbnd_id)
    }

    pub fn add_vecitf_cnst(&mut self, itf_id: usize, value: (f64, f64)) -> Result<usize, VFVChemError> {
        let vecitf_id = self.vecitf.len();
        let itf = &self.itf[itf_id];
        let vecitf = VectorInterface::new_cnst(vecitf_id, itf, value)?;
        self.vecitf.push(vecitf);
        Ok(vecitf_id)
    }

    pub fn add_vecitf_func(&mut self, itf_id: usize, func: VecItfFunc, unkn_a_id: Vec<usize>, unkn_b_id: Vec<usize>) -> Result<usize, VFVChemError> {
        let vecitf_id = self.vecitf.len();
        let itf = &self.itf[itf_id];
        let vecitf = VectorInterface::new_func(vecitf_id, itf, &self.scldom, func, unkn_a_id, unkn_b_id)?;
        self.vecitf.push(vecitf);
        Ok(vecitf_id)
    }

}
