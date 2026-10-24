use crate::base::prelude::*;
use crate::operator::extend::pore_helper;
use crate::operator::prelude::*;
use crate::steady::prelude::*;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::HashMap;

#[derive(Default)]
pub struct SteadyHeatPore {
    // fluid heat internal data
    // indexed by domain
    pub itr_templ: HashMap<usize, usize>,  // fluid temperature
    pub itr_vlcpl: HashMap<usize, usize>,  // fluid volumetric heat capacity (rho * cp)_f
    pub itr_condl: HashMap<usize, usize>,  // fluid thermal conductivity
    pub itr_hprdl: HashMap<usize, usize>,  // fluid heat production coefficient; per fluid volume
    pub itr_hsrcl: HashMap<usize, usize>,  // fluid heat source; per fluid volume

    // solid heat internal data
    // indexed by domain
    pub itr_temps: HashMap<usize, usize>,  // solid temperature
    pub itr_vlcps: HashMap<usize, usize>,  // solid volumetric heat capacity (rho * c)_s
    pub itr_conds: HashMap<usize, usize>,  // solid thermal conductivity
    pub itr_hprds: HashMap<usize, usize>,  // solid heat production coefficient; per solid volume
    pub itr_hsrcs: HashMap<usize, usize>,  // solid heat source; per solid volume
    pub itr_aprds: HashMap<usize, usize>,  // surface heat production coefficient; per solid surface area
    pub itr_asrcs: HashMap<usize, usize>,  // surface heat source; per solid surface area

    // heat matrix data
    // indexed by domain
    pub itr_htc: HashMap<usize, usize>,    // interstitial heat transfer coefficient
    pub itr_area: HashMap<usize, usize>,   // specific surface area
    pub itr_displ: HashMap<usize, usize>,  // longitudinal thermal dispersivity
    pub itr_dispt: HashMap<usize, usize>,  // transverse thermal dispersivity

    // fluid heat boundary data
    pub templ_bnd: Vec<usize>,               // boundary with fluid temperature BC
    pub templ_templ: HashMap<usize, usize>,  // fluid temperature
    pub hflxl_bnd: Vec<usize>,               // boundary with fluid heat flux BC
    pub hflxl_hflxl: HashMap<usize, usize>,  // fluid heat flux; conductive plus dispersive
    pub houtl_bnd: Vec<usize>,               // boundary with fluid heat outflow BC

    // solid heat boundary data
    pub temps_bnd: Vec<usize>,               // boundary with solid temperature BC
    pub temps_temps: HashMap<usize, usize>,  // solid temperature
    pub hflxs_bnd: Vec<usize>,               // boundary with solid heat flux BC
    pub hflxs_hflxs: HashMap<usize, usize>,  // solid conductive heat flux

    // heat interface data
    pub hcntl_itf: Vec<usize>,  // interface with fluid temperature and heat flux continuity IC
    pub hcnts_itf: Vec<usize>,  // interface with solid temperature and heat flux continuity IC

    // flow internal data
    pub itr_dom: Vec<usize>,              // domain
    pub itr_vel: HashMap<usize, usize>,   // superficial velocity
    pub itr_pres: HashMap<usize, usize>,  // pressure
    pub itr_den: HashMap<usize, usize>,   // density
    pub itr_visc: HashMap<usize, usize>,  // viscosity
    pub itr_fce: HashMap<usize, usize>,   // body force
    pub itr_epsl: HashMap<usize, usize>,  // porosity
    pub itr_perm: HashMap<usize, usize>,  // permeability
    pub itr_forch: HashMap<usize, usize>, // Forchheimer coefficient
    pub itr_flw: HashMap<usize, usize>,   // flow data; created in initial_start

    // flow boundary data
    pub vel_bnd: Vec<usize>,              // boundary with velocity BC
    pub vel_vel: HashMap<usize, usize>,   // velocity
    pub pres_bnd: Vec<usize>,             // boundary with pressure BC
    pub pres_pres: HashMap<usize, usize>, // pressure
    pub trac_bnd: Vec<usize>,             // boundary with traction BC
    pub trac_pres: HashMap<usize, usize>, // pressure

    // flow interface data
    pub fcnt_itf: Vec<usize>,             // interface with continuity IC

    // reference pressure
    pub pref_dom: Vec<usize>,             // domain with reference pressure
    pub pref_nid: HashMap<usize, usize>,  // node with reference pressure
    pub pref_pres: HashMap<usize, f64>,   // reference pressure

    // fluid heat operators
    pub oper_itr_heatl: Vec<(PoreOpSclDomHeatDiffusion, PoreOpSclDomDispersion, PoreOpSclDomHeatAdvection, PoreOpSclDomProduction, PoreOpSclDomSource)>,
    pub oper_itr_hexcl: Vec<PoreOpSclDomHeatExchange>,
    pub oper_bnd_templ: Vec<OpSclBndDirichlet>,
    pub oper_bnd_hflxl: Vec<(PoreOpSclBndHeatDiffusion, PoreOpSclBndDispersion, PoreOpSclBndHeatAdvection, PoreOpSclBndProduction, PoreOpSclBndSource, OpSclBndNeumann)>,
    pub oper_bnd_houtl: Vec<(PoreOpSclBndHeatDiffusion, PoreOpSclBndDispersion, PoreOpSclBndHeatAdvection, PoreOpSclBndProduction, PoreOpSclBndSource)>,
    pub oper_bnd_hexcl: Vec<PoreOpSclBndHeatExchange>,
    pub oper_itf_hcntl: Vec<(OpSclItfContinuity, PoreOpSclItfHeatDiffusion, PoreOpSclItfProduction, PoreOpSclItfSource)>,
    pub oper_itf_hadvl: Vec<(PoreOpSclItfDispersion, PoreOpSclItfHeatAdvection)>,

    // solid heat operators
    pub oper_itr_heats: Vec<(PoreOpSclDomSolidDiffusion, PoreOpSclDomSolidProduction, PoreOpSclDomSolidSource, PoreOpSclDomProduction, PoreOpSclDomSource, PoreOpSclDomHeatExchange)>,
    pub oper_bnd_temps: Vec<OpSclBndDirichlet>,
    pub oper_bnd_hflxs: Vec<(PoreOpSclBndSolidDiffusion, PoreOpSclBndSolidProduction, PoreOpSclBndSolidSource, PoreOpSclBndProduction, PoreOpSclBndSource, PoreOpSclBndHeatExchange, OpSclBndNeumann)>,
    pub oper_itf_hcnts: Vec<OpSclItfContinuity>,
    pub oper_itf_heats: Vec<(PoreOpSclItfSolidDiffusion, PoreOpSclItfSolidProduction, PoreOpSclItfSolidSource, PoreOpSclItfSurfaceProduction, PoreOpSclItfSurfaceSource, PoreOpSclItfHeatExchange, PoreOpSclItfHeatExchange)>,

    // flow operators
    pub oper_itr: Vec<(PoreOpVecDomAdvection, OpVecDomGradient, PoreOpVecDomDiffusion, PoreOpVecDomTranspose, PoreOpVecDomDilatatory, PoreOpVecDomDrag, OpVecDomSource, PoreOpSclDomDivergence, OpSclDomPressureStabilizer)>,
    pub oper_bnd_vel: Vec<(OpVecBndDirichlet, PoreOpSclBndDivergence, OpSclBndPressureStabilizer)>,
    pub oper_bnd_pres: Vec<(OpSclBndDirichlet, PoreOpVecBndAdvection, OpVecBndGradient, PoreOpVecBndDiffusion, PoreOpVecBndTranspose, PoreOpVecBndDilatatory, PoreOpVecBndDrag, OpVecBndSource)>,
    pub oper_bnd_trac: Vec<(OpVecBndTraction, PoreOpVecBndAdvection, PoreOpVecBndDiffusion, PoreOpVecBndTranspose, PoreOpVecBndDilatatory, PoreOpVecBndDrag, OpVecBndSource, PoreOpSclBndDivergence, OpSclBndPressureStabilizer)>,
    pub oper_itf_fcnt: Vec<(OpVecItfContinuity, OpSclItfContinuity, PoreOpVecItfAdvection, OpVecItfGradient, PoreOpVecItfDiffusion, PoreOpVecItfTranspose, PoreOpVecItfDilatatory, PoreOpVecItfDrag, OpVecItfSource, PoreOpSclItfDivergence, OpSclItfPressureStabilizer)>,
    pub oper_pref: Vec<OpSclDomPressureReference>,
}

impl SteadyHeatPore {
    pub fn new() -> Result<SteadyHeatPore, VFVChemError> {
        // local thermal non-equilibrium heat transfer in a porous medium; fluid (l) and solid (s) temperatures
        // fluid: div( (rho cp)_f * vel * T_f ) = div( (epsl * k_f * I + K_d) * grad(T_f) ) + h * a_v * (T_s - T_f) + epsl * (s_f * T_f + q_f)
        // solid: 0 = div( (1 - epsl) * k_s * grad(T_s) ) - h * a_v * (T_s - T_f) + (1 - epsl) * (s_s * T_s + q_s) + a_v * (s_a * T_s + q_a)
        // vel is the superficial velocity; K_d = (rho cp)_f * (alpha_T * |v| * I + (alpha_L - alpha_T) * v v^T / |v|) with v = vel / epsl
        // q_f and q_s are per unit phase volume; q_a is per unit solid surface area; the exchange is implicit
        // every fluid heat domain must be a pore domain with matrix data; a solid heat domain is optional and needs fluid heat data
        // at vertices where epsl = 1 on the vertex and its neighbors the solid equation is replaced with T_s = T_f
        // a free fluid is a pore domain with epsl = 1, K = inf, F = 0, and alpha_L = alpha_T = 0, without a solid heat domain
        // interfaces join the fluid temperatures; the solid temperatures are joined if both sides have a solid, otherwise the solid surface is adiabatic
        Ok(SteadyHeatPore::default())
    }

    pub fn add_heatl_dom(&mut self, dom_id: usize, templ_id: usize, vlcpl_id: usize, condl_id: usize, hprdl_id: usize, hsrcl_id: usize) -> Result<(), VFVChemError> {
        self.itr_templ.insert(dom_id, templ_id);
        self.itr_vlcpl.insert(dom_id, vlcpl_id);
        self.itr_condl.insert(dom_id, condl_id);
        self.itr_hprdl.insert(dom_id, hprdl_id);
        self.itr_hsrcl.insert(dom_id, hsrcl_id);
        Ok(())
    }

    pub fn add_heats_dom(&mut self, dom_id: usize, temps_id: usize, vlcps_id: usize, conds_id: usize, hprds_id: usize, hsrcs_id: usize, aprds_id: usize, asrcs_id: usize) -> Result<(), VFVChemError> {
        // the solid volumetric heat capacity is not used in the steady solver; kept for a common API
        self.itr_temps.insert(dom_id, temps_id);
        self.itr_vlcps.insert(dom_id, vlcps_id);
        self.itr_conds.insert(dom_id, conds_id);
        self.itr_hprds.insert(dom_id, hprds_id);
        self.itr_hsrcs.insert(dom_id, hsrcs_id);
        self.itr_aprds.insert(dom_id, aprds_id);
        self.itr_asrcs.insert(dom_id, asrcs_id);
        Ok(())
    }

    pub fn add_matrix_dom(&mut self, dom_id: usize, htc_id: usize, area_id: usize, displ_id: usize, dispt_id: usize) -> Result<(), VFVChemError> {
        self.itr_htc.insert(dom_id, htc_id);
        self.itr_area.insert(dom_id, area_id);
        self.itr_displ.insert(dom_id, displ_id);
        self.itr_dispt.insert(dom_id, dispt_id);
        Ok(())
    }

    pub fn add_templ_bnd(&mut self, bnd_id: usize, templ_id: usize) -> Result<(), VFVChemError> {
        self.templ_bnd.push(bnd_id);
        self.templ_templ.insert(bnd_id, templ_id);
        Ok(())
    }

    pub fn add_hflxl_bnd(&mut self, bnd_id: usize, hflxl_id: usize) -> Result<(), VFVChemError> {
        self.hflxl_bnd.push(bnd_id);
        self.hflxl_hflxl.insert(bnd_id, hflxl_id);
        Ok(())
    }

    pub fn add_houtl_bnd(&mut self, bnd_id: usize) -> Result<(), VFVChemError> {
        self.houtl_bnd.push(bnd_id);
        Ok(())
    }

    pub fn add_hcntl_itf(&mut self, itf_id: usize) -> Result<(), VFVChemError> {
        self.hcntl_itf.push(itf_id);
        Ok(())
    }

    pub fn add_temps_bnd(&mut self, bnd_id: usize, temps_id: usize) -> Result<(), VFVChemError> {
        self.temps_bnd.push(bnd_id);
        self.temps_temps.insert(bnd_id, temps_id);
        Ok(())
    }

    pub fn add_hflxs_bnd(&mut self, bnd_id: usize, hflxs_id: usize) -> Result<(), VFVChemError> {
        self.hflxs_bnd.push(bnd_id);
        self.hflxs_hflxs.insert(bnd_id, hflxs_id);
        Ok(())
    }

    pub fn add_hcnts_itf(&mut self, itf_id: usize) -> Result<(), VFVChemError> {
        self.hcnts_itf.push(itf_id);
        Ok(())
    }

    pub fn add_pore_dom(&mut self, dom_id: usize, vel_id: usize, pres_id: usize, den_id: usize, visc_id: usize, fce_id: usize, epsl_id: usize, perm_id: usize, forch_id: usize) -> Result<(), VFVChemError> {
        self.itr_dom.push(dom_id);
        self.itr_vel.insert(dom_id, vel_id);
        self.itr_pres.insert(dom_id, pres_id);
        self.itr_den.insert(dom_id, den_id);
        self.itr_visc.insert(dom_id, visc_id);
        self.itr_fce.insert(dom_id, fce_id);
        self.itr_epsl.insert(dom_id, epsl_id);
        self.itr_perm.insert(dom_id, perm_id);
        self.itr_forch.insert(dom_id, forch_id);
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

    pub fn add_trac_bnd(&mut self, bnd_id: usize, pres_id: usize) -> Result<(), VFVChemError> {
        self.trac_bnd.push(bnd_id);
        self.trac_pres.insert(bnd_id, pres_id);
        Ok(())
    }

    pub fn add_fcnt_itf(&mut self, itf_id: usize) -> Result<(), VFVChemError> {
        self.fcnt_itf.push(itf_id);
        Ok(())
    }

    pub fn set_pres_ref(&mut self, dom_id: usize, nid: usize, pres: f64) -> Result<(), VFVChemError> {
        self.pref_dom.push(dom_id);
        self.pref_nid.insert(dom_id, nid);
        self.pref_pres.insert(dom_id, pres);
        Ok(())
    }

}

impl SteadyBase for SteadyHeatPore {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknown vectors and scalars
        let mut xid = 0;  // index in solution vector
        for (&dom_id, &vel_id) in &self.itr_vel {
            vars.vecdom[vel_id].vecdom_type = VecDomType::Unknown { start: xid };
            xid += 2 * vars.dom[dom_id].num_tot();
        }
        for (&dom_id, &pres_id) in &self.itr_pres {
            vars.scldom[pres_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }
        for (&dom_id, &templ_id) in &self.itr_templ {
            vars.scldom[templ_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }
        for (&dom_id, &temps_id) in &self.itr_temps {
            vars.scldom[temps_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }

        // create flow data
        // replace in place on a re-solve to avoid orphaned flow data
        for &dom_id in &self.itr_dom {
            let vel_id = self.itr_vel[&dom_id];
            let pres_id = self.itr_pres[&dom_id];
            let den_id = self.itr_den[&dom_id];
            let visc_id = self.itr_visc[&dom_id];
            let fce_id = self.itr_fce[&dom_id];
            match self.itr_flw.get(&dom_id) {
                Some(&flw_id) => { vars.sclflw[flw_id] = ScalarFlow::new(flw_id, &vars.dom[dom_id], vel_id, pres_id, den_id, visc_id, fce_id)?; }
                None => { self.itr_flw.insert(dom_id, vars.add_sclflw(dom_id, vel_id, pres_id, den_id, visc_id, fce_id)?); }
            }
        }

        // last xid is the matrix size
        Ok(xid)
    }

    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError> {
        // purge operators from a previous solve
        self.oper_itr_heatl.clear();
        self.oper_itr_hexcl.clear();
        self.oper_bnd_templ.clear();
        self.oper_bnd_hflxl.clear();
        self.oper_bnd_houtl.clear();
        self.oper_bnd_hexcl.clear();
        self.oper_itf_hcntl.clear();
        self.oper_itf_hadvl.clear();
        self.oper_itr_heats.clear();
        self.oper_bnd_temps.clear();
        self.oper_bnd_hflxs.clear();
        self.oper_itf_hcnts.clear();
        self.oper_itf_heats.clear();
        self.oper_itr.clear();
        self.oper_bnd_vel.clear();
        self.oper_bnd_pres.clear();
        self.oper_bnd_trac.clear();
        self.oper_itf_fcnt.clear();
        self.oper_pref.clear();

        // every fluid heat domain needs flow and matrix data; every solid heat domain needs fluid heat data
        for &dom_id in self.itr_templ.keys() {
            if !self.itr_flw.contains_key(&dom_id) || !self.itr_htc.contains_key(&dom_id) { return Err(VFVChemError::HeatporeMissingDomain { dom_id }) }
        }
        for &dom_id in self.itr_temps.keys() {
            if !self.itr_templ.contains_key(&dom_id) { return Err(VFVChemError::HeatporeMissingDomain { dom_id }) }
        }

        // fluid heat internal
        for (&dom_id, &templ_id) in &self.itr_templ {
            let vlcpl_id = self.itr_vlcpl[&dom_id];
            let condl_id = self.itr_condl[&dom_id];
            let hprdl_id = self.itr_hprdl[&dom_id];
            let hsrcl_id = self.itr_hsrcl[&dom_id];
            let vel_id = self.itr_vel[&dom_id];
            let den_id = self.itr_den[&dom_id];
            let epsl_id = self.itr_epsl[&dom_id];
            let displ_id = self.itr_displ[&dom_id];
            let dispt_id = self.itr_dispt[&dom_id];
            let flw_id = self.itr_flw[&dom_id];
            let oper_diff = PoreOpSclDomHeatDiffusion::new(vars, dom_id, templ_id, condl_id, epsl_id)?;
            let oper_disp = PoreOpSclDomDispersion::new(vars, dom_id, templ_id, vlcpl_id, vel_id, epsl_id, displ_id, dispt_id)?;
            let oper_adv = PoreOpSclDomHeatAdvection::new(vars, dom_id, templ_id, vlcpl_id, condl_id, epsl_id, den_id, vel_id, displ_id, dispt_id, flw_id)?;
            let oper_hprd = PoreOpSclDomProduction::new(vars, dom_id, templ_id, hprdl_id, epsl_id)?;
            let oper_hsrc = PoreOpSclDomSource::new(vars, dom_id, templ_id, hsrcl_id, epsl_id)?;
            self.oper_itr_heatl.push((oper_diff, oper_disp, oper_adv, oper_hprd, oper_hsrc));

            // exchange with the solid
            if let Some(&temps_id) = self.itr_temps.get(&dom_id) {
                let htc_id = self.itr_htc[&dom_id];
                let area_id = self.itr_area[&dom_id];
                self.oper_itr_hexcl.push(PoreOpSclDomHeatExchange::new(vars, dom_id, templ_id, temps_id, htc_id, area_id, epsl_id)?);
            }
        }

        // fluid heat boundary
        for &bnd_id in &self.templ_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_templ_id = self.itr_templ[&dom_id];
            let bnd_templ_id = self.templ_templ[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(vars, bnd_id, dom_templ_id, bnd_templ_id)?;
            self.oper_bnd_templ.push(oper_dir);
        }
        for &bnd_id in &self.hflxl_bnd {
            // the Neumann value is the total non-advective (conductive plus dispersive) fluid heat flux on the boundary half-edges
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_templ_id = self.itr_templ[&dom_id];
            let dom_vlcpl_id = self.itr_vlcpl[&dom_id];
            let dom_condl_id = self.itr_condl[&dom_id];
            let dom_hprdl_id = self.itr_hprdl[&dom_id];
            let dom_hsrcl_id = self.itr_hsrcl[&dom_id];
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let dom_epsl_id = self.itr_epsl[&dom_id];
            let dom_displ_id = self.itr_displ[&dom_id];
            let dom_dispt_id = self.itr_dispt[&dom_id];
            let dom_flw_id = self.itr_flw[&dom_id];
            let bnd_hflxl_id = self.hflxl_hflxl[&bnd_id];
            let oper_diff = PoreOpSclBndHeatDiffusion::new(vars, bnd_id, dom_templ_id, dom_condl_id, dom_epsl_id)?;
            let oper_disp = PoreOpSclBndDispersion::new(vars, bnd_id, dom_templ_id, dom_vlcpl_id, dom_vel_id, dom_epsl_id, dom_displ_id, dom_dispt_id)?;
            let oper_adv = PoreOpSclBndHeatAdvection::new(vars, bnd_id, dom_templ_id, dom_vlcpl_id, dom_condl_id, dom_epsl_id, dom_den_id, dom_vel_id, dom_displ_id, dom_dispt_id, dom_flw_id)?;
            let oper_hprd = PoreOpSclBndProduction::new(vars, bnd_id, dom_templ_id, dom_hprdl_id, dom_epsl_id)?;
            let oper_hsrc = PoreOpSclBndSource::new(vars, bnd_id, dom_templ_id, dom_hsrcl_id, dom_epsl_id)?;
            let oper_neu = OpSclBndNeumann::new(vars, bnd_id, dom_templ_id, bnd_hflxl_id)?;
            self.oper_bnd_hflxl.push((oper_diff, oper_disp, oper_adv, oper_hprd, oper_hsrc, oper_neu));
        }
        for &bnd_id in &self.houtl_bnd {
            // zero conductive and dispersive fluid heat flux on the boundary half-edges; heat leaves by advection
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_templ_id = self.itr_templ[&dom_id];
            let dom_vlcpl_id = self.itr_vlcpl[&dom_id];
            let dom_condl_id = self.itr_condl[&dom_id];
            let dom_hprdl_id = self.itr_hprdl[&dom_id];
            let dom_hsrcl_id = self.itr_hsrcl[&dom_id];
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let dom_epsl_id = self.itr_epsl[&dom_id];
            let dom_displ_id = self.itr_displ[&dom_id];
            let dom_dispt_id = self.itr_dispt[&dom_id];
            let dom_flw_id = self.itr_flw[&dom_id];
            let oper_diff = PoreOpSclBndHeatDiffusion::new(vars, bnd_id, dom_templ_id, dom_condl_id, dom_epsl_id)?;
            let oper_disp = PoreOpSclBndDispersion::new(vars, bnd_id, dom_templ_id, dom_vlcpl_id, dom_vel_id, dom_epsl_id, dom_displ_id, dom_dispt_id)?;
            let oper_adv = PoreOpSclBndHeatAdvection::new(vars, bnd_id, dom_templ_id, dom_vlcpl_id, dom_condl_id, dom_epsl_id, dom_den_id, dom_vel_id, dom_displ_id, dom_dispt_id, dom_flw_id)?;
            let oper_hprd = PoreOpSclBndProduction::new(vars, bnd_id, dom_templ_id, dom_hprdl_id, dom_epsl_id)?;
            let oper_hsrc = PoreOpSclBndSource::new(vars, bnd_id, dom_templ_id, dom_hsrcl_id, dom_epsl_id)?;
            self.oper_bnd_houtl.push((oper_diff, oper_disp, oper_adv, oper_hprd, oper_hsrc));
        }
        for &bnd_id in self.hflxl_bnd.iter().chain(self.houtl_bnd.iter()) {
            // exchange with the solid on the boundary nodes
            let dom_id = vars.bnd[bnd_id].dom_id;
            if let Some(&dom_temps_id) = self.itr_temps.get(&dom_id) {
                let dom_templ_id = self.itr_templ[&dom_id];
                let dom_htc_id = self.itr_htc[&dom_id];
                let dom_area_id = self.itr_area[&dom_id];
                let dom_epsl_id = self.itr_epsl[&dom_id];
                self.oper_bnd_hexcl.push(PoreOpSclBndHeatExchange::new(vars, bnd_id, dom_templ_id, dom_temps_id, dom_htc_id, dom_area_id, dom_epsl_id)?);
            }
        }

        // fluid heat interface
        // rows of domain a set continuous fluid temperature; rows of domain b balance fluid heat over both halves
        for &itf_id in &self.hcntl_itf {
            let dom_a_id = vars.itf[itf_id].dom_a_id;
            let dom_b_id = vars.itf[itf_id].dom_b_id;

            // both sides need fluid heat, flow, and matrix data
            for dom_id in [dom_a_id, dom_b_id] {
                if !self.itr_templ.contains_key(&dom_id) || !self.itr_flw.contains_key(&dom_id) || !self.itr_htc.contains_key(&dom_id) { return Err(VFVChemError::HeatporeMissingDomain { dom_id }) }
            }

            let (templ_a_id, templ_b_id) = (self.itr_templ[&dom_a_id], self.itr_templ[&dom_b_id]);
            let (condl_a_id, condl_b_id) = (self.itr_condl[&dom_a_id], self.itr_condl[&dom_b_id]);
            let (hprdl_a_id, hprdl_b_id) = (self.itr_hprdl[&dom_a_id], self.itr_hprdl[&dom_b_id]);
            let (hsrcl_a_id, hsrcl_b_id) = (self.itr_hsrcl[&dom_a_id], self.itr_hsrcl[&dom_b_id]);
            let (epsl_a_id, epsl_b_id) = (self.itr_epsl[&dom_a_id], self.itr_epsl[&dom_b_id]);
            let oper_cont = OpSclItfContinuity::new(vars, itf_id, templ_a_id, templ_b_id)?;
            let oper_diff = PoreOpSclItfHeatDiffusion::new(vars, itf_id, templ_a_id, condl_a_id, epsl_a_id, templ_b_id, condl_b_id, epsl_b_id)?;
            let oper_hprd = PoreOpSclItfProduction::new(vars, itf_id, templ_a_id, hprdl_a_id, epsl_a_id, templ_b_id, hprdl_b_id, epsl_b_id)?;
            let oper_hsrc = PoreOpSclItfSource::new(vars, itf_id, templ_a_id, hsrcl_a_id, epsl_a_id, templ_b_id, hsrcl_b_id, epsl_b_id)?;
            self.oper_itf_hcntl.push((oper_cont, oper_diff, oper_hprd, oper_hsrc));

            // dispersion and advection on each side
            for dom_id in [dom_a_id, dom_b_id] {
                let vlcpl_id = self.itr_vlcpl[&dom_id];
                let condl_id = self.itr_condl[&dom_id];
                let vel_id = self.itr_vel[&dom_id];
                let den_id = self.itr_den[&dom_id];
                let epsl_id = self.itr_epsl[&dom_id];
                let displ_id = self.itr_displ[&dom_id];
                let dispt_id = self.itr_dispt[&dom_id];
                let flw_id = self.itr_flw[&dom_id];
                let oper_disp = PoreOpSclItfDispersion::new(vars, itf_id, templ_a_id, templ_b_id, vlcpl_id, vel_id, epsl_id, displ_id, dispt_id)?;
                let oper_adv = PoreOpSclItfHeatAdvection::new(vars, itf_id, templ_a_id, templ_b_id, vlcpl_id, condl_id, epsl_id, den_id, vel_id, displ_id, dispt_id, flw_id)?;
                self.oper_itf_hadvl.push((oper_disp, oper_adv));
            }
        }

        // solid heat internal
        for (&dom_id, &temps_id) in &self.itr_temps {
            let templ_id = self.itr_templ[&dom_id];
            let conds_id = self.itr_conds[&dom_id];
            let hprds_id = self.itr_hprds[&dom_id];
            let hsrcs_id = self.itr_hsrcs[&dom_id];
            let aprds_id = self.itr_aprds[&dom_id];
            let asrcs_id = self.itr_asrcs[&dom_id];
            let epsl_id = self.itr_epsl[&dom_id];
            let htc_id = self.itr_htc[&dom_id];
            let area_id = self.itr_area[&dom_id];
            let oper_diff = PoreOpSclDomSolidDiffusion::new(vars, dom_id, temps_id, conds_id, epsl_id)?;
            let oper_hprd = PoreOpSclDomSolidProduction::new(vars, dom_id, temps_id, hprds_id, epsl_id)?;
            let oper_hsrc = PoreOpSclDomSolidSource::new(vars, dom_id, temps_id, hsrcs_id, epsl_id)?;
            let oper_aprd = PoreOpSclDomProduction::new(vars, dom_id, temps_id, aprds_id, area_id)?;
            let oper_asrc = PoreOpSclDomSource::new(vars, dom_id, temps_id, asrcs_id, area_id)?;
            let oper_hexc = PoreOpSclDomHeatExchange::new(vars, dom_id, temps_id, templ_id, htc_id, area_id, epsl_id)?;
            self.oper_itr_heats.push((oper_diff, oper_hprd, oper_hsrc, oper_aprd, oper_asrc, oper_hexc));
        }

        // solid heat boundary
        for &bnd_id in &self.temps_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temps_id = self.itr_temps[&dom_id];
            let bnd_temps_id = self.temps_temps[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(vars, bnd_id, dom_temps_id, bnd_temps_id)?;
            self.oper_bnd_temps.push(oper_dir);
        }
        for &bnd_id in &self.hflxs_bnd {
            // the Neumann value is the solid conductive heat flux on the boundary half-edges; zero for adiabatic
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_templ_id = self.itr_templ[&dom_id];
            let dom_temps_id = self.itr_temps[&dom_id];
            let dom_conds_id = self.itr_conds[&dom_id];
            let dom_hprds_id = self.itr_hprds[&dom_id];
            let dom_hsrcs_id = self.itr_hsrcs[&dom_id];
            let dom_aprds_id = self.itr_aprds[&dom_id];
            let dom_asrcs_id = self.itr_asrcs[&dom_id];
            let dom_epsl_id = self.itr_epsl[&dom_id];
            let dom_htc_id = self.itr_htc[&dom_id];
            let dom_area_id = self.itr_area[&dom_id];
            let bnd_hflxs_id = self.hflxs_hflxs[&bnd_id];
            let oper_diff = PoreOpSclBndSolidDiffusion::new(vars, bnd_id, dom_temps_id, dom_conds_id, dom_epsl_id)?;
            let oper_hprd = PoreOpSclBndSolidProduction::new(vars, bnd_id, dom_temps_id, dom_hprds_id, dom_epsl_id)?;
            let oper_hsrc = PoreOpSclBndSolidSource::new(vars, bnd_id, dom_temps_id, dom_hsrcs_id, dom_epsl_id)?;
            let oper_aprd = PoreOpSclBndProduction::new(vars, bnd_id, dom_temps_id, dom_aprds_id, dom_area_id)?;
            let oper_asrc = PoreOpSclBndSource::new(vars, bnd_id, dom_temps_id, dom_asrcs_id, dom_area_id)?;
            let oper_hexc = PoreOpSclBndHeatExchange::new(vars, bnd_id, dom_temps_id, dom_templ_id, dom_htc_id, dom_area_id, dom_epsl_id)?;
            let oper_neu = OpSclBndNeumann::new(vars, bnd_id, dom_temps_id, bnd_hflxs_id)?;
            self.oper_bnd_hflxs.push((oper_diff, oper_hprd, oper_hsrc, oper_aprd, oper_asrc, oper_hexc, oper_neu));
        }

        // solid heat interface
        // with continuity, rows of domain a set continuous solid temperature and rows of domain b balance solid heat over both halves
        // without continuity, each side with a solid balances its own half; the solid surface on the interface is adiabatic
        for &itf_id in &self.hcnts_itf {
            let dom_a_id = vars.itf[itf_id].dom_a_id;
            let dom_b_id = vars.itf[itf_id].dom_b_id;

            // both sides need solid heat data and a fluid heat interface
            for dom_id in [dom_a_id, dom_b_id] {
                if !self.itr_temps.contains_key(&dom_id) || !self.hcntl_itf.contains(&itf_id) { return Err(VFVChemError::HeatporeMissingDomain { dom_id }) }
            }

            let oper_cont = OpSclItfContinuity::new(vars, itf_id, self.itr_temps[&dom_a_id], self.itr_temps[&dom_b_id])?;
            self.oper_itf_hcnts.push(oper_cont);
        }
        for &itf_id in &self.hcntl_itf {
            let dom_a_id = vars.itf[itf_id].dom_a_id;
            let dom_b_id = vars.itf[itf_id].dom_b_id;
            let templ_b_id = self.itr_templ[&dom_b_id];
            let hcnts = self.hcnts_itf.contains(&itf_id);

            // solid and exchange on each side with a solid
            for dom_id in [dom_a_id, dom_b_id] {
                let Some(&temps_id) = self.itr_temps.get(&dom_id) else { continue };
                let row_id = if hcnts { self.itr_temps[&dom_b_id] } else { temps_id };
                let templ_id = self.itr_templ[&dom_id];
                let conds_id = self.itr_conds[&dom_id];
                let hprds_id = self.itr_hprds[&dom_id];
                let hsrcs_id = self.itr_hsrcs[&dom_id];
                let aprds_id = self.itr_aprds[&dom_id];
                let asrcs_id = self.itr_asrcs[&dom_id];
                let epsl_id = self.itr_epsl[&dom_id];
                let htc_id = self.itr_htc[&dom_id];
                let area_id = self.itr_area[&dom_id];
                let oper_diff = PoreOpSclItfSolidDiffusion::new(vars, itf_id, temps_id, row_id, conds_id, epsl_id)?;
                let oper_hprd = PoreOpSclItfSolidProduction::new(vars, itf_id, temps_id, row_id, hprds_id, epsl_id)?;
                let oper_hsrc = PoreOpSclItfSolidSource::new(vars, itf_id, temps_id, row_id, hsrcs_id, epsl_id)?;
                let oper_aprd = PoreOpSclItfSurfaceProduction::new(vars, itf_id, temps_id, row_id, aprds_id, area_id)?;
                let oper_asrc = PoreOpSclItfSurfaceSource::new(vars, itf_id, temps_id, row_id, asrcs_id, area_id)?;
                let oper_hexcs = PoreOpSclItfHeatExchange::new(vars, itf_id, temps_id, row_id, templ_id, htc_id, area_id, epsl_id)?;
                let oper_hexcl = PoreOpSclItfHeatExchange::new(vars, itf_id, templ_id, templ_b_id, temps_id, htc_id, area_id, epsl_id)?;
                self.oper_itf_heats.push((oper_diff, oper_hprd, oper_hsrc, oper_aprd, oper_asrc, oper_hexcs, oper_hexcl));
            }
        }

        // flow internal
        for &dom_id in &self.itr_dom {
            let vel_id = self.itr_vel[&dom_id];
            let pres_id = self.itr_pres[&dom_id];
            let den_id = self.itr_den[&dom_id];
            let visc_id = self.itr_visc[&dom_id];
            let fce_id = self.itr_fce[&dom_id];
            let epsl_id = self.itr_epsl[&dom_id];
            let perm_id = self.itr_perm[&dom_id];
            let forch_id = self.itr_forch[&dom_id];
            let flw_id = self.itr_flw[&dom_id];
            let oper_adv = PoreOpVecDomAdvection::new(vars, dom_id, vel_id, den_id, visc_id, epsl_id, flw_id)?;
            let oper_grad = OpVecDomGradient::new(vars, dom_id, vel_id, pres_id)?;
            let oper_diff = PoreOpVecDomDiffusion::new(vars, dom_id, vel_id, visc_id, epsl_id)?;
            let oper_dtrn = PoreOpVecDomTranspose::new(vars, dom_id, vel_id, visc_id, epsl_id)?;
            let oper_ddil = PoreOpVecDomDilatatory::new(vars, dom_id, vel_id, visc_id, epsl_id)?;
            let oper_drag = PoreOpVecDomDrag::new(vars, dom_id, vel_id, den_id, visc_id, perm_id, forch_id)?;
            let oper_src = OpVecDomSource::new(vars, dom_id, vel_id, fce_id)?;
            let oper_div = PoreOpSclDomDivergence::new(vars, dom_id, vel_id, pres_id, den_id, epsl_id)?;
            let oper_pres = OpSclDomPressureStabilizer::new(vars, dom_id, pres_id, den_id, fce_id, flw_id)?;
            self.oper_itr.push((oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_drag, oper_src, oper_div, oper_pres));
        }

        // flow boundary
        for &bnd_id in &self.vel_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_pres_id = self.itr_pres[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let dom_fce_id = self.itr_fce[&dom_id];
            let dom_epsl_id = self.itr_epsl[&dom_id];
            let dom_flw_id = self.itr_flw[&dom_id];
            let bnd_vel_id = self.vel_vel[&bnd_id];
            let oper_dir = OpVecBndDirichlet::new(vars, bnd_id, dom_vel_id, bnd_vel_id)?;
            let oper_div = PoreOpSclBndDivergence::new(vars, bnd_id, dom_vel_id, dom_pres_id, dom_den_id, dom_epsl_id)?;
            let oper_pres = OpSclBndPressureStabilizer::new(vars, bnd_id, dom_pres_id, dom_den_id, dom_fce_id, dom_flw_id)?;
            self.oper_bnd_vel.push((oper_dir, oper_div, oper_pres));
        }
        for &bnd_id in &self.pres_bnd {
            // momentum over the whole node control volumes; mu_eff * grad(vel) . n is do-nothing on the boundary half-edges
            // continuity is replaced with the boundary pressure
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_pres_id = self.itr_pres[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let dom_visc_id = self.itr_visc[&dom_id];
            let dom_fce_id = self.itr_fce[&dom_id];
            let dom_epsl_id = self.itr_epsl[&dom_id];
            let dom_perm_id = self.itr_perm[&dom_id];
            let dom_forch_id = self.itr_forch[&dom_id];
            let dom_flw_id = self.itr_flw[&dom_id];
            let bnd_pres_id = self.pres_pres[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(vars, bnd_id, dom_pres_id, bnd_pres_id)?;
            let oper_adv = PoreOpVecBndAdvection::new(vars, bnd_id, dom_vel_id, dom_den_id, dom_visc_id, dom_epsl_id, dom_flw_id)?;
            let oper_grad = OpVecBndGradient::new(vars, bnd_id, dom_vel_id, dom_pres_id)?;
            let oper_diff = PoreOpVecBndDiffusion::new(vars, bnd_id, dom_vel_id, dom_visc_id, dom_epsl_id)?;
            let oper_dtrn = PoreOpVecBndTranspose::new(vars, bnd_id, dom_vel_id, dom_visc_id, dom_epsl_id)?;
            let oper_ddil = PoreOpVecBndDilatatory::new(vars, bnd_id, dom_vel_id, dom_visc_id, dom_epsl_id)?;
            let oper_drag = PoreOpVecBndDrag::new(vars, bnd_id, dom_vel_id, dom_den_id, dom_visc_id, dom_perm_id, dom_forch_id)?;
            let oper_src = OpVecBndSource::new(vars, bnd_id, dom_vel_id, dom_fce_id)?;
            self.oper_bnd_pres.push((oper_dir, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_drag, oper_src));
        }
        for &bnd_id in &self.trac_bnd {
            // momentum over the whole node control volumes; the boundary pressure sets the traction on the boundary half-edges
            // continuity is kept, so pressure stays unknown on the boundary
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_pres_id = self.itr_pres[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let dom_visc_id = self.itr_visc[&dom_id];
            let dom_fce_id = self.itr_fce[&dom_id];
            let dom_epsl_id = self.itr_epsl[&dom_id];
            let dom_perm_id = self.itr_perm[&dom_id];
            let dom_forch_id = self.itr_forch[&dom_id];
            let dom_flw_id = self.itr_flw[&dom_id];
            let bnd_pres_id = self.trac_pres[&bnd_id];
            let oper_trac = OpVecBndTraction::new(vars, bnd_id, dom_vel_id, dom_pres_id, bnd_pres_id)?;
            let oper_adv = PoreOpVecBndAdvection::new(vars, bnd_id, dom_vel_id, dom_den_id, dom_visc_id, dom_epsl_id, dom_flw_id)?;
            let oper_diff = PoreOpVecBndDiffusion::new(vars, bnd_id, dom_vel_id, dom_visc_id, dom_epsl_id)?;
            let oper_dtrn = PoreOpVecBndTranspose::new(vars, bnd_id, dom_vel_id, dom_visc_id, dom_epsl_id)?;
            let oper_ddil = PoreOpVecBndDilatatory::new(vars, bnd_id, dom_vel_id, dom_visc_id, dom_epsl_id)?;
            let oper_drag = PoreOpVecBndDrag::new(vars, bnd_id, dom_vel_id, dom_den_id, dom_visc_id, dom_perm_id, dom_forch_id)?;
            let oper_src = OpVecBndSource::new(vars, bnd_id, dom_vel_id, dom_fce_id)?;
            let oper_div = PoreOpSclBndDivergence::new(vars, bnd_id, dom_vel_id, dom_pres_id, dom_den_id, dom_epsl_id)?;
            let oper_pres = OpSclBndPressureStabilizer::new(vars, bnd_id, dom_pres_id, dom_den_id, dom_fce_id, dom_flw_id)?;
            self.oper_bnd_trac.push((oper_trac, oper_adv, oper_diff, oper_dtrn, oper_ddil, oper_drag, oper_src, oper_div, oper_pres));
        }

        // flow interface
        // rows of domain a set continuous velocity and pressure; rows of domain b balance momentum and mass over both halves
        for &itf_id in &self.fcnt_itf {
            let dom_a_id = vars.itf[itf_id].dom_a_id;
            let dom_b_id = vars.itf[itf_id].dom_b_id;

            // both sides need flow data
            for dom_id in [dom_a_id, dom_b_id] {
                if !self.itr_flw.contains_key(&dom_id) { return Err(VFVChemError::HeatporeMissingDomain { dom_id }) }
            }

            let (vel_a_id, vel_b_id) = (self.itr_vel[&dom_a_id], self.itr_vel[&dom_b_id]);
            let (pres_a_id, pres_b_id) = (self.itr_pres[&dom_a_id], self.itr_pres[&dom_b_id]);
            let (den_a_id, den_b_id) = (self.itr_den[&dom_a_id], self.itr_den[&dom_b_id]);
            let (visc_a_id, visc_b_id) = (self.itr_visc[&dom_a_id], self.itr_visc[&dom_b_id]);
            let (fce_a_id, fce_b_id) = (self.itr_fce[&dom_a_id], self.itr_fce[&dom_b_id]);
            let (epsl_a_id, epsl_b_id) = (self.itr_epsl[&dom_a_id], self.itr_epsl[&dom_b_id]);
            let (perm_a_id, perm_b_id) = (self.itr_perm[&dom_a_id], self.itr_perm[&dom_b_id]);
            let (forch_a_id, forch_b_id) = (self.itr_forch[&dom_a_id], self.itr_forch[&dom_b_id]);
            let (flw_a_id, flw_b_id) = (self.itr_flw[&dom_a_id], self.itr_flw[&dom_b_id]);
            let oper_vcnt = OpVecItfContinuity::new(vars, itf_id, vel_a_id, vel_b_id)?;
            let oper_pcnt = OpSclItfContinuity::new(vars, itf_id, pres_a_id, pres_b_id)?;
            let oper_adv = PoreOpVecItfAdvection::new(vars, itf_id, vel_a_id, den_a_id, visc_a_id, epsl_a_id, flw_a_id, vel_b_id, den_b_id, visc_b_id, epsl_b_id, flw_b_id)?;
            let oper_grad = OpVecItfGradient::new(vars, itf_id, vel_a_id, pres_a_id, vel_b_id, pres_b_id)?;
            let oper_diff = PoreOpVecItfDiffusion::new(vars, itf_id, vel_a_id, visc_a_id, epsl_a_id, vel_b_id, visc_b_id, epsl_b_id)?;
            let oper_dtrn = PoreOpVecItfTranspose::new(vars, itf_id, vel_a_id, visc_a_id, epsl_a_id, vel_b_id, visc_b_id, epsl_b_id)?;
            let oper_ddil = PoreOpVecItfDilatatory::new(vars, itf_id, vel_a_id, visc_a_id, epsl_a_id, vel_b_id, visc_b_id, epsl_b_id)?;
            let oper_drag = PoreOpVecItfDrag::new(vars, itf_id, vel_a_id, den_a_id, visc_a_id, perm_a_id, forch_a_id, vel_b_id, den_b_id, visc_b_id, perm_b_id, forch_b_id)?;
            let oper_src = OpVecItfSource::new(vars, itf_id, vel_a_id, fce_a_id, vel_b_id, fce_b_id)?;
            let oper_div = PoreOpSclItfDivergence::new(vars, itf_id, vel_a_id, pres_a_id, den_a_id, epsl_a_id, vel_b_id, pres_b_id, den_b_id, epsl_b_id)?;
            let oper_pres = OpSclItfPressureStabilizer::new(vars, itf_id, pres_a_id, den_a_id, fce_a_id, flw_a_id, pres_b_id, den_b_id, fce_b_id, flw_b_id)?;
            self.oper_itf_fcnt.push((oper_vcnt, oper_pcnt, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_drag, oper_src, oper_div, oper_pres));
        }

        // reference pressure
        for &dom_id in &self.pref_dom {
            let pres_id = self.itr_pres[&dom_id];
            let nid = self.pref_nid[&dom_id];
            let pres = self.pref_pres[&dom_id];
            let oper_pref = OpSclDomPressureReference::new(vars, dom_id, pres_id, nid, pres)?;
            self.oper_pref.push(oper_pref);
        }

        Ok(())
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // flow internal
        for (oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_drag, oper_src, oper_div, oper_pres) in &self.oper_itr {
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_drag.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }

        // flow boundary
        for (oper_dir, oper_div, oper_pres) in &self.oper_bnd_vel {
            oper_dir.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_dir, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_drag, oper_src) in &self.oper_bnd_pres {
            oper_dir.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_drag.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_trac, oper_adv, oper_diff, oper_dtrn, oper_ddil, oper_drag, oper_src, oper_div, oper_pres) in &self.oper_bnd_trac {
            oper_trac.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_drag.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }

        // flow interface
        for (oper_vcnt, oper_pcnt, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_drag, oper_src, oper_div, oper_pres) in &self.oper_itf_fcnt {
            oper_vcnt.apply(vars, a_triplet, b_vec)?;
            oper_pcnt.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_drag.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }

        // fluid heat internal
        for (oper_diff, oper_disp, oper_adv, oper_hprd, oper_hsrc) in &self.oper_itr_heatl {
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_disp.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_hprd.apply(vars, a_triplet, b_vec)?;
            oper_hsrc.apply(vars, a_triplet, b_vec)?;
        }
        for oper_hexc in &self.oper_itr_hexcl {
            oper_hexc.apply(vars, a_triplet, b_vec)?;
        }

        // fluid heat boundary
        for oper_dir in &self.oper_bnd_templ {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_diff, oper_disp, oper_adv, oper_hprd, oper_hsrc, oper_neu) in &self.oper_bnd_hflxl {
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_disp.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_hprd.apply(vars, a_triplet, b_vec)?;
            oper_hsrc.apply(vars, a_triplet, b_vec)?;
            oper_neu.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_diff, oper_disp, oper_adv, oper_hprd, oper_hsrc) in &self.oper_bnd_houtl {
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_disp.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_hprd.apply(vars, a_triplet, b_vec)?;
            oper_hsrc.apply(vars, a_triplet, b_vec)?;
        }
        for oper_hexc in &self.oper_bnd_hexcl {
            oper_hexc.apply(vars, a_triplet, b_vec)?;
        }

        // fluid heat interface
        for (oper_cont, oper_diff, oper_hprd, oper_hsrc) in &self.oper_itf_hcntl {
            oper_cont.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_hprd.apply(vars, a_triplet, b_vec)?;
            oper_hsrc.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_disp, oper_adv) in &self.oper_itf_hadvl {
            oper_disp.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
        }

        // solid heat internal
        for (oper_diff, oper_hprd, oper_hsrc, oper_aprd, oper_asrc, oper_hexc) in &self.oper_itr_heats {
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_hprd.apply(vars, a_triplet, b_vec)?;
            oper_hsrc.apply(vars, a_triplet, b_vec)?;
            oper_aprd.apply(vars, a_triplet, b_vec)?;
            oper_asrc.apply(vars, a_triplet, b_vec)?;
            oper_hexc.apply(vars, a_triplet, b_vec)?;
        }

        // solid heat boundary
        for oper_dir in &self.oper_bnd_temps {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_diff, oper_hprd, oper_hsrc, oper_aprd, oper_asrc, oper_hexc, oper_neu) in &self.oper_bnd_hflxs {
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_hprd.apply(vars, a_triplet, b_vec)?;
            oper_hsrc.apply(vars, a_triplet, b_vec)?;
            oper_aprd.apply(vars, a_triplet, b_vec)?;
            oper_asrc.apply(vars, a_triplet, b_vec)?;
            oper_hexc.apply(vars, a_triplet, b_vec)?;
            oper_neu.apply(vars, a_triplet, b_vec)?;
        }

        // solid heat interface
        for oper_cont in &self.oper_itf_hcnts {
            oper_cont.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_diff, oper_hprd, oper_hsrc, oper_aprd, oper_asrc, oper_hexcs, oper_hexcl) in &self.oper_itf_heats {
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_hprd.apply(vars, a_triplet, b_vec)?;
            oper_hsrc.apply(vars, a_triplet, b_vec)?;
            oper_aprd.apply(vars, a_triplet, b_vec)?;
            oper_asrc.apply(vars, a_triplet, b_vec)?;
            oper_hexcs.apply(vars, a_triplet, b_vec)?;
            oper_hexcl.apply(vars, a_triplet, b_vec)?;
        }

        // reference pressure; replaces equations, so it is applied last
        for oper_pref in &self.oper_pref {
            oper_pref.apply(vars, a_triplet, b_vec)?;
        }

        // solid mask; replaces the solid equation with T_s = T_f where there is no solid
        for (&dom_id, &temps_id) in &self.itr_temps {
            pore_helper::mask_solid(vars, a_triplet, b_vec, dom_id, self.itr_epsl[&dom_id], temps_id, self.itr_templ[&dom_id])?;
        }

        Ok(())
    }

    fn initial_dirichlet(&self, vars: &mut Variables) {
        // boundary
        for oper_dir in &self.oper_bnd_templ {
            oper_dir.apply_initial(vars);
        }
        for oper_dir in &self.oper_bnd_temps {
            oper_dir.apply_initial(vars);
        }
        for (oper_dir, ..) in &self.oper_bnd_vel {
            oper_dir.apply_initial(vars);
        }
        for (oper_dir, ..) in &self.oper_bnd_pres {
            oper_dir.apply_initial(vars);
        }
    }

    fn update_flow(&self, vars: &mut Variables) {
        // Rhie-Chow face mass flux of a porous medium
        for &dom_id in &self.itr_dom {
            pore_helper::update_flow(vars, self.itr_flw[&dom_id], self.itr_epsl[&dom_id], self.itr_perm[&dom_id], self.itr_forch[&dom_id]);
        }
    }
}
