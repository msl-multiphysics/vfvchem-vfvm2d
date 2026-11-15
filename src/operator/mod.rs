// base class
pub mod prelude;
pub mod operator_base;

// generic scalar operators
pub mod scl_dom_lin;
pub mod scl_dom_src;
pub mod scl_bnd_dir;
pub mod scl_bnd_neu;
pub mod scl_bnd_out;
pub mod scl_itf_cont;

// generic vector operators
pub mod vec_dom_src;
pub mod vec_bnd_dir;
pub mod vec_itf_cont;

// heat transfer operators
pub mod heat_dom_adv;
pub mod heat_dom_diff;
pub mod heat_dom_time;

// mass transfer operators
pub mod mole_dom_adv;
pub mod mole_dom_diff;
pub mod mole_dom_time;

// momentum transfer operators
pub mod flow_dom_adv;
pub mod flow_dom_ddil;
pub mod flow_dom_diff;
pub mod flow_dom_div;
pub mod flow_dom_dtrn;
pub mod flow_dom_grad;
pub mod flow_dom_rc;
pub mod flow_dom_tcon;
pub mod flow_dom_time;
pub mod flow_bnd_ddil;
pub mod flow_bnd_div;
pub mod flow_bnd_dtrn;
pub mod flow_nds_ref;
