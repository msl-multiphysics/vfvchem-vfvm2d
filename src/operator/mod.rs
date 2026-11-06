// base class
pub mod prelude;
pub mod operator_base;

// scalar domain operators
pub mod opscl_dom_diff;
pub mod opscl_dom_div;
pub mod opscl_dom_lin;
pub mod opscl_dom_src;
pub mod opscl_dom_rc;
pub mod opscl_dom_time;
pub mod opscl_dom_twgt;

// scalar boundary operators
pub mod opscl_bnd_dir;
pub mod opscl_bnd_div;
pub mod opscl_bnd_neu;

// scalar interface operators
pub mod opscl_itf_cont;

// scalar point operators
pub mod opscl_pts_ref;

// vector domain operators
pub mod opvec_dom_adv;
pub mod opvec_dom_ddil;
pub mod opvec_dom_diff;
pub mod opvec_dom_dtrn;
pub mod opvec_dom_grad;
pub mod opvec_dom_src;
pub mod opvec_dom_time;

// vector boundary operators
pub mod opvec_bnd_dir;

// vector interface operators
pub mod opvec_itf_cont;
