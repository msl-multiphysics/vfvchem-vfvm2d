// base class
pub mod prelude;
pub mod operator_base;

// scalar domain operators
pub mod opscl_dom_diff;
pub mod opscl_dom_div;
pub mod opscl_dom_pref;
pub mod opscl_dom_pstab;
pub mod opscl_dom_src;

// scalar boundary operators
pub mod opscl_bnd_diff;
pub mod opscl_bnd_dir;
pub mod opscl_bnd_div;
pub mod opscl_bnd_neu;
pub mod opscl_bnd_pstab;
pub mod opscl_bnd_src;

// scalar interface operators

// vector domain operators
pub mod opvec_dom_adv;
pub mod opvec_dom_ddil;
pub mod opvec_dom_diff;
pub mod opvec_dom_dtrn;
pub mod opvec_dom_grad;
pub mod opvec_dom_src;

// vector boundary operators
pub mod opvec_bnd_adv;
pub mod opvec_bnd_ddil;
pub mod opvec_bnd_diff;
pub mod opvec_bnd_dir;
pub mod opvec_bnd_dtrn;
pub mod opvec_bnd_grad;
pub mod opvec_bnd_src;

// vector interface operators
