// helpers; not in prelude, called by path
pub mod pore_helper;

// porous scalar domain operators
pub mod pore_opscl_dom_adv;
pub mod pore_opscl_dom_diff;
pub mod pore_opscl_dom_disp;
pub mod pore_opscl_dom_div;
pub mod pore_opscl_dom_hadvl;
pub mod pore_opscl_dom_hdiffl;
pub mod pore_opscl_dom_hdiffs;
pub mod pore_opscl_dom_hexc;
pub mod pore_opscl_dom_hprods;
pub mod pore_opscl_dom_hsrcs;
pub mod pore_opscl_dom_htimeprevs;
pub mod pore_opscl_dom_prod;
pub mod pore_opscl_dom_src;
pub mod pore_opscl_dom_timeprev;

// porous scalar boundary operators
pub mod pore_opscl_bnd_adv;
pub mod pore_opscl_bnd_diff;
pub mod pore_opscl_bnd_disp;
pub mod pore_opscl_bnd_div;
pub mod pore_opscl_bnd_hadvl;
pub mod pore_opscl_bnd_hdiffl;
pub mod pore_opscl_bnd_hdiffs;
pub mod pore_opscl_bnd_hexc;
pub mod pore_opscl_bnd_hprods;
pub mod pore_opscl_bnd_hsrcs;
pub mod pore_opscl_bnd_htimeprevs;
pub mod pore_opscl_bnd_prod;
pub mod pore_opscl_bnd_src;
pub mod pore_opscl_bnd_timeprev;

// porous scalar interface operators
pub mod pore_opscl_itf_adv;
pub mod pore_opscl_itf_aprods;
pub mod pore_opscl_itf_asrcs;
pub mod pore_opscl_itf_diff;
pub mod pore_opscl_itf_disp;
pub mod pore_opscl_itf_div;
pub mod pore_opscl_itf_hadvl;
pub mod pore_opscl_itf_hdiffl;
pub mod pore_opscl_itf_hdiffs;
pub mod pore_opscl_itf_hexc;
pub mod pore_opscl_itf_hprods;
pub mod pore_opscl_itf_hsrcs;
pub mod pore_opscl_itf_htimeprevs;
pub mod pore_opscl_itf_prod;
pub mod pore_opscl_itf_src;
pub mod pore_opscl_itf_timeprev;

// porous vector domain operators
pub mod pore_opvec_dom_adv;
pub mod pore_opvec_dom_ddil;
pub mod pore_opvec_dom_diff;
pub mod pore_opvec_dom_drag;
pub mod pore_opvec_dom_dtrn;
pub mod pore_opvec_dom_time;

// porous vector boundary operators
pub mod pore_opvec_bnd_adv;
pub mod pore_opvec_bnd_ddil;
pub mod pore_opvec_bnd_diff;
pub mod pore_opvec_bnd_drag;
pub mod pore_opvec_bnd_dtrn;
pub mod pore_opvec_bnd_time;

// porous vector interface operators
pub mod pore_opvec_itf_adv;
pub mod pore_opvec_itf_ddil;
pub mod pore_opvec_itf_diff;
pub mod pore_opvec_itf_drag;
pub mod pore_opvec_itf_dtrn;
pub mod pore_opvec_itf_time;
