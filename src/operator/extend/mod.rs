// helpers; not in prelude, called by path
pub mod pore_helper;

// porous scalar domain operators
pub mod pore_opscl_dom_div;

// porous scalar boundary operators
pub mod pore_opscl_bnd_div;

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
