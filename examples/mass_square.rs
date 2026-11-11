use vfvchem_vfvm2d::*;
use std::collections::HashMap;

// cargo run --release --example mass_square
// Steady diffusion with the analytical solution conc = 1 - x.
fn main() -> Result<(), VFVChemError> {
    // variables
    let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, 0.02, 0.05, "examples/output_mass_square")?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknowns
    let conc = vars.add_scldom_unkn(dom, 0.5, "conc")?;

    // constant domain values for the diffusion benchmark
    let den = vars.add_scldom_cnst(dom, 1.0, "den")?;
    let diff = vars.add_scldom_cnst(dom, 1.0, "diff")?;
    let nlin = vars.add_scldom_cnst(dom, 0.0, "nlin")?;
    let nsrc = vars.add_scldom_cnst(dom, 0.0, "nsrc")?;  // positive if source; negative if sink

    // boundary values
    let conc_l = vars.add_sclbnd_cnst(bnd_l, 1.0)?;
    let conc_r = vars.add_sclbnd_cnst(bnd_r, 0.0)?;
    let nflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0)?;
    let nflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0)?;  // positive if outward; negative if inward

    // physics
    let comp = 0;
    let mut phys = SteadyMole::new(1)?;
    phys.add_mole_dom(comp, dom, conc, den, HashMap::from([(comp, diff)]), nlin, nsrc)?;  // arguments: component, domain, c, rho, D, L, R
    phys.add_conc_bnd(comp, bnd_l, conc_l)?;  // arguments: component, boundary, c
    phys.add_conc_bnd(comp, bnd_r, conc_r)?;  // arguments: component, boundary, c
    phys.add_hflx_bnd(comp, bnd_b, nflx_b)?;  // arguments: component, boundary, N
    phys.add_hflx_bnd(comp, bnd_t, nflx_t)?;  // arguments: component, boundary, N

    // solver
    let mut lin = RslabLu::new(4)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-3, 1.0)?;

    Ok(())
}
