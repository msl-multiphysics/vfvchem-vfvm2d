Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
2D Vertex-Centered Finite Volume Method Module - Source Code

Overview
--------
VFVChem is a vertex-centered finite volume method (FVM) solver for heat, mass, and/or momentum transfer problems in chemical engineering.
It supports steady-state and transient simulations with non-constant properties and multiple material domains.

Nodes are placed along boundaries and interfaces, and points are placed inside domains via Poisson-disk sampling.
Each domain is triangulated with a constrained Delaunay triangulation, and a control volume is built around every vertex.
Currently, only backward Euler time stepping is implemented, with constant or adaptive time steps.
Advective transport is discretized via the power law scheme.
Pressure-velocity coupling in the Navier-Stokes equation is handled via Rhie-Chow interpolation.

The source code is organized into the following folders:
- base      : core data structures (e.g., geometry, vertices, scalars, vectors, writers)
- linear    : linear matrix solvers (e.g., LU, GMRES)
- operator  : applies discretized conservation equations into global matrix equation
- steady    : steady-state solvers
- transient : transient solvers

The operator, steady, and transient folders are split into two subfolders:
- core   : heat, mass, and momentum transfer in single-phase materials
- extend : additional physics built on the core (e.g., porous medium flow)
The base classes (operator_base, steady_base, transient_base, transient_adaptive)
and the preludes stay at the top of each folder. The preludes re-export both
core and extend types, so `use vfvchem_vfvm2d::*` covers both.

Base
----
- error    : defines errors returned by the library
- geom_bnd : stores 1D boundary geometry
- geom_dom : stores 2D domain geometry
- geom_itf : stores 1D interface geometry; interfaces join two domains
- geom_tri : stores the triangulation and control volume data of a domain
- nds      : stores nodes on boundaries and interfaces
- pts      : stores points inside domains; placed via Poisson-disk sampling
- read_obj : imports boundary, interface, and domain polylines from an obj file
- scaling  : per-field reference scales for the Picard residual and the adaptive error norm
- scl_bnd  : stores a scalar field defined on a boundary
- scl_dom  : stores a scalar field defined on a domain
- scl_flw  : stores the Rhie-Chow face mass flux shared by the momentum and continuity equations; filled by the solvers
- scl_itf  : stores a scalar field defined on an interface
- vars     : owns all geometry, scalar, and vector data
- vec_bnd  : stores a 2D vector field defined on a boundary
- vec_dom  : stores a 2D vector field defined on a domain
- vec_itf  : stores a 2D vector field defined on an interface
- write_csv: writes scalar and vector values to CSV files
- write_vtu: writes scalar and vector values to VTU files

Vertices are either points or nodes. Within a domain, vertex ids are ordered with
points first, then nodes, i.e., vid = pid for points and vid = num_pts + nid for nodes.

Field types are registered through Variables helpers:
- *_cnst    : constant value
- *_func    : function of time, position [x, y], and unknown scalar values;
              interface functions receive time and unknown values from both domains
- *_unkn    : unknown field (solved for)

Flow data is registered through add_sclflw and is used by the momentum and
continuity operators so that both use the same face mass flux.
SteadyBase and TransientBase call update_func, update_grad, and update_flow
before each matrix assembly so that functions, gradients, and face mass fluxes
stay consistent with the unknowns. update_flow is a required method of each
solver because the Rhie-Chow coefficient depends on the momentum equation being
solved: flow solvers call flow_helper::update_flow, porous medium solvers call
pore_helper::update_flow, and solvers without flow leave it empty.

Operator
--------
Operator names use the form op + field type + location + term. For example,
opscl_dom_diff is a scalar domain diffusion operator. Files ending in unity
are variants whose weighting coefficient is fixed at one.
Domain operators apply to the whole domain, boundary operators apply to the
control volumes of boundary nodes, and interface operators apply to the
combined control volumes of interface nodes.

Helpers hold code shared by several solvers. They are not re-exported by the
prelude and are called by path (e.g., flow_helper::update_flow).

- operator_base       : defines the common operator trait and global matrix assembly helpers
- prelude             : re-exports operator types for convenient internal imports

Core operators (operator/core)
- flow_helper         : computes the Rhie-Chow coefficient and face mass flux of a free fluid
- opscl_bnd_adv       : adds weighted advection to the control volumes of boundary nodes
- opscl_bnd_advunity  : adds advection with unit weight to the control volumes of boundary nodes
- opscl_bnd_diff      : adds diffusion to the control volumes of boundary nodes
- opscl_bnd_dir       : enforces a scalar Dirichlet boundary condition
- opscl_bnd_div       : adds velocity divergence to the control volumes of boundary nodes
- opscl_bnd_neu       : applies a prescribed scalar outward flux
- opscl_bnd_prod      : adds a production term (s * c) to the control volumes of boundary nodes
- opscl_bnd_pstab     : adds the Rhie-Chow correction to the control volumes of boundary nodes
- opscl_bnd_src       : adds a source to the control volumes of boundary nodes
- opscl_bnd_time      : adds a weighted time derivative to the control volumes of boundary nodes
- opscl_bnd_timeunity : adds a time derivative with unit weight to the control volumes of boundary nodes
- opscl_dom_adv       : adds weighted advection to a scalar transport equation
- opscl_dom_advunity  : adds advection with unit weight to a scalar transport equation
- opscl_dom_diff      : adds diffusion to a scalar transport equation
- opscl_dom_div       : adds velocity divergence to the continuity equation
- opscl_dom_pref      : sets the pressure at one reference node
- opscl_dom_prod      : adds a production term (s * c) to a scalar transport equation
- opscl_dom_pstab     : adds the Rhie-Chow correction to the continuity equation
- opscl_dom_src       : adds a source to a scalar transport equation
- opscl_dom_time      : adds a weighted time derivative to a scalar transport equation
- opscl_dom_timeunity : adds a time derivative with unit weight to a scalar transport equation
- opscl_itf_adv       : adds weighted advection to the combined control volumes of interface nodes
- opscl_itf_advunity  : adds advection with unit weight to the combined control volumes of interface nodes
- opscl_itf_cont      : enforces scalar continuity between two domains
- opscl_itf_diff      : adds diffusion to the combined control volumes of interface nodes
- opscl_itf_div       : adds velocity divergence to the combined control volumes of interface nodes
- opscl_itf_prod      : adds a production term (s * c) to the combined control volumes of interface nodes
- opscl_itf_pstab     : adds the Rhie-Chow correction to the combined control volumes of interface nodes
- opscl_itf_src       : adds a source to the combined control volumes of interface nodes
- opscl_itf_time      : adds a weighted time derivative to the combined control volumes of interface nodes
- opscl_itf_timeunity : adds a time derivative with unit weight to the combined control volumes of interface nodes
- opvec_bnd_adv       : adds momentum advection to the control volumes of boundary nodes
- opvec_bnd_ddil      : adds the dilatatory viscous term to the control volumes of boundary nodes
- opvec_bnd_diff      : adds the viscous diffusion term to the control volumes of boundary nodes
- opvec_bnd_dir       : enforces a vector Dirichlet boundary condition
- opvec_bnd_dtrn      : adds the transpose viscous term to the control volumes of boundary nodes
- opvec_bnd_grad      : adds the pressure gradient to the control volumes of boundary nodes
- opvec_bnd_src       : adds a body force to the control volumes of boundary nodes
- opvec_bnd_time      : adds a momentum time derivative to the control volumes of boundary nodes
- opvec_dom_adv       : adds advection to the momentum equation
- opvec_dom_ddil      : adds the dilatatory viscous term to the momentum equation
- opvec_dom_diff      : adds the viscous diffusion term to the momentum equation
- opvec_dom_dtrn      : adds the transpose viscous term to the momentum equation
- opvec_dom_grad      : adds the pressure gradient to the momentum equation
- opvec_dom_src       : adds a body force to the momentum equation
- opvec_dom_time      : adds a time derivative to the momentum equation
- opvec_itf_adv       : adds momentum advection to the combined control volumes of interface nodes
- opvec_itf_cont      : enforces velocity continuity between two domains
- opvec_itf_ddil      : adds the dilatatory viscous term to the combined control volumes of interface nodes
- opvec_itf_diff      : adds the viscous diffusion term to the combined control volumes of interface nodes
- opvec_itf_dtrn      : adds the transpose viscous term to the combined control volumes of interface nodes
- opvec_itf_grad      : adds the pressure gradient to the combined control volumes of interface nodes
- opvec_itf_src       : adds a body force to the combined control volumes of interface nodes
- opvec_itf_time      : adds a momentum time derivative to the combined control volumes of interface nodes

Extended operators (operator/extend)
Porous medium operators solve the Darcy-Brinkman-Forchheimer equation for the
superficial velocity. They take the porosity (epsl) and, where needed, the
permeability (K) and Forchheimer coefficient (F) as user-input scalars; K and F are
computed outside the code (e.g., from the Ergun relation). The viscosity (mu) is that of the liquid;
the operators compute the effective (Brinkman) viscosity mu / epsl. The drag is decided
by K only: K = inf gives 1 / K = 0 and no drag, and K <= 0 is an error. A free (clear)
fluid is a pore domain with epsl = 1, K = inf, and F = 0.

Porous mass transfer operators solve the species equation for the concentration in the
pore fluid. The effective diffusivity is (epsl / tau) * rho * D, the mechanical dispersion
is rho * D_disp with D_disp = alpha_T |v| I + (alpha_L - alpha_T) v v^T / |v| and v = vel / epsl,
and surface terms (s_s, R_s) are multiplied by the specific surface area a_v. For a free
fluid, set tau = 1, a_v = 0, and alpha_L = alpha_T = 0.

Porous heat transfer operators solve the two-temperature (local thermal non-equilibrium) equations
for the fluid temperature T_f and the solid temperature T_s. Names ending in l act on the fluid (liquid)
and names ending in s act on the solid. The fluid conductance is epsl * k_f, the solid conductance is
(1 - epsl) * k_s, and the exchange h * a_v * (T_s - T_f) couples the two temperatures. The fluid
reuses the mass transfer operators for dispersion, production, source, and time derivative with
(rho * cp)_f in place of rho.

Porous interface operators join two pore domains (e.g., a packed bed and a free fluid) through
the combined control volumes of interface nodes. Each side uses its own properties. One-sided
operators (pore_opscl_itf_adv, pore_opscl_itf_disp, pore_opscl_itf_hadvl) are applied once per side.
The solid interface operators are one-sided with an explicit row unknown: the solid temperature
of domain b if the solid temperature is continuous, or the own solid temperature if the solid
surface is adiabatic.
- pore_helper            : computes the Rhie-Chow coefficient and face mass flux of a porous medium; effective viscosity,
                           effective diffusivity, weighted scalars (w * s), the face dispersion tensor,
                           solid-weighted scalars ((1 - epsl) * s), and the solid mask
- pore_opscl_bnd_adv     : adds porous advection with dispersion-aware power law weighting to the control volumes of boundary nodes
- pore_opscl_bnd_diff    : adds effective diffusion to the control volumes of boundary nodes
- pore_opscl_bnd_disp    : adds mechanical dispersion to the control volumes of boundary nodes
- pore_opscl_bnd_div     : adds velocity divergence and porous storage to the control volumes of boundary nodes
- pore_opscl_bnd_hadvl   : adds fluid heat advection with dispersion-aware power law weighting to the control volumes of boundary nodes
- pore_opscl_bnd_hdiffl  : adds fluid heat conduction, epsl * k_f, to the control volumes of boundary nodes
- pore_opscl_bnd_hdiffs  : adds solid heat conduction, (1 - epsl) * k_s, to the control volumes of boundary nodes
- pore_opscl_bnd_hexc    : adds fluid-solid heat exchange, h * a_v * (T_drv - T), to the control volumes of boundary nodes
- pore_opscl_bnd_hprods  : adds a solid production term ((1 - epsl) * s * T) to the control volumes of boundary nodes
- pore_opscl_bnd_hsrcs   : adds a solid source ((1 - epsl) * q) to the control volumes of boundary nodes
- pore_opscl_bnd_htimeprevs: adds a solid time derivative with the previous solid heat capacity to the control volumes of boundary nodes
- pore_opscl_bnd_prod    : adds a weighted production term (w * s * c) to the control volumes of boundary nodes
- pore_opscl_bnd_src     : adds a weighted source (w * R) to the control volumes of boundary nodes
- pore_opscl_bnd_timeprev: adds a porous time derivative with the previous fluid mass to the control volumes of boundary nodes
- pore_opscl_dom_adv     : adds porous advection with dispersion-aware power law weighting to a species equation
- pore_opscl_dom_diff    : adds effective diffusion, (epsl / tau) * rho * D, to a species equation
- pore_opscl_dom_disp    : adds mechanical dispersion to a species equation
- pore_opscl_dom_div     : adds velocity divergence and porous storage (epsl * rho) to the continuity equation
- pore_opscl_dom_hadvl   : adds fluid heat advection, (rho * cp)_f * vel . grad(T_f), with dispersion-aware power law weighting
- pore_opscl_dom_hdiffl  : adds fluid heat conduction, epsl * k_f, to the fluid heat equation
- pore_opscl_dom_hdiffs  : adds solid heat conduction, (1 - epsl) * k_s, to the solid heat equation
- pore_opscl_dom_hexc    : adds fluid-solid heat exchange, h * a_v * (T_drv - T), to the fluid or solid heat equation
- pore_opscl_dom_hprods  : adds a solid production term ((1 - epsl) * s * T) to the solid heat equation
- pore_opscl_dom_hsrcs   : adds a solid source ((1 - epsl) * q) to the solid heat equation
- pore_opscl_dom_htimeprevs: adds a solid time derivative, ((1 - epsl) * (rho * c)_s)_old * (T_s - T_s_old) / dt
- pore_opscl_dom_prod    : adds a weighted production term (w * s * c) to a species equation; w = epsl or a_v
- pore_opscl_dom_src     : adds a weighted source (w * R) to a species equation; w = epsl or a_v
- pore_opscl_dom_timeprev: adds a porous time derivative, (epsl * rho)_old * (c - c_old) / dt, to a species equation
- pore_opscl_itf_adv     : adds porous advection from one side to the combined control volumes of interface nodes
- pore_opscl_itf_aprods  : adds a surface production term (a_v * s * T) from one side to the solid rows of interface nodes
- pore_opscl_itf_asrcs   : adds a surface source (a_v * q) from one side to the solid rows of interface nodes
- pore_opscl_itf_diff    : adds effective diffusion to the combined control volumes of interface nodes
- pore_opscl_itf_disp    : adds mechanical dispersion from one side to the combined control volumes of interface nodes
- pore_opscl_itf_div     : adds velocity divergence and porous storage to the combined control volumes of interface nodes
- pore_opscl_itf_hadvl   : adds fluid heat advection from one side to the combined control volumes of interface nodes
- pore_opscl_itf_hdiffl  : adds fluid heat conduction to the combined control volumes of interface nodes
- pore_opscl_itf_hdiffs  : adds solid heat conduction from one side to the solid rows of interface nodes
- pore_opscl_itf_hexc    : adds fluid-solid heat exchange from one side to the fluid or solid rows of interface nodes
- pore_opscl_itf_hprods  : adds a solid production term from one side to the solid rows of interface nodes
- pore_opscl_itf_hsrcs   : adds a solid source from one side to the solid rows of interface nodes
- pore_opscl_itf_htimeprevs: adds a solid time derivative from one side to the solid rows of interface nodes
- pore_opscl_itf_prod    : adds a weighted production term to the combined control volumes of interface nodes
- pore_opscl_itf_src     : adds a weighted source to the combined control volumes of interface nodes
- pore_opscl_itf_timeprev: adds a porous time derivative with the previous fluid mass to the combined control volumes of interface nodes
- pore_opvec_bnd_adv  : adds porous momentum advection to the control volumes of boundary nodes
- pore_opvec_bnd_ddil : adds the dilatatory viscous term with mu_eff to the control volumes of boundary nodes
- pore_opvec_bnd_diff : adds the viscous diffusion term with mu_eff to the control volumes of boundary nodes
- pore_opvec_bnd_drag : adds Darcy and Forchheimer drag to the control volumes of boundary nodes
- pore_opvec_bnd_dtrn : adds the transpose viscous term with mu_eff to the control volumes of boundary nodes
- pore_opvec_bnd_time : adds a porous momentum time derivative (rho / epsl) to the control volumes of boundary nodes
- pore_opvec_dom_adv  : adds porous advection, (rho / epsl) * vel . grad(vel / epsl), to the momentum equation
- pore_opvec_dom_ddil : adds the dilatatory viscous term with mu_eff to the momentum equation
- pore_opvec_dom_diff : adds the viscous diffusion term with mu_eff to the momentum equation
- pore_opvec_dom_drag : adds Darcy and Forchheimer drag to the momentum equation; K and F are user inputs
- pore_opvec_dom_dtrn : adds the transpose viscous term with mu_eff to the momentum equation
- pore_opvec_dom_time : adds a porous time derivative (rho / epsl) to the momentum equation
- pore_opvec_itf_adv  : adds porous momentum advection to the combined control volumes of interface nodes
- pore_opvec_itf_ddil : adds the dilatatory viscous term with mu_eff to the combined control volumes of interface nodes
- pore_opvec_itf_diff : adds the viscous diffusion term with mu_eff to the combined control volumes of interface nodes
- pore_opvec_itf_drag : adds Darcy and Forchheimer drag to the combined control volumes of interface nodes
- pore_opvec_itf_dtrn : adds the transpose viscous term with mu_eff to the combined control volumes of interface nodes
- pore_opvec_itf_time : adds a porous momentum time derivative (rho / epsl) to the combined control volumes of interface nodes

Linear
------
- linear_base  : defines the common linear solver trait
- linear_gmres : solves sparse systems with the Jacobi-preconditioned restarted GMRES method
- linear_lu    : solves sparse systems with an LU factorization

Steady
------
- steady_base         : implements the common steady nonlinear solution loop

Core solvers (steady/core)
- steady_flow         : solves steady momentum transfer
- steady_heat         : solves steady heat transfer
- steady_heatflow     : solves coupled steady heat and momentum transfer
- steady_heatmassflow : solves coupled steady heat, mass, and momentum transfer
- steady_mass         : solves steady mass transfer
- steady_massflow     : solves coupled steady mass and momentum transfer

Extended solvers (steady/extend)
- steady_heatpore     : solves coupled steady two-temperature heat transfer and flow in porous media
- steady_masspore     : solves coupled steady mass transfer and flow in porous media
- steady_pore         : solves steady porous medium flow (Darcy-Brinkman-Forchheimer)

Transient
---------
- transient_base         : implements the common time-stepping and nonlinear solution loops
- transient_adaptive     : implements adaptive time stepping for every transient solver

Core solvers (transient/core)
- transient_flow         : solves transient momentum transfer
- transient_heat         : solves transient heat transfer
- transient_heatflow     : solves coupled transient heat and momentum transfer
- transient_heatmassflow : solves coupled transient heat, mass, and momentum transfer
- transient_mass         : solves transient mass transfer
- transient_massflow     : solves coupled transient mass and momentum transfer

Extended solvers (transient/extend)
- transient_heatpore     : solves coupled transient two-temperature heat transfer and flow in porous media
- transient_masspore     : solves coupled transient mass transfer and flow in porous media
- transient_pore         : solves transient porous medium flow (Darcy-Brinkman-Forchheimer)

Porous solvers
--------------
The pore solvers take each domain through add_pore_dom (v, p, rho, mu, f, epsl, K, F).
The masspore solvers also take add_matrix_dom (tau, a_v, alpha_L, alpha_T) and
add_mass_dom (component, c, rho * D, s_v, R_v, s_s, R_s) for every domain; volumetric
terms are weighted by epsl and surface terms by a_v. Every mass domain must be a pore
domain with matrix data.

Boundary conditions of the masspore solvers:
- add_conc_bnd : prescribed concentration
- add_mflx_bnd : prescribed total non-advective (diffusive plus dispersive) outward flux on the
                 boundary half-edges; mass also crosses the boundary by advection
- add_mout_bnd : outflow; zero diffusive and dispersive flux on the boundary half-edges, mass leaves by advection

Interfaces join two pore domains with continuous velocity, pressure, and concentration:
- add_fcnt_itf          : flow continuity (pore and masspore solvers)
- add_mcnt_itf          : mass continuity of one component (masspore solvers); use with add_fcnt_itf on the same interface
Both domains of an interface must be registered (add_pore_dom, and add_matrix_dom for the masspore solvers). A free fluid next to a porous medium is a
second pore domain with epsl = 1, K = inf, F = 0, tau = 1, a_v = 0, and alpha_L = alpha_T = 0.

The heatpore solvers solve local thermal non-equilibrium heat transfer with a fluid temperature
T_f and a solid temperature T_s. Names ending in l refer to the fluid (liquid) and names ending
in s refer to the solid (e.g., templ and temps, hcntl and hcnts); the exchange (hexc) has no suffix.
  fluid: epsl * (rho cp)_f * dT_f/dt + div( (rho cp)_f * vel * T_f ) = div( (epsl * k_f * I + K_d) * grad(T_f) )
         + h * a_v * (T_s - T_f) + epsl * (s_f * T_f + q_f)
  solid: (1 - epsl) * (rho c)_s * dT_s/dt = div( (1 - epsl) * k_s * grad(T_s) ) - h * a_v * (T_s - T_f)
         + (1 - epsl) * (s_s * T_s + q_s) + a_v * (s_a * T_s + q_a)
vel is the superficial velocity and K_d = (rho cp)_f * D_disp is the thermal dispersion tensor with the
dispersivities alpha_L and alpha_T. q_f and q_s are per unit phase volume, and q_a is per unit solid
surface area. The exchange is implicit in both temperatures.

The heatpore solvers take add_pore_dom as the pore solvers and:
- add_heatl_dom  (T_f, (rho cp)_f, k_f, s_f, q_f)            : fluid heat data
- add_heats_dom  (T_s, (rho c)_s, k_s, s_s, q_s, s_a, q_a)    : solid heat data; (rho c)_s is not used by the steady solver
- add_matrix_dom (h, a_v, alpha_L, alpha_T)                  : exchange and dispersion data
Every fluid heat domain must be a pore domain with matrix data. The solid heat domain is optional
and needs fluid heat data. A free fluid is a pore domain with epsl = 1, K = inf, F = 0, and
alpha_L = alpha_T = 0, registered with add_heatl_dom and add_matrix_dom but without add_heats_dom.
Inside a solid heat domain, the solid equation of a vertex is replaced with T_s = T_f if epsl >= 1 on
the vertex and all its neighbors; the mask is applied last and overrides solid boundary and interface
rows. Production and source coefficients given as functions of T_s (e.g., a reaction heat) are
evaluated at the previous iterate (lagged); the Picard loop converges them.

In the transient solver, the time derivatives use the heat capacity of the previous time step:
(epsl * (rho cp)_f)_old * (T_f - T_f_old) / dt for the fluid and ((1 - epsl) * (rho c)_s)_old *
(T_s - T_s_old) / dt for the solid.

Boundary conditions of the heatpore solvers are set per phase; fluxes are outward:
- add_templ_bnd : prescribed fluid temperature
- add_hflxl_bnd : prescribed total non-advective (conductive plus dispersive) outward fluid heat flux on the
                  boundary half-edges; heat also crosses the boundary by advection
- add_houtl_bnd : fluid outflow; zero conductive and dispersive flux on the boundary half-edges, heat leaves by advection
- add_temps_bnd : prescribed solid temperature
- add_hflxs_bnd : prescribed conductive outward solid heat flux on the boundary half-edges
Every boundary of a solid heat domain needs a solid boundary condition (add_hflxs_bnd with zero for adiabatic).

Interfaces join two pore domains:
- add_fcnt_itf  : flow continuity
- add_hcntl_itf : fluid temperature and total fluid heat flux continuity; use with add_fcnt_itf on the same interface
- add_hcnts_itf : solid temperature and solid conductive flux continuity; use with add_hcntl_itf on the same interface
Both domains of a fluid heat interface must be registered (add_pore_dom, add_matrix_dom, and add_heatl_dom).
The solid temperature is joined only with add_hcnts_itf, which needs a solid heat domain on both sides.
Without add_hcnts_itf, each side with a solid heat domain has an adiabatic solid surface at the interface,
and a side without a solid heat domain (e.g., a free fluid) has no solid and no exchange there.
