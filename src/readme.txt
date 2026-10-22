__     _________     ______ _                    
\ \   / /  ___\ \   / / ___| |__   ___ _ __ ___  
 \ \ / /| |_   \ \ / / |   | '_ \ / _ \ '_ ` _ \ 
  \ V / |  _|   \ V /| |___| | | |  __/ | | | | |
   \_/  |_|      \_/  \____|_| |_|\___|_| |_| |_|

Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
2D Vertex-Centered Finite Volume Method Module - Source Code

Overview
--------
VFVChem is a vertex-centered finite volume method (FVM) solver for heat, mass, and/or momentum transfer problems in chemical engineering.
It supports steady-state and transient simulations with non-constant properties and multiple material domains.

Nodes are placed along boundaries and interfaces, and points are placed inside domains via Poisson-disk sampling.
Each domain is triangulated with a constrained Delaunay triangulation, and a control volume is built around every vertex.
Currently, only backward Euler time stepping with constant time steps is implemented.
Advective transport is discretized via the power law scheme.
Pressure-velocity coupling in the Navier-Stokes equation is handled via Rhie-Chow interpolation.

The source code is organized into the following folders:
- base      : core data structures (e.g., geometry, vertices, scalars, vectors, writers)
- linear    : linear matrix solvers (e.g., LU, GMRES)
- operator  : applies discretized conservation equations into global matrix equation
- steady    : steady-state solvers
- transient : transient solvers

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
- scl_bnd  : stores a scalar field defined on a boundary
- scl_dom  : stores a scalar field defined on a domain
- scl_flw  : stores the Rhie-Chow face mass flux shared by the momentum and continuity equations
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
stay consistent with the unknowns.

Operator
--------
Operator names use the form op + field type + location + term. For example,
opscl_dom_diff is a scalar domain diffusion operator. Files ending in unity
are variants whose weighting coefficient is fixed at one.
Domain operators apply to the whole domain, boundary operators apply to the
control volumes of boundary nodes, and interface operators apply to the
combined control volumes of interface nodes.

- operator_base       : defines the common operator trait and global matrix assembly helpers
- opscl_bnd_adv       : adds weighted advection to the control volumes of boundary nodes
- opscl_bnd_advunity  : adds advection with unit weight to the control volumes of boundary nodes
- opscl_bnd_diff      : adds diffusion to the control volumes of boundary nodes
- opscl_bnd_dir       : enforces a scalar Dirichlet boundary condition
- opscl_bnd_div       : adds velocity divergence to the control volumes of boundary nodes
- opscl_bnd_neu       : applies a prescribed scalar outward flux
- opscl_bnd_pstab     : adds the Rhie-Chow correction to the control volumes of boundary nodes
- opscl_bnd_src       : adds a source to the control volumes of boundary nodes
- opscl_bnd_time      : adds a weighted time derivative to the control volumes of boundary nodes
- opscl_bnd_timeunity : adds a time derivative with unit weight to the control volumes of boundary nodes
- opscl_dom_adv       : adds weighted advection to a scalar transport equation
- opscl_dom_advunity  : adds advection with unit weight to a scalar transport equation
- opscl_dom_diff      : adds diffusion to a scalar transport equation
- opscl_dom_div       : adds velocity divergence to the continuity equation
- opscl_dom_pref      : sets the pressure at one reference node
- opscl_dom_pstab     : adds the Rhie-Chow correction to the continuity equation
- opscl_dom_src       : adds a source to a scalar transport equation
- opscl_dom_time      : adds a weighted time derivative to a scalar transport equation
- opscl_dom_timeunity : adds a time derivative with unit weight to a scalar transport equation
- opscl_itf_adv       : adds weighted advection to the combined control volumes of interface nodes
- opscl_itf_advunity  : adds advection with unit weight to the combined control volumes of interface nodes
- opscl_itf_cont      : enforces scalar continuity between two domains
- opscl_itf_diff      : adds diffusion to the combined control volumes of interface nodes
- opscl_itf_div       : adds velocity divergence to the combined control volumes of interface nodes
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
- prelude             : re-exports operator types for convenient internal imports

Linear
------
- linear_base  : defines the common linear solver trait
- linear_gmres : solves sparse systems with the Jacobi-preconditioned restarted GMRES method
- linear_lu    : solves sparse systems with an LU factorization

Steady
------
- steady_base         : implements the common steady nonlinear solution loop
- steady_flow         : solves steady momentum transfer
- steady_heat         : solves steady heat transfer
- steady_heatflow     : solves coupled steady heat and momentum transfer
- steady_heatmassflow : solves coupled steady heat, mass, and momentum transfer
- steady_mass         : solves steady mass transfer
- steady_massflow     : solves coupled steady mass and momentum transfer

Transient
---------
- transient_base         : implements the common time-stepping and nonlinear solution loops
- transient_flow         : solves transient momentum transfer
- transient_heat         : solves transient heat transfer
- transient_heatflow     : solves coupled transient heat and momentum transfer
- transient_heatmassflow : solves coupled transient heat, mass, and momentum transfer
- transient_mass         : solves transient mass transfer
- transient_massflow     : solves coupled transient mass and momentum transfer
