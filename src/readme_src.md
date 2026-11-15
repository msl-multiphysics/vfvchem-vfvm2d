## Numerical Formulation
### Geometry
Domains are discretized into vertices, which are connected via a simple Delaunay triangulation. Edges outside the domain are then removed. Vertices are positioned so that boundaries and interfaces become Delaunay edges without needing constrained Delaunay triangulation (CDT). This would be convenient later on when extending the code to 3D where CDT is a lot more difficult.

Vertices within a domain are referred to as "points" and those along a boundary or interface are referred to as "nodes." Nodes are not placed on corners to simplify implementation (i.e., no ad-hoc rules about which boundary condition applies at the corner). Instead, two nodes are placed with a very small offset from the corner.

### Advection
Advection operators use a non-conservative formulation for numerical stability and ease of implementation. Non-conservative formulations enable the direct solution of primary variables (e.g., temperature instead of enthalpy), which simplifies the handling of non-constant properties. It also preserves uniform transported fields even when the current mass fluxes do not exactly satisfy the discrete continuity equation, which may occur during nonlinear iterations.

Highly advective flows are stabilized using the Venkatakrishnan limiter. This scheme applies naturally to unstructured meshes. It also has less numerical damping than other schemes, which allows for the use of coarser meshes.

### Pressure Stabilization
A modified Rhie-Chow correction is used to stabilize the velocity-pressure coupling. In transient problems, the Rhie–Chow coefficient is formed only from the steady part of the momentum diagonal. This simplifies the implementation by not needing to store corrections from the previous time step. In addition, scalar advection operators estimate fluxes by interpolating velocity instead of the Rhie-Chow corrected flux.

### Interfaces
Nodes along an interface implement an unknown equality and a control volume balance between the adjacent domains. After the domain equations are assembled, one balance row is replaced by an equation equating the two primary variables. The original balance rows are then summed into the remaining row, which enforces continuity of the normal flux.

This design affects the implementation of boundary terms. For divergence terms (e.g., viscosity and continuity), the balance rows are written without their boundary terms $F_a$ and $F_b$ so that summing these rows yields the correct balance. The exact boundary terms satisfy $F_a + F_b = 0$, but this is not guaranteed if $F_a$ or $F_b$ are approximated using one-sided properties, thus they are omitted on interfaces and applied separately on boundaries. Likewise, for volumetric terms (e.g., sources and non-conservative advection), the balance rows are written in full along with any boundary terms that arise.

## Software Architecture
The source code is organized into the following folders:
- `base` - core data structures (e.g, geometry, fields, I/O, etc.)
- `linear` - linear matrix solvers
- `operator` - discretized PDE terms assembled into the global matrix equation
- `steady` - steady-state solvers
- `transient` - transient solvers

### `base` Code
| Name | Description |
|------|-------------|
| `error` | Library error types |
| `geom_bnd` | Boundary geometry |
| `geom_dom` | Domain geometry |
| `geom_itf` | Interface geometry joining two domains |
| `geom_tri` | Triangulation and Voronoi control-volume data of a domain |
| `nds` | Nodes along boundaries and interfaces |
| `pts` | Interior points of a domain |
| `scl_bnd` | Scalar field defined on a boundary |
| `scl_dom` | Scalar field defined on a domain |
| `scl_itf` | Scalar field defined on an interface |
| `vec_bnd` | Vector field defined on a boundary |
| `vec_dom` | Vector field defined on a domain |
| `vec_itf` | Vector field defined on an interface |
| `vars` | Owns all geometry and field data for a simulation |
| `read_obj` | Imports boundary, interface, and domain polylines from an OBJ file |
| `write_vtu` | Writes scalar and vector fields to VTU files |

### `linear` Code
| Name | Description |
|------|-------------|
| `linear_base` | Common linear solver trait |
| `rslab_lu` | Sparse LU factorization through the `rslab` crate |

### `operator` Code
| Name | Description |
|------|-------------|
| `operator_base` | Common operator trait, matrix assembly helpers, and input checks |
| `scl_dom_lin` | Adds a linear source term to a scalar transport equation |
| `scl_dom_src` | Adds a source term to a scalar transport equation |
| `scl_bnd_dir` | Enforces a scalar Dirichlet boundary condition |
| `scl_bnd_neu` | Applies a prescribed scalar outward flux |
| `scl_bnd_out` | Marks an outflow boundary (zero normal diffusive flux) |
| `scl_itf_cont` | Enforces scalar value and flux continuity across an interface |
| `vec_dom_src` | Adds a body force to the momentum equation |
| `vec_bnd_dir` | Enforces a vector Dirichlet boundary condition |
| `vec_itf_cont` | Enforces vector value and flux continuity across an interface |
| `heat_dom_adv` | Adds non-conservative advection to the heat equation |
| `heat_dom_diff` | Adds conduction to the heat equation |
| `heat_dom_time` | Adds the heat capacity time derivative |
| `mole_dom_adv` | Adds non-conservative advection to a species equation |
| `mole_dom_diff` | Adds diffusion to a species equation |
| `mole_dom_time` | Adds the species time derivative |
| `flow_dom_adv` | Adds non-conservative advection to the momentum equation |
| `flow_dom_diff` | Adds viscous diffusion to the momentum equation |
| `flow_dom_dtrn` | Adds the transpose viscous term |
| `flow_dom_ddil` | Adds the dilatatory viscous term |
| `flow_dom_grad` | Adds the pressure gradient to the momentum equation |
| `flow_dom_div` | Adds velocity divergence to the continuity equation |
| `flow_dom_rc` | Adds the Rhie-Chow correction to the continuity equation |
| `flow_dom_time` | Adds the momentum time derivative |
| `flow_dom_tcon` | Adds the density time derivative to the continuity equation |
| `flow_bnd_div` | Adds the boundary contribution of velocity divergence |
| `flow_bnd_dtrn` | Adds the boundary contribution of the transpose viscous term |
| `flow_bnd_ddil` | Adds the boundary contribution of the dilatatory viscous term |
| `flow_nds_ref` | Fixes pressure at a reference node |

### `steady` Code
| Name | Description |
|------|-------------|
| `steady_base` | Common steady nonlinear solution loop |
| `steady_heat` | Steady heat transfer |
| `steady_mole` | Steady mass transfer |
| `steady_flow` | Steady momentum transfer |
| `steady_heatflow` | Coupled steady heat and momentum transfer |
| `steady_moleflow` | Coupled steady mass and momentum transfer |
| `steady_heatmoleflow` | Coupled steady heat, mass, and momentum transfer |

### `transient` Code
| Name | Description |
|------|-------------|
| `transient_base` | Common time-stepping and nonlinear solution loops |
| `transient_adapt` | Adaptive time stepping shared by all transient solvers |
| `transient_heat` | Transient heat transfer |
| `transient_mole` | Transient mass transfer |
| `transient_flow` | Transient momentum transfer |
| `transient_heatflow` | Coupled transient heat and momentum transfer |
| `transient_moleflow` | Coupled transient mass and momentum transfer |
| `transient_heatmoleflow` | Coupled transient heat, mass, and momentum transfer |
