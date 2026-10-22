__     _________     ______ _                    
\ \   / /  ___\ \   / / ___| |__   ___ _ __ ___  
 \ \ / /| |_   \ \ / / |   | '_ \ / _ \ '_ ` _ \ 
  \ V / |  _|   \ V /| |___| | | |  __/ | | | | |
   \_/  |_|      \_/  \____|_| |_|\___|_| |_| |_|

Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
2D Vertex-Centered Finite Volume Method Module - Examples

How to Use Examples
-------------------
First-time users are encouraged to run quickstart.rs

Change the directory to vfvchem_vfvm2d then run quickstart.rs:
`cargo run --release --example quickstart`

The other examples are registered in `Cargo.toml` and can be run in place:
`cargo run --release --example name_of_file`

Geometries used by the examples are stored as `.obj` files in the `examples/obj` folder.
The Python scripts in the same folder generate these files with VFVChem geom2d.

List of Examples
----------------
Below is a list of examples and the concepts they introduce:

Heat Transfer
- quickstart.rs     : introduction to VFVChem code
- heat_transient.rs : transient heat transfer simulation
- heat_tworeg.rs    : importing `.obj` geometries; multiple domains joined by an interface

Mass Transfer
- mass_square.rs : steady-state multicomponent diffusion with cross-diffusion

Momentum Transfer
- flow_lid.rs    : lid-driven cavity flow
- flow_step.rs   : transient flow past a backward-facing step
- flow_tworeg.rs : lid-driven cavity flow through two domains joined by an interface

Heat and Momentum Transfer
- heatflow_fourcyl.rs : flow through a channel with four heated cylinders
- heatflow_heater.rs  : heat transfer from a solid to a moving fluid
- heatflow_step.rs    : transient flow past a backward-facing step with heat transfer

Mass and Momentum Transfer
- massflow_step.rs : transient reactive flow past a backward-facing step

Heat, Mass, and Momentum Transfer
- heatmassflow_step.rs : transient flow past a backward-facing step with an exothermic reaction
