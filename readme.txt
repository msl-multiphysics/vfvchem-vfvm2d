__     _________     ______ _                    
\ \   / /  ___\ \   / / ___| |__   ___ _ __ ___  
 \ \ / /| |_   \ \ / / |   | '_ \ / _ \ '_ ` _ \ 
  \ V / |  _|   \ V /| |___| | | |  __/ | | | | |
   \_/  |_|      \_/  \____|_| |_|\___|_| |_| |_|

Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
2D Vertex-Centered Finite Volume Method Module
Copyright (c) 2026 VFVChem Development Team

Overview
--------
VFVChem is a vertex-centered finite volume method (FVM) solver for heat, mass, and/or momentum transfer problems in chemical engineering.
It supports steady-state and transient simulations with non-constant properties and multiple material domains.
It is a pure Rust library. Simulations are configured as Rust code, and no separate installation or build system is required.

Note: The project is under active development. The public API and numerical models may change between releases.

Features
--------
- Steady or transient simulations of heat, mass, or momentum transport and couplings thereof
- Non-constant inputs, such as temperature and concentration-dependent properties and reactions
- Multicomponent mass transfer with cross-diffusion
- Multiple domains, boundaries, and interfaces, which can be imported from `.obj` geometry files (e.g., from VFVChem geom2d)
- Automatic placement of vertices via Poisson-disk sampling; no separate mesh generation is required
- Control volumes constructed from a constrained Delaunay triangulation of the vertices
- Power law scheme for advective transport
- Rhie-Chow interpolation for pressure-velocity coupling in the Navier-Stokes equation
- Direct LU linear solver and an iterative Jacobi-preconditioned GMRES solver
- Outputs simulation results as CSV or VTU files for post-processing

Requirements
------------
- A Rust toolchain with Rust 2024 edition support
- VFVChem geom2d (Python) for custom geometries (optional)
- ParaView for visualizing output (optional)

Quick Start
-----------
Clone or download the project, change to the vfvchem_vfvm2d directory, and run `cargo run --release --example quickstart`.
This example solves a steady heat transfer problem on a square domain and writes the temperature field to `examples/output_quickstart/temp_0000.vtu`.
Other examples can be found in the `examples` folder. An overview of the examples is provided in `examples/readme.txt`.

Project Structure
-----------------
A guide to the source code is provided in `src/readme.txt`.

Contributing
------------
Due to limited maintenance capacity, we are not currently accepting external contributions.
However, users are welcome to fork and modify the repository in accordance with the MIT License.

License
-------
VFVChem is available under the MIT License. See license.txt for details.
