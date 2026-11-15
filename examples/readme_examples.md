## Running Examples
Change to the `vfvchem_vfvm2d` directory, and run any of the examples (e.g., `name_of_file.rs`):
```sh
cargo run --release --example name_of_file
```

## Modifying Examples
To create your own examples, you may create one in `examples` (e.g., `examples/my_file.rs`) and run:
```sh
cargo run --release --example my_file
```

If you create an example within a subfolder (e.g., `examples/my_physics/my_file.rs`), you must register it first in `Cargo.toml`. Add the following in `Cargo.toml`:
```
[[example]]
name = "my_file"
path = "examples/my_physics/my_file.rs"
```

After registering your example, run it as usual:
```sh
cargo run --release --example my_file
```

## List of Examples
### Heat Transfer
| Example | Description |
|---------|-------------|
| `quickstart` | Steady heat transfer in a square domain |
| `heat_transient` | Transient heat transfer in a square domain |
| `heat_diag` | Steady heat transfer in a divided domain |
| `heat_adapt` | Transient heat transfer with adaptive time stepping |

### Mass Transfer
| Example | Description |
|---------|-------------|
| `mole_turing` | Turing pattern formation (Brusselator) in a square |

### Momentum Transfer
| Example | Description |
|---------|-------------|
| `flow_lid` | Steady lid-driven cavity flow in a square |
| `flow_step` | Steady flow past a backward-facing step |
| `flow_diag` | Steady lid-driven cavity flow in a divided domain |

### Heat and Momentum Transfer
| Example | Description |
|---------|-------------|
| `heatflow_channel` | Steady heat and momentum transfer through a channel |
| `heatflow_step` | Steady heat and momentum transfer past a backward-facing step |

### Mass and Momentum Transfer
| Example | Description |
|---------|-------------|
| `moleflow_bouy` | Transient buoyancy-driven flow in a square cavity |
| `moleflow_step` | Steady reactive mass and momentum transfer past a backward-facing step |

### Heat, Mass, and Momentum Transfer
| Example | Description |
|---------|-------------|
| `heatmoleflow_step` | Transient exothermic reacting flow past a backward-facing step |

### Benchmarks
| Example | Description |
|---------|-------------|
| `bench_vert` | Composite-wall heat conduction with discontinuous conductivity |
| `bench_heat` | Heat conduction with non-constant thermal conductivity |
| `bench_flow` | Poiseuille flow with non-constant viscosity |
| `bench_cavity` | Natural convection in a differentially heated cavity (de Vahl Davis) |
| `bench_graetz` | Graetz problem: thermal entry length of plane Poiseuille flow |
| `bench_tayloraris` | Taylor-Aris dispersion of a solute pulse in plane Poiseuille flow |
| `bench_vortex` | Schäfer-Turek flow around a cylinder (Re = 100, vortex shedding) |
