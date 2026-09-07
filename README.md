## Thermoelastic DG

This repository implements a 2D linear thermoelastic solver using the
discontinuous Galerkin method.

### Install and check the environment

Install Julia 1.12 or newer, then run these commands from the repository root:

```sh
julia --project=. check_env.jl
```

The check activates this project, installs the packages recorded in
`Manifest.toml`, loads the local `DG` module, and runs a small exported
calculation. A successful run prints:

```text
Environment check passed: dependencies installed and DG loaded.
```

### Run the solver

After the environment check succeeds:

```sh
julia --project=. main.jl
```

The solver writes VTK snapshots, animations, and a final plot under `output/`.