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

### Run cases

`main.jl` is the repository driver. With no argument it runs every registered
case in order:

```sh
julia --project=. main.jl
```

Run one case or list the available cases:

```sh
julia --project=. main.jl case_1
julia --project=. main.jl case_2
julia --project=. main.jl --list
```

The individual case implementations are `case_1_carcione.jl` and
`case_2_carcione.jl`. Generated results are written to `output_figure6/` and
`output_figure7/`; both directories are ignored by git.
