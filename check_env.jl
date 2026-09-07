using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

push!(LOAD_PATH, joinpath(@__DIR__, "src"))
using DG

@assert isfinite(compute_cfl(1.0, 0.01, (; J = [1.0]), (; N = 1)))
println("Environment check passed: dependencies installed and DG loaded.")