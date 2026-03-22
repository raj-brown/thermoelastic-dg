using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using StartUpDG
using OrdinaryDiffEq
using LinearAlgebra
using OrdinaryDiffEq

include("dg.jl")

println("running simulation")
# Set up the problem
μ = 1.0
λ = 1.0
c11 = λ + 2μ
c12 = λ
c22 = λ + 2μ


N=3
K1D = 8

rd = RefElemData(Tri(), N)
VXY, EToV = uniform_mesh(Tri(), K1D)
md = MeshData(VXY, EToV, rd)

(; x, y) = md
u = initial_condition(rd.Np, md.K)
params = (; rd, md)
du = similar(u)
rhs!(du, u, params, 0.0)
tspan = (0.0, 3.0)
ode = ODEProblem(rhs!, u, tspan, params)

sol = solve(ode, Tsit5(), saveat=LinRange(tspan..., 10), 
            abstol=1e-10, reltol=1e-10)