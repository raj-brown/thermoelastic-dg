using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using StartUpDG
using OrdinaryDiffEq
using LinearAlgebra
using OrdinaryDiffEq

include("dg.jl")

println("running simulation")
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

sol = solve(ode, Tsit5())