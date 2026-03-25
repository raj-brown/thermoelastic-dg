using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using StartUpDG
using OrdinaryDiffEq
using LinearAlgebra
using OrdinaryDiffEq
using Plots

include("dg.jl")

println("running simulation")
# Set up the problem

gaussian(x, y, s) = exp(-(x^2 + y^2) / (2*s*s))

N=2
K1D = 32

rd = RefElemData(Tri(), N)
VXY, EToV = uniform_mesh(Tri(), K1D)
md = MeshData(VXY, EToV, rd)
md = make_periodic(md)


(; x, y) = md
u = zeros(rd.Np, md.K, 5)
u[:, :, 5] = gaussian.(x, y, 0.01) 

xp, yp = vec(rd.Vp * x), vec(rd.Vp * y)  

scatter(xp, yp, zcolor=vec(rd.Vp*u[:,:,5]), markersize=2, markerstrokewidth=0, legend=false, clims=(-.25, .25))


params = (; rd, md)
du = similar(u)
rhs!(du, u, params, 0.0)
tspan = (0.0, 1.0)
ode = ODEProblem(rhs!, u, tspan, params)

sol = solve(ode, Tsit5(), saveat=LinRange(tspan..., 10), 
            abstol=1e-14, reltol=1e-7)

xp, yp = vec(rd.Vp * x), vec(rd.Vp * y)   
@gif for i in eachindex(sol.u)    
    scatter(xp, yp, zcolor=vec(rd.Vp * sol.u[i][:,:,5]), 
            markersize=2, markerstrokewidth=0, legend=false, clims=(-.25, .25))
end