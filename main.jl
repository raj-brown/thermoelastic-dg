using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using StartUpDG
using OrdinaryDiffEq
using LinearAlgebra
using OrdinaryDiffEq
using Plots
using Revise
using Logging
using TerminalLoggers

global_logger(TerminalLogger())


include("dg.jl")

println("running simulation")
# Set up the problem

#gaussian(x, y, s) = exp(-(x^2 + y^2) / (2 * s * s))

N = 4
K1D = 32

rd = RefElemData(Quad(), N)
VXY, EToV = uniform_mesh(Quad(), K1D)
md = MeshData(VXY, EToV, rd)
md = make_periodic(md)


(; x, y) = md
u = zeros(rd.Np, md.K, 5)
#u[:, :, 5] = gaussian.(x, y, 0.01)

pt_src = zeros(rd.Np, md.K)
x0 = 0.0
y0 = 0.0
ids = findall(abs.(x[:] .- x0) .+ abs.(y[:] .- y0) .< 1e-8)
pt_src[ids] .= 1.0

pt_src = (rd.VDM * rd.VDM') * pt_src

xp, yp = vec(rd.Vp * x), vec(rd.Vp * y)

#scatter(xp, yp, zcolor=vec(rd.Vp * pt_src), markersize=2, markerstrokewidth=0, legend=false, clims=(-0.25, 0.25))

cache = RHSCache(u, rd)

params = (; rd, md, pt_src, cache)
du = similar(u)
tspan = (0.0, 1.0)
ode = ODEProblem(rhs!, u, tspan, params)


sol = solve(ode, RK4(), saveat=LinRange(tspan..., 10), progress=true, progress_steps=1, verbose=true)

xp = vec(rd.Vp * x)
yp = vec(rd.Vp * y)

frames = [vec(rd.Vp * ui[:, :, 5]) for ui in sol.u]

xmin, xmax = extrema(xp)
ymin, ymax = extrema(yp)

xp = vec(rd.Vp * x)
yp = vec(rd.Vp * y)

@gif for i in eachindex(sol.u)
    scatter(xp, yp, zcolor=vec(rd.Vp * sol.u[i][:, :, 1]),
        markersize=2, markerstrokewidth=0, legend=false)
end fps = 1