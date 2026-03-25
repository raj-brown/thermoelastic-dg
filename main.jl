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
K1D = 64

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
#rhs!(du, u, params, 0.0)
tspan = (0.0, 1.0)
ode = ODEProblem(rhs!, u, tspan, params)

sol = solve(ode, RK4(), saveat=LinRange(tspan..., 50), progress=true)

using Plots

xp = vec(rd.Vp * x)
yp = vec(rd.Vp * y)

frames = [vec(rd.Vp * ui[:, :, 5]) for ui in sol.u]

xmin, xmax = extrema(xp)
ymin, ymax = extrema(yp)

# Auto scale (better than fixed -0.25, 0.25)
cmax = maximum(abs, reduce(vcat, frames))
clims = (-cmax, cmax)

anim = @animate for i in eachindex(frames)
    scatter(
        xp, yp;
        zcolor = frames[i],
        markersize = 3,
        markerstrokewidth = 0,
        marker = :circle,
        legend = false,
        color = :balance,
        clims = clims,
        colorbar = true,
        xlims = (xmin, xmax),
        ylims = (ymin, ymax),
        aspect_ratio = :equal,
        xlabel = "x",
        ylabel = "y",
        title = "Wave field (u₅), step $i"
    )
end

gif(anim, "solution.gif", fps = 5)