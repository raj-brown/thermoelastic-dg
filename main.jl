using StartUpDG
using OrdinaryDiffEq
using OrdinaryDiffEqSSPRK
using LinearAlgebra
using Plots
using Revise
using Dates

push!(LOAD_PATH, "src")
using DG

gr()
ENV["GKSwstype"] = "100"

# ============================================================
# Output folder + autoclean old runs
# ============================================================
function make_output_dir(; root="output", keep_last=5)
    mkpath(root)

    run_dir = joinpath(root, Dates.format(now(), "yyyy-mm-dd_HHMMSS"))
    mkpath(run_dir)

    runs = filter(name -> isdir(joinpath(root, name)), readdir(root))
    sort!(runs)

    old_runs = runs[1:max(0, length(runs) - keep_last)]

    for r in old_runs
        rm(joinpath(root, r); recursive=true, force=true)
    end

    return run_dir
end

outdir = make_output_dir(; root="output", keep_last=5)

println("Running thermoelastic LDG solver")
println("Writing outputs to: ", outdir)

# ============================================================
# Mesh
# ============================================================
N = 4
K1D = 32

rd = RefElemData(Quad(), N)
VXY, EToV = uniform_mesh(Quad(), K1D)
md = MeshData(VXY, EToV, rd)
md = make_periodic(md)

(; x, y) = md

# ============================================================
# State: [σxx, σyy, σxy, vx, vy, T, ψ]
# ============================================================
u = zeros(rd.Np, md.K, 7)

T0 = 10.0
u[:, :, 6] .= T0
u[:, :, 7] .= 0.0

#gaussian(x, y, σ) = exp(-((x)^2 + (y)^2) / (2σ^2))
#u[:, :, 6] .+= 10.0 .* gaussian.(x, y, 0.08)

# ============================================================
# Source
# ============================================================
pt_src = zeros(rd.Np, md.K)

ids = findall(abs.(x[:]) .+ abs.(y[:]) .< 1e-8)
pt_src[ids] .= 1.0
pt_src = (rd.VDM * rd.VDM') * pt_src

# ============================================================
# Cache + params
# ============================================================
cache = RHSCacheTE(u, rd)
params = (; rd, md, pt_src, cache)

# ============================================================
# Material parameters
# ============================================================
ρ = 1.0
μ = 1.0
λ = 1.0
τ = 1.0

c11 = λ + 2μ

dt, h, c_elastic = compute_dt(c11, ρ, τ, md)
CFL = compute_cfl(c_elastic, dt, md, rd)

tspan = (0.0, 1.0)

# ============================================================
# Solve
# ============================================================
prob = ODEProblem(rhs_thermoelastic_cldg!, u, tspan, params)
cb = make_progress_callback(tspan, dt, CFL)

sol = solve(
    prob,
    SSPRK43(),
    dt = dt,
    adaptive = false,
    saveat = LinRange(tspan[1], tspan[2], 40),
    progress = true,
    callback = cb
)

println("Solve complete")

# ============================================================
# Visualization
# ============================================================
println("Writing VTK snapshots...")

xplot = rd.Vp * x
yplot = rd.Vp * y

for (i, t) in enumerate(sol.t)
    uplot = similar(sol.u[i], size(rd.Vp, 1), md.K, 7)

    for fld in 1:7
        uplot[:, :, fld] = rd.Vp * sol.u[i][:, :, fld]
    end

    filename = joinpath(outdir, "vtk_fields_$(lpad(i, 4, '0')).vtu")
    export_quad_subcells_vtu(filename, xplot, yplot, uplot)
end

xp = vec(xplot)
yp = vec(yplot)

println("Writing animations...")

# Velocity animation
anim_v = @animate for i in eachindex(sol.u)
    zi = vec(rd.Vp * sol.u[i][:, :, 5])

    scatter(
        xp, yp,
        zcolor = zi,
        markersize = 2,
        markerstrokewidth = 0,
        legend = false,
        title = "v₂, t=$(round(sol.t[i], sigdigits=3))"
    )
end

mp4(anim_v, joinpath(outdir, "v2.mp4"), fps=10)

# Temperature animation
anim_T = @animate for i in eachindex(sol.u)
    zi = vec(rd.Vp * sol.u[i][:, :, 7])

    scatter(
        xp, yp,
        zcolor = zi,
        markersize = 2,
        markerstrokewidth = 0,
        legend = false,
        clims = (T0 - 5, T0 + 5),
        title = "T, t=$(round(sol.t[i], sigdigits=3))"
    )
end

mp4(anim_T, joinpath(outdir, "T.mp4"), fps = 10)

println("Animations saved")

# ============================================================
# Final field plot
# ============================================================
p = plot(layout = (2, 4), size = (1400, 700))

field_names = ["σxx", "σyy", "σxy", "v₁", "v₂", "T", "ψ"]

for i in 1:7
    zi = vec(rd.Vp * sol.u[end][:, :, i])
    scatter!(
        p[i],
        xp, yp,
        zcolor = zi,
        markersize = 1,
        markerstrokewidth = 0,
        legend = false,
        title = field_names[i]
    )
end

plot!(p[8], framestyle = :none)

savefig(p, joinpath(outdir, "final.png"))

println("Saved final.png")
println("All outputs saved in: ", outdir)
