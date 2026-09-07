# Carcione et al. (2019), Figure 6 homogeneous heat-source experiment.
# DOI: 10.1190/geo2018-0448.1
#
# Mechanical penalty structure follows Shukla, Chan, de Hoop, and Jaiswal
# (JCP 403, 109061, 2020), DOI: 10.1016/j.jcp.2019.109061.
#
# Proposed DG discretization:
#   * elastic stress: weak form (evaluated by the SBP strong-equivalent form),
#   * particle velocity: strong form,
#   * Chan--Shukla scalar penalty flux using total thermoelastic stress,
#   * temperature: weak form,
#   * Cattaneo heat flux: strong form,
#   * BR1 arithmetic-average thermal traces.
#
# The state is
#   [s_xx, s_zz, s_xz, v_x, v_z, theta, q_x, q_z],
# where s is elastic stress and total stress is Sigma = s - b*theta.

ENV["GKSwstype"] = "100"

using Dates
using LinearAlgebra
using Printf
using StartUpDG
using OrdinaryDiffEq
using OrdinaryDiffEqSSPRK
using Plots

include(joinpath(@__DIR__, "src", "DG.jl"))
using .DG

gr()

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
function make_output_dir(; root=joinpath(@__DIR__, "output"), keep_last=5)
    mkpath(root)

    stamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS_sss")
    run_dir = joinpath(root, stamp)
    mkpath(run_dir)

    timestamp_pattern = r"^\d{4}-\d{2}-\d{2}_\d{6}_\d{3}$"
    runs = filter(
        name -> occursin(timestamp_pattern, name) && isdir(joinpath(root, name)),
        readdir(root),
    )
    sort!(runs)

    for name in runs[1:max(0, length(runs) - keep_last)]
        rm(joinpath(root, name); recursive=true, force=true)
    end

    return run_dir
end

parse_env(::Type{T}, name::AbstractString, default) where {T} =
    parse(T, get(ENV, name, string(default)))

function unit_range(field)
    scale = maximum(abs, field)
    if !isfinite(scale) || scale <= eps(eltype(field))
        return zero.(field), zero(scale)
    end
    return field ./ scale, scale
end

# -----------------------------------------------------------------------------
# Carcione Figure 6 experiment settings
# -----------------------------------------------------------------------------
# The paper uses a 231 x 231 grid with dx = dz = 0.1 mm on a 2.3 cm square,
# a centered heat source of central frequency 3.5 MHz, and a snapshot at 3 us.
const DOMAIN_LENGTH = 2.3e-2       # m
const REFERENCE_DX = 1.0e-4        # m
const PUBLISHED_FINAL_TIME = 3.0e-6  # s
const PUBLISHED_DT = 1.0e-8          # s; Carcione's Fourier/splitting solver
const SOURCE_FREQUENCY = 3.5e6       # Hz, used exactly in Carcione's equation (35)
const SOURCE_DELAY = 3 / (2 * SOURCE_FREQUENCY)

final_time = parse_env(Float64, "DG_FINAL_TIME", PUBLISHED_FINAL_TIME)
N = parse_env(Int, "DG_N", 4)
# For N=4, 58 elements give a nominal within-element spacing close to 0.1 mm.
default_K1D = ceil(Int, (DOMAIN_LENGTH / REFERENCE_DX) / N)
K1D = parse_env(Int, "DG_K1D", default_K1D)
CFL = parse_env(Float64, "DG_CFL", 0.20)
penalty_scale = parse_env(Float64, "DG_PENALTY_SCALE", 1.0)
br1_penalty_scale = parse_env(Float64, "DG_BR1_PENALTY_SCALE", 0.0)
source_sigma = parse_env(Float64, "DG_SOURCE_SIGMA", 1.5 * REFERENCE_DX)
source_amplitude = parse_env(Float64, "DG_SOURCE_AMPLITUDE", 1.0)
nsaves = parse_env(Int, "DG_NSAVES", 31)
nvtk = parse_env(Int, "DG_NVTK", 7)

final_time > 0 || error("DG_FINAL_TIME must be positive")
N >= 1 || error("DG_N must be at least 1")
K1D >= 1 || error("DG_K1D must be at least 1")
nsaves >= 2 || error("DG_NSAVES must be at least 2")
nvtk >= 1 || error("DG_NVTK must be at least 1")

outdir = make_output_dir()
println("Carcione Figure 6 thermoelastic DG experiment")
println("Output directory: ", outdir)

# -----------------------------------------------------------------------------
# SBP/DGSEM quadrilateral mesh
# -----------------------------------------------------------------------------
# The SBP reference element makes the weak and strong forms algebraically
# equivalent under the discrete summation-by-parts identity used by the proof.
rd = RefElemData(Quad(), SBP(), N)
VXY_ref, EToV = uniform_mesh(Quad(), K1D)
half_length = DOMAIN_LENGTH / 2
VXY = (half_length .* VXY_ref[1], half_length .* VXY_ref[2])
md = MeshData(VXY, EToV, rd)

# Carcione uses Fourier differentiation for the homogeneous snapshots.  The
# Fourier grid is periodic; we use the corresponding periodic DG connectivity.
md = make_periodic(md)
(; x, y) = md  # y is the paper's z coordinate.

# -----------------------------------------------------------------------------
# Material and characteristic speeds
# -----------------------------------------------------------------------------
material = carcione_reference_material(
    penalty_scale=penalty_scale,
    br1_penalty_scale=br1_penalty_scale,
)
properties = carcione_reference_properties(material)

@printf("Polynomial degree N              = %d\n", N)
@printf("Elements per coordinate         = %d\n", K1D)
@printf("Total quadrilateral elements    = %d\n", md.K)
@printf("Nominal nodal spacing           = %.4e m\n", DOMAIN_LENGTH / (K1D * N))
@printf("Isothermal P speed V_I          = %.2f m/s\n", properties.vI)
@printf("S-wave speed V_S                = %.2f m/s\n", properties.vS)
@printf("Adiabatic P speed V_A           = %.2f m/s\n", properties.vA)
@printf("High-frequency E speed V_Einf   = %.2f m/s\n", properties.vEinf)
@printf("High-frequency T speed V_Tinf   = %.2f m/s\n", properties.vTinf)
@printf("Thermal relaxation time tau     = %.4e s\n", material.tau)
@printf("Carcione reported time step     = %.4e s\n", PUBLISHED_DT)
@printf("Mechanical penalty scale        = %.3f\n", penalty_scale)
@printf("BR1 jump-penalty scale          = %.3f\n", br1_penalty_scale)
effective_propagation_time = max(final_time - SOURCE_DELAY, 0.0)
@printf(
    "Expected E-front radius range   = %.3f--%.3f mm\n",
    1e3 * properties.vA * effective_propagation_time,
    1e3 * properties.vEinf * effective_propagation_time,
)
@printf(
    "High-frequency T-front radius   = %.3f mm\n",
    1e3 * properties.vTinf * effective_propagation_time,
)

# -----------------------------------------------------------------------------
# Initial state and centered regularized heat source
# -----------------------------------------------------------------------------
u0 = zeros(eltype(x), size(x, 1), size(x, 2), NSTATE)

source_shape = normalized_gaussian_source(
    rd,
    md;
    x0=0.0,
    z0=0.0,
    sigma=source_sigma,
)

# Carcione et al. do not report an absolute source amplitude for Figure 6.
# The default is therefore a unit line-source strength (per unit out-of-plane
# thickness).  Because the PDE is linear, changing DG_SOURCE_AMPLITUDE rescales
# all fields without changing wave speeds or arrival times.
source_filter_dt = parse_env(
    Float64,
    "DG_SOURCE_FILTER_DT",
    min(material.tau / 20, inv(200 * SOURCE_FREQUENCY)),
)

source = CarcioneHeatSource(
    source_shape;
    tau=material.tau,
    tmax=final_time,
    amplitude=source_amplitude,
    f0=SOURCE_FREQUENCY,
    t0=SOURCE_DELAY,
    filter_dt=source_filter_dt,
)

cache = RHSCacheTE(u0, rd)
params = (; rd, md, material, cache, source)

# Confirm discrete source normalization.
source_integral = sum(md.wJq .* (rd.Vq * source_shape))
@printf("Discrete spatial source integral = %.16e\n", source_integral)


# Save both the source history printed by Carcione et al. and the equivalent
# first-order supply used by the Cattaneo system.  They satisfy
#     tau * r_t + r = -q_paper.
source_times = collect(range(0.0, final_time; length=max(1001, 20 * nsaves)))
open(joinpath(outdir, "source_history.csv"), "w") do io
    println(io, "time_s,paper_q,equivalent_first_order_supply_r")
    for t in source_times
        @printf(
            io,
            "%.16e,%.16e,%.16e\n",
            t,
            paper_heat_source(source, t),
            first_order_heat_supply(source, t),
        )
    end
end

# -----------------------------------------------------------------------------
# Explicit time step and solve
# -----------------------------------------------------------------------------
dt, dtinfo = estimate_dt(material, rd, md; CFL=CFL, relaxation_fraction=0.5)
reported_cfl = compute_cfl(dt, material, rd, md)

@printf("Element width estimate h         = %.4e m\n", dtinfo.h)
@printf("Coupled maximum speed            = %.2f m/s\n", dtinfo.cmax)
@printf("Wave-limited dt                  = %.4e s\n", dtinfo.dt_wave)
@printf("Relaxation-limited dt            = %.4e s\n", dtinfo.dt_relax)
@printf("Selected dt                      = %.4e s (%s limited)\n", dt, string(dtinfo.limiting))
@printf("Reported DG CFL                  = %.4f\n", reported_cfl)
@printf("Estimated number of steps        = %d\n", ceil(Int, final_time / dt))

save_times = collect(range(0.0, final_time; length=nsaves))
prob = ODEProblem(rhs_thermoelastic_br1!, u0, (0.0, final_time), params)

println("Starting SSPRK54 solve...")
sol = solve(
    prob,
    SSPRK54();
    dt=dt,
    adaptive=false,
    saveat=save_times,
    save_start=true,
    save_end=true,
    save_everystep=false,
    progress=true,
)
println("Solve complete with ", length(sol.t), " stored snapshots")

# -----------------------------------------------------------------------------
# VTK snapshots
# -----------------------------------------------------------------------------
xplot = rd.Vp * x
zplot = rd.Vp * y
vtk_ids = unique(round.(Int, range(1, length(sol.u); length=min(nvtk, length(sol.u)))))

println("Writing VTK snapshots...")
for (counter, index) in enumerate(vtk_ids)
    uplot = interpolate_state(rd, sol.u[index])
    filename = joinpath(outdir, @sprintf("carcione_fig6_%04d.vtu", counter))
    export_quad_subcells_vtu(filename, xplot, zplot, uplot, material)
end

# -----------------------------------------------------------------------------
# Figure 6-style final snapshot: vertical velocity and temperature increment
# -----------------------------------------------------------------------------
ufinal = interpolate_state(rd, sol.u[end])
vz_final = vec(@view ufinal[:, :, IVZ])
theta_final = vec(@view ufinal[:, :, ITHETA])
vz_normalized, vz_scale = unit_range(vz_final)
theta_normalized, theta_scale = unit_range(theta_final)

x_mm = 1e3 .* vec(xplot)
z_mm = 1e3 .* vec(zplot)

p_vz = scatter(
    x_mm,
    z_mm;
    marker_z=vz_normalized,
    markersize=1.7,
    markerstrokewidth=0,
    color=:grays,
    clims=(-1, 1),
    legend=false,
    colorbar=true,
    aspect_ratio=:equal,
    yflip=true,
    xlabel="x (mm)",
    ylabel="z (mm)",
    title=@sprintf("v_z / max|v_z| at %.3f us\nmax|v_z| = %.3e m/s", 1e6 * sol.t[end], vz_scale),
)

p_theta = scatter(
    x_mm,
    z_mm;
    marker_z=theta_normalized,
    markersize=1.7,
    markerstrokewidth=0,
    color=:grays,
    clims=(-1, 1),
    legend=false,
    colorbar=true,
    aspect_ratio=:equal,
    yflip=true,
    xlabel="x (mm)",
    ylabel="z (mm)",
    title=@sprintf("theta / max|theta| at %.3f us\nmax|theta| = %.3e K", 1e6 * sol.t[end], theta_scale),
)

snapshot_plot = plot(p_vz, p_theta; layout=(1, 2), size=(1300, 560))
savefig(snapshot_plot, joinpath(outdir, "carcione_figure6_DG.png"))

# -----------------------------------------------------------------------------
# Discrete energy history
# -----------------------------------------------------------------------------
energy = [energy_components(U, rd, md, material) for U in sol.u]
energy_table = hcat(
    sol.t,
    getproperty.(energy, :elastic),
    getproperty.(energy, :kinetic),
    getproperty.(energy, :temperature),
    getproperty.(energy, :heat_flux),
    getproperty.(energy, :total),
    getproperty.(energy, :relaxation),
)

open(joinpath(outdir, "energy.csv"), "w") do io
    println(io, "time_s,elastic_J,kinetic_J,temperature_J,heat_flux_J,total_J,relaxation_W")
    for row in eachrow(energy_table)
        println(io, join((@sprintf("%.16e", value) for value in row), ","))
    end
end

energy_scale = maximum(getproperty.(energy, :total))
energy_scale = energy_scale > 0 ? energy_scale : 1.0
p_energy = plot(
    1e6 .* sol.t,
    getproperty.(energy, :total) ./ energy_scale;
    linewidth=2,
    label="total",
    xlabel="time (us)",
    ylabel="energy / max energy",
    title="Discrete thermoelastic energy",
)
plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :elastic) ./ energy_scale; label="elastic")
plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :kinetic) ./ energy_scale; label="kinetic")
plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :temperature) ./ energy_scale; label="temperature")
plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :heat_flux) ./ energy_scale; label="heat flux")
savefig(p_energy, joinpath(outdir, "energy_history.png"))

# -----------------------------------------------------------------------------
# Instantaneous semi-discrete energy-balance audit
# -----------------------------------------------------------------------------
du_audit = similar(u0)
balance = map(eachindex(sol.t)) do i
    t = sol.t[i]
    U = sol.u[i]
    rhs_thermoelastic_br1!(du_audit, U, params, t)
    energy_balance(
        U,
        du_audit,
        rd,
        md,
        material;
        source_supply=first_order_heat_supply(source, t),
        source_spatial=source.spatial,
    )
end

balance_scale = map(balance) do item
    max(
        abs(item.denergy_dt) + item.relaxation + item.mechanical_face +
        item.thermal_face + abs(item.source_power),
        eps(Float64),
    )
end
relative_balance_residual = abs.(getproperty.(balance, :residual)) ./ balance_scale

open(joinpath(outdir, "energy_balance.csv"), "w") do io
    println(
        io,
        "time_s,dE_dt_W,D_rel_W,D_mech_W,D_BR1_W,source_power_W,residual_W,relative_residual",
    )
    for i in eachindex(sol.t)
        item = balance[i]
        @printf(
            io,
            "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
            sol.t[i],
            item.denergy_dt,
            item.relaxation,
            item.mechanical_face,
            item.thermal_face,
            item.source_power,
            item.residual,
            relative_balance_residual[i],
        )
    end
end

@printf(
    "Maximum relative semi-discrete energy residual = %.3e\n",
    maximum(relative_balance_residual),
)

# -----------------------------------------------------------------------------
# Run metadata
# -----------------------------------------------------------------------------
open(joinpath(outdir, "run_summary.txt"), "w") do io
    println(io, "Carcione et al. (2019), Figure 6-inspired heat-source experiment")
    println(io, "DOI: 10.1190/GEO2018-0448.1")
    println(io)
    println(io, "DG scheme: weak elastic stress / strong velocity")
    println(io, "Mechanical flux: Chan--Shukla scalar penalty flux")
    println(io, "Thermal scheme: weak temperature / strong Cattaneo heat flux, BR1 traces")
    println(io)
    @printf(io, "N = %d\n", N)
    @printf(io, "K1D = %d\n", K1D)
    @printf(io, "domain length = %.16e m\n", DOMAIN_LENGTH)
    @printf(io, "final time = %.16e s\n", final_time)
    @printf(io, "Carcione reported dt = %.16e s\n", PUBLISHED_DT)
    @printf(io, "DG selected dt = %.16e s\n", dt)
    @printf(io, "CFL = %.16e\n", reported_cfl)
    @printf(io, "source f0 = %.16e Hz\n", SOURCE_FREQUENCY)
    @printf(io, "source t0 = %.16e s\n", SOURCE_DELAY)
    @printf(io, "source sigma = %.16e m\n", source_sigma)
    @printf(io, "source amplitude = %.16e\n", source_amplitude)
    @printf(io, "source spatial integral = %.16e\n", source_integral)
    @printf(io, "source filter dt = %.16e s\n", source.filter_dt)
    @printf(io, "effective propagation time = %.16e s\n", effective_propagation_time)
    @printf(io, "adiabatic E-front radius = %.16e m\n", properties.vA * effective_propagation_time)
    @printf(io, "high-frequency E-front radius = %.16e m\n", properties.vEinf * effective_propagation_time)
    @printf(io, "high-frequency T-front radius = %.16e m\n", properties.vTinf * effective_propagation_time)
    @printf(io, "maximum relative energy-balance residual = %.16e\n", maximum(relative_balance_residual))
    println(io)
    println(io, "The paper specifies the source waveform and location but not its absolute amplitude.")
    println(io, "The printed second-order heat source q_paper is represented exactly in the")
    println(io, "first-order Cattaneo system by the causal supply r satisfying")
    println(io, "tau*r_t + r = -q_paper, r(0)=0.")
    println(io, "The spatial delta is regularized by a normalized Gaussian; linearity makes")
    println(io, "the unspecified source amplitude a simple scaling parameter.")
end

println("Saved Figure 6-style image: ", joinpath(outdir, "carcione_figure6_DG.png"))
println("Saved energy history: ", joinpath(outdir, "energy.csv"))
println("Saved energy balance: ", joinpath(outdir, "energy_balance.csv"))
println("Saved source history: ", joinpath(outdir, "source_history.csv"))
println("All outputs saved in: ", outdir)
