# Standalone driver for Carcione et al. (2019), Figure 6.
#
# Experiment: centred heat source in a 2.3 cm x 2.3 cm periodic homogeneous
# medium, thermal conductivity 10.5 W/(m K), source parameter 3.5 MHz, and a
# snapshot at 3 microseconds.
#
# This file does not call or modify main.jl.  It requires the revised
# eight-field src/DG.jl used by the corrected Figure 7 driver:
#
#   [s_xx, s_zz, s_xz, v_x, v_z, theta, q_x, q_z].
#
# The default pulse mode is :printed, matching Carcione's equation (35)
# literally. Set FIG6_WAVELET_MODE=cyclic_hz only for the alternative
# interpretation in which f0 is treated as cycles per second.

module CarcioneFigure6Driver

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using Dates
using LinearAlgebra
using ProgressMeter
using Printf
using StartUpDG
using OrdinaryDiffEq
using OrdinaryDiffEqSSPRK
using Plots


"""
    make_fixed_step_progress_callback(tspan, dt;
        cfl=NaN,
        description="Solving",
        max_updates=200,
    )

Create a throttled `DiscreteCallback` for a fixed-step ODE solve. The
callback advances the meter only after time has moved forward and updates the
terminal at most `max_updates` times. It does not create extra saved states.
"""
function make_fixed_step_progress_callback(
    tspan::Tuple{<:Real,<:Real},
    dt::Real;
    cfl::Real=NaN,
    description::AbstractString="Solving",
    max_updates::Integer=200,
)
    t_initial = Float64(tspan[1])
    t_final = Float64(tspan[2])
    dt_value = Float64(dt)

    dt_value > 0 || throw(ArgumentError("dt must be positive"))
    t_final > t_initial ||
        throw(ArgumentError("tspan must satisfy t_final > t_initial"))
    max_updates > 0 ||
        throw(ArgumentError("max_updates must be positive"))

    estimated_steps = max(
        1,
        ceil(Int, (t_final - t_initial) / dt_value),
    )
    update_interval = max(1, cld(estimated_steps, Int(max_updates)))

    progress = Progress(
        estimated_steps;
        desc=String(description),
        showspeed=true,
        dt=0.25,
    )

    start_wall_time = time()
    step_counter = Ref(0)
    last_displayed_step = Ref(0)
    last_callback_time = Ref(t_initial)
    progress_finished = Ref(false)

    # The callback is active only once time has advanced. This prevents an
    # initialization check, or a repeated check at the same time, from
    # incrementing the progress meter twice.
    condition(u, t, integrator) = t > last_callback_time[]

    function affect!(integrator)
        current_time = Float64(integrator.t)
        last_callback_time[] = current_time
        step_counter[] += 1

        displayed_step = min(step_counter[], estimated_steps)
        fraction = clamp(
            (current_time - t_initial) / (t_final - t_initial),
            0.0,
            1.0,
        )

        elapsed = time() - start_wall_time
        eta = fraction > 0 ? elapsed * (1 - fraction) / fraction : Inf

        final_tolerance =
            100 * eps(Float64) * max(1.0, abs(t_final))
        is_final = current_time >= t_final - final_tolerance

        should_display =
            displayed_step - last_displayed_step[] >= update_interval ||
            is_final

        if should_display
            values = [
                (:time, @sprintf("%.6e / %.6e", current_time, t_final)),
                (:step, @sprintf("%d / %d", displayed_step, estimated_steps)),
                (:dt, @sprintf("%.3e", Float64(integrator.dt))),
                (:CFL, isfinite(cfl) ? @sprintf("%.3e", cfl) : "not supplied"),
                (:elapsed, @sprintf("%.1f s", elapsed)),
                (:ETA, isfinite(eta) ? @sprintf("%.1f s", eta) : "estimating"),
            ]

            ProgressMeter.update!(
                progress,
                is_final ? estimated_steps : displayed_step;
                showvalues=values,
            )
            last_displayed_step[] = displayed_step
        end

        if is_final && !progress_finished[]
            ProgressMeter.finish!(progress)
            progress_finished[] = true
        end

        return nothing
    end

    return DiscreteCallback(
        condition,
        affect!;
        save_positions=(false, false),
    )
end


include(joinpath(@__DIR__, "src", "DG.jl"))
using .DG

gr()

const DOMAIN_LENGTH = 2.3e-2
const REFERENCE_DX = 1.0e-4
const PUBLISHED_FINAL_TIME = 3.0e-6
const PUBLISHED_DT = 1.0e-8
const SOURCE_FREQUENCY = 3.5e6
const SOURCE_DELAY = 3 / (2 * SOURCE_FREQUENCY)
const FIGURE6_CONDUCTIVITY = 10.5

parse_env(::Type{T}, name::AbstractString, default) where {T} =
    parse(T, get(ENV, name, string(default)))

function parse_env_alias(
    ::Type{T},
    primary::AbstractString,
    legacy::AbstractString,
    default,
) where {T}
    raw = haskey(ENV, primary) ? ENV[primary] : get(ENV, legacy, string(default))
    return parse(T, raw)
end

function parse_bool_value(name::AbstractString, raw)
    value = lowercase(strip(string(raw)))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$name must be true/false, yes/no, on/off, or 1/0"))
end

function parse_bool_env(name::AbstractString, default::Bool)
    return parse_bool_value(name, get(ENV, name, string(default)))
end

function parse_bool_env_alias(
    primary::AbstractString,
    legacy::AbstractString,
    default::Bool,
)
    raw = haskey(ENV, primary) ? ENV[primary] : get(ENV, legacy, string(default))
    return parse_bool_value(primary, raw)
end

function string_env_alias(
    primary::AbstractString,
    legacy::AbstractString,
    default::AbstractString,
)
    return haskey(ENV, primary) ? ENV[primary] : get(ENV, legacy, default)
end

function make_output_dir(; root=joinpath(@__DIR__, "output_figure6"), keep_last=8)
    mkpath(root)
    stamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS_sss")
    run_dir = joinpath(root, stamp)
    mkpath(run_dir)

    pattern = r"^\d{4}-\d{2}-\d{2}_\d{6}_\d{3}$"
    runs = filter(
        name -> occursin(pattern, name) && isdir(joinpath(root, name)),
        readdir(root),
    )
    sort!(runs)
    for name in runs[1:max(0, length(runs) - keep_last)]
        rm(joinpath(root, name); recursive=true, force=true)
    end
    return run_dir
end

function normalized_field(field)
    scale = maximum(abs, field)
    if !isfinite(scale) || scale <= eps(Float64)
        return zero.(field), zero(scale)
    end
    return field ./ scale, scale
end

function volume_values(Vq, field)
    if Vq isa UniformScaling
        return Vq.λ .* field
    end
    return Vq * field
end


"""
    broken_physical_derivatives(field, rd, md)

Return the elementwise physical derivatives `(dfield_dx, dfield_dz)` at the
solution nodes. The derivatives are the broken DG derivatives, which are the
appropriate quantities for diagnosing compressional and rotational content
inside the elements. Face-jump diagnostics are reported separately by the
energy audit.
"""
function broken_physical_derivatives(field, rd, md)
    field_r = rd.Dr * field
    field_s = rd.Ds * field

    field_x = @. (md.rxJ * field_r + md.sxJ * field_s) / md.J
    field_z = @. (md.ryJ * field_r + md.syJ * field_s) / md.J
    return field_x, field_z
end

"""
    longitudinal_shear_diagnostics(U, rd, md, material; rmin=0, rmax=Inf)

Distinguish longitudinal E/T content from transverse S-wave content.

For a centred heat source in a homogeneous isotropic medium, both the E-wave
and T-wave are longitudinal. Their velocity is radial, while an S-wave is
transverse. We therefore report

* the tangential-to-radial velocity L2 ratio,
* the tangential kinetic-energy fraction, and
* the vorticity-to-dilatation L2 ratio.

The ratios are evaluated on the annulus `rmin <= r <= rmax`. Exact heat-source
solutions have zero tangential velocity and zero vorticity; nonzero values are
numerical contamination and should decrease under refinement.
"""
function longitudinal_shear_diagnostics(
    U,
    rd,
    md,
    material;
    rmin::Real=0.0,
    rmax::Real=Inf,
)
    rmin >= 0 || throw(ArgumentError("rmin must be nonnegative"))
    rmax > rmin || throw(ArgumentError("rmax must exceed rmin"))

    @views begin
        vx = U[:, :, IVX]
        vz = U[:, :, IVZ]
    end

    xq = volume_values(rd.Vq, md.x)
    zq = volume_values(rd.Vq, md.y)
    vxq = volume_values(rd.Vq, vx)
    vzq = volume_values(rd.Vq, vz)

    radius = @. sqrt(xq^2 + zq^2)
    inv_radius = similar(radius)
    radius_tol = sqrt(eps(Float64)) * max(1.0, maximum(radius))
    @inbounds for i in eachindex(radius)
        inv_radius[i] = radius[i] > radius_tol ? inv(radius[i]) : 0.0
    end

    radial_velocity = @. (xq * vxq + zq * vzq) * inv_radius
    tangential_velocity = @. (-zq * vxq + xq * vzq) * inv_radius

    vx_x, vx_z = broken_physical_derivatives(vx, rd, md)
    vz_x, vz_z = broken_physical_derivatives(vz, rd, md)
    dilatation_q = volume_values(rd.Vq, @. vx_x + vz_z)
    vorticity_q = volume_values(rd.Vq, @. vz_x - vx_z)

    mask = @. (radius >= rmin) & (radius <= rmax)
    weights = md.wJq

    radial_norm_sq = 0.0
    tangential_norm_sq = 0.0
    dilatation_norm_sq = 0.0
    vorticity_norm_sq = 0.0

    @inbounds for i in eachindex(weights)
        if mask[i]
            wi = weights[i]
            radial_norm_sq += wi * radial_velocity[i]^2
            tangential_norm_sq += wi * tangential_velocity[i]^2
            dilatation_norm_sq += wi * dilatation_q[i]^2
            vorticity_norm_sq += wi * vorticity_q[i]^2
        end
    end

    radial_norm = sqrt(max(radial_norm_sq, 0.0))
    tangential_norm = sqrt(max(tangential_norm_sq, 0.0))
    dilatation_norm = sqrt(max(dilatation_norm_sq, 0.0))
    vorticity_norm = sqrt(max(vorticity_norm_sq, 0.0))

    velocity_scale = max(radial_norm, eps(Float64))
    derivative_scale = max(dilatation_norm, eps(Float64))
    kinetic_total = radial_norm_sq + tangential_norm_sq

    return (;
        rmin=Float64(rmin),
        rmax=Float64(rmax),
        radial_velocity_norm=radial_norm,
        tangential_velocity_norm=tangential_norm,
        tangential_to_radial=tangential_norm / velocity_scale,
        tangential_kinetic_fraction=(
            kinetic_total > eps(Float64) ? tangential_norm_sq / kinetic_total : 0.0
        ),
        dilatation_norm=dilatation_norm,
        vorticity_norm=vorticity_norm,
        vorticity_to_dilatation=vorticity_norm / derivative_scale,
    )
end

function write_mode_diagnostics(path, global_diag, candidate_diag, candidate_radius)
    open(path, "w") do io
        println(io, "Longitudinal-versus-shear diagnostic")
        println(io)
        println(io, "For a centered heat source in a homogeneous isotropic medium,")
        println(io, "E and T modes are longitudinal. S-wave content is measured by")
        println(io, "tangential velocity and out-of-plane vorticity.")
        println(io)
        @printf(io, "candidate_S_or_T_radius_m = %.16e\n", candidate_radius)
        println(io)
        println(io, "global_annulus")
        @printf(io, "rmin_m = %.16e\n", global_diag.rmin)
        @printf(io, "rmax_m = %.16e\n", global_diag.rmax)
        @printf(io, "tangential_to_radial_L2 = %.16e\n", global_diag.tangential_to_radial)
        @printf(io, "tangential_kinetic_fraction = %.16e\n", global_diag.tangential_kinetic_fraction)
        @printf(io, "vorticity_to_dilatation_L2 = %.16e\n", global_diag.vorticity_to_dilatation)
        println(io)
        println(io, "candidate_ring_annulus")
        @printf(io, "rmin_m = %.16e\n", candidate_diag.rmin)
        @printf(io, "rmax_m = %.16e\n", candidate_diag.rmax)
        @printf(io, "tangential_to_radial_L2 = %.16e\n", candidate_diag.tangential_to_radial)
        @printf(io, "tangential_kinetic_fraction = %.16e\n", candidate_diag.tangential_kinetic_fraction)
        @printf(io, "vorticity_to_dilatation_L2 = %.16e\n", candidate_diag.vorticity_to_dilatation)
    end
end

function carrier_angular_frequency(f0::Real, wavelet_mode::Symbol)
    wavelet_mode === :cyclic_hz && return 2pi * f0
    wavelet_mode === :printed && return f0
    throw(ArgumentError("wavelet_mode must be :printed or :cyclic_hz"))
end

function positive_real_sqrt(value::Complex)
    root = sqrt(value)
    return real(root) < 0 ? -root : root
end

"""
Return the finite-frequency E- and T-mode phase velocities and attenuation
coefficients from Carcione's coupled compressional-wave dispersion relation.
"""
function carcione_mode_diagnostics(
    material::ThermoelasticMaterial,
    angular_frequency::Real,
)
    angular_frequency > 0 || throw(ArgumentError("angular_frequency must be positive"))

    properties = carcione_reference_properties(material)
    diffusivity = material.K[1, 1] / material.heat_capacity
    kernel = im * angular_frequency * diffusivity /
             (1 + im * angular_frequency * material.tau)

    discriminant =
        (properties.vA^2 + kernel)^2 -
        4 * kernel * properties.vI^2

    velocity_squared = (
        0.5 * (properties.vA^2 + kernel + sqrt(discriminant)),
        0.5 * (properties.vA^2 + kernel - sqrt(discriminant)),
    )

    modes = map(velocity_squared) do value
        complex_velocity = positive_real_sqrt(complex(value))
        slowness = inv(complex_velocity)
        phase_velocity = inv(real(slowness))
        attenuation = max(0.0, -angular_frequency * imag(slowness))
        attenuation_length = attenuation > 0 ? inv(attenuation) : Inf
        return (;
            complex_velocity,
            phase_velocity,
            attenuation,
            attenuation_length,
        )
    end

    ordered = sort(collect(modes); by=mode -> mode.phase_velocity)
    return (; thermal=ordered[1], elastic=ordered[end])
end

function figure6_material(
    conductivity::Real;
    penalty_scale::Real=1.0,
    br1_penalty_scale::Real=0.0,
)
    rho = 2650.0
    vI = 2457.0
    vS = 1505.0
    heat_capacity = 117.0
    Tref = 300.0
    beta = 7.92e4

    mu = rho * vS^2
    lambda = rho * vI^2 - 2mu
    tau = conductivity / (heat_capacity * vI^2)

    C = [
        lambda + 2mu  lambda          0.0
        lambda        lambda + 2mu    0.0
        0.0           0.0             mu
    ]
    b = [beta, beta, 0.0]
    K = [conductivity 0.0; 0.0 conductivity]

    return make_material(
        rho=rho,
        C=C,
        b=b,
        heat_capacity=heat_capacity,
        Tref=Tref,
        K=K,
        tau=tau,
        penalty_scale=penalty_scale,
        br1_penalty_scale=br1_penalty_scale,
    )
end

function write_source_history(path, source, final_time, nsaves)
    times = collect(range(0.0, final_time; length=max(1001, 20 * nsaves)))
    open(path, "w") do io
        println(io, "time_s,paper_q,equivalent_first_order_supply_r")
        for t in times
            @printf(
                io,
                "%.16e,%.16e,%.16e\n",
                t,
                paper_heat_source(source, t),
                first_order_heat_supply(source, t),
            )
        end
    end
end

function write_energy_history(path, sol, rd, md, material)
    energy = [energy_components(U, rd, md, material) for U in sol.u]
    open(path, "w") do io
        println(io, "time_s,elastic_J,kinetic_J,temperature_J,heat_flux_J,total_J,relaxation_W")
        for i in eachindex(sol.t)
            item = energy[i]
            @printf(
                io,
                "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                sol.t[i],
                item.elastic,
                item.kinetic,
                item.temperature,
                item.heat_flux,
                item.total,
                item.relaxation,
            )
        end
    end
    return energy
end

function write_energy_balance(path, sol, params, rd, md, material, source)
    du = similar(sol.u[1])
    balances = map(eachindex(sol.t)) do i
        t = sol.t[i]
        U = sol.u[i]
        rhs_thermoelastic_br1!(du, U, params, t)
        energy_balance(
            U,
            du,
            rd,
            md,
            material;
            source_supply=first_order_heat_supply(source, t),
            source_spatial=source.spatial,
        )
    end

    relative = similar(sol.t)
    for i in eachindex(balances)
        item = balances[i]
        scale = max(
            abs(item.denergy_dt) + item.relaxation + item.mechanical_face +
            item.thermal_face + abs(item.source_power),
            eps(Float64),
        )
        relative[i] = abs(item.residual) / scale
    end

    open(path, "w") do io
        println(io, "time_s,dE_dt_W,D_rel_W,D_mech_W,D_BR1_W,source_power_W,residual_W,relative_residual")
        for i in eachindex(sol.t)
            item = balances[i]
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
                relative[i],
            )
        end
    end

    return maximum(relative)
end

function radial_rms_profile(x, z, field; nbins=240)
    (size(x) == size(z) && size(z) == size(field)) ||
        throw(DimensionMismatch("x, z, and field must have identical sizes"))

    radius = sqrt.(x .^ 2 .+ z .^ 2)

    # Restrict the profile to the largest full circle contained in the square.
    # Bins beyond this radius would contain only corner samples and would not
    # represent a genuine angular average.
    rmax = min(maximum(abs, x), maximum(abs, z))
    edges = collect(range(0.0, rmax; length=nbins + 1))
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])
    sums = zeros(Float64, nbins)
    counts = zeros(Int, nbins)

    for i in eachindex(radius)
        radius[i] <= rmax || continue
        bin = clamp(searchsortedlast(edges, radius[i]), 1, nbins)
        sums[bin] += abs2(field[i])
        counts[bin] += 1
    end

    profile = zeros(Float64, nbins)
    for i in eachindex(profile)
        profile[i] = counts[i] == 0 ? NaN : sqrt(sums[i] / counts[i])
    end
    return centers, profile
end

function main()
    isdefined(DG, :rhs_thermoelastic_br1!) || error(
        "src/DG.jl does not export rhs_thermoelastic_br1!. " *
        "This driver requires the revised eight-field solver.",
    )
    isdefined(DG, :carcione_pulse) || error(
        "src/DG.jl does not provide the corrected Carcione pulse API. " *
        "Use the same revised DG.jl used by the fixed Figure 7 driver.",
    )
    DG.NSTATE == 8 || error("Figure 6 driver expects the eight-field 2D state")

    final_time = parse_env_alias(
        Float64,
        "FIG6_FINAL_TIME",
        "DG_FINAL_TIME",
        PUBLISHED_FINAL_TIME,
    )
    N = parse_env_alias(Int, "FIG6_N", "DG_N", 4)

    # Match the paper's nominal 0.1 mm sampling.  Make the element count odd
    # so that the source lies at the centre of a single element rather than at
    # the intersection of four elements.
    nominal_elements = ceil(Int, (DOMAIN_LENGTH / REFERENCE_DX) / N)
    default_K1D = isodd(nominal_elements) ? nominal_elements : nominal_elements + 1
    K1D = parse_env_alias(Int, "FIG6_K1D", "DG_K1D", default_K1D)

    cfl_target = parse_env_alias(Float64, "FIG6_CFL", "DG_CFL", 0.20)
    penalty_scale = parse_env_alias(
        Float64,
        "FIG6_PENALTY_SCALE",
        "DG_PENALTY_SCALE",
        1.0,
    )
    br1_penalty_scale = parse_env_alias(
        Float64,
        "FIG6_BR1_PENALTY_SCALE",
        "DG_BR1_PENALTY_SCALE",
        0.0,
    )
    conductivity = parse_env_alias(
        Float64,
        "FIG6_CONDUCTIVITY",
        "DG_CONDUCTIVITY",
        FIGURE6_CONDUCTIVITY,
    )
    source_sigma = parse_env_alias(
        Float64,
        "FIG6_SOURCE_SIGMA",
        "DG_SOURCE_SIGMA",
        1.5 * REFERENCE_DX,
    )
    source_amplitude = parse_env_alias(
        Float64,
        "FIG6_SOURCE_AMPLITUDE",
        "DG_SOURCE_AMPLITUDE",
        1.0,
    )
    nsaves = parse_env_alias(Int, "FIG6_NSAVES", "DG_NSAVES", 41)
    nvtk = parse_env_alias(Int, "FIG6_NVTK", "DG_NVTK", 7)
    write_vtk = parse_bool_env_alias("FIG6_WRITE_VTK", "DG_WRITE_VTK", true)
    display_clip = parse_env_alias(
        Float64,
        "FIG6_DISPLAY_CLIP",
        "DG_DISPLAY_CLIP",
        1.0,
    )
    progress_updates = parse_env(Int, "FIG6_PROGRESS_UPDATES", 200)
    radial_bins = parse_env(Int, "FIG6_RADIAL_BINS", 240)
    relaxation_fraction = parse_env(Float64, "FIG6_RELAXATION_FRACTION", 0.5)

    wavelet_name = lowercase(strip(string_env_alias(
        "FIG6_WAVELET_MODE",
        "DG_WAVELET_MODE",
        "printed",
    )))
    wavelet_mode = if wavelet_name == "printed"
        :printed
    elseif wavelet_name in ("cyclic_hz", "hz", "frequency")
        :cyclic_hz
    else
        error("FIG6_WAVELET_MODE must be printed or cyclic_hz")
    end

    final_time > 0 || error("FIG6_FINAL_TIME must be positive")
    N >= 1 || error("FIG6_N must be at least 1")
    K1D >= 1 || error("FIG6_K1D must be at least 1")
    cfl_target > 0 || error("FIG6_CFL must be positive")
    penalty_scale >= 0 || error("FIG6_PENALTY_SCALE must be nonnegative")
    br1_penalty_scale >= 0 || error("FIG6_BR1_PENALTY_SCALE must be nonnegative")
    conductivity > 0 || error("FIG6_CONDUCTIVITY must be positive")
    source_sigma > 0 || error("FIG6_SOURCE_SIGMA must be positive")
    nsaves >= 2 || error("FIG6_NSAVES must be at least 2")
    nvtk >= 1 || error("FIG6_NVTK must be at least 1")
    radial_bins >= 16 || error("FIG6_RADIAL_BINS must be at least 16")
    progress_updates >= 1 || error("FIG6_PROGRESS_UPDATES must be positive")
    relaxation_fraction > 0 || error("FIG6_RELAXATION_FRACTION must be positive")
    0 < display_clip <= 1 || error("FIG6_DISPLAY_CLIP must lie in (0,1]")

    outdir = make_output_dir()
    println("Carcione Figure 6 standalone DG experiment")
    println("Output directory: ", outdir)

    # ---------------------------------------------------------------------
    # Periodic SBP/DGSEM mesh
    # ---------------------------------------------------------------------
    rd = RefElemData(Quad(), SBP(), N)
    VXY_ref, EToV = uniform_mesh(Quad(), K1D)
    half_length = DOMAIN_LENGTH / 2
    VXY = (half_length .* VXY_ref[1], half_length .* VXY_ref[2])
    md = MeshData(VXY, EToV, rd)
    md = make_periodic(md)
    (; x, y) = md

    # ---------------------------------------------------------------------
    # Low-conductivity Figure 6 material and dispersion diagnostics
    # ---------------------------------------------------------------------
    material = figure6_material(
        conductivity;
        penalty_scale=penalty_scale,
        br1_penalty_scale=br1_penalty_scale,
    )
    properties = carcione_reference_properties(material)

    omega_carrier = carrier_angular_frequency(SOURCE_FREQUENCY, wavelet_mode)
    carrier_frequency_hz = omega_carrier / (2pi)
    modes = carcione_mode_diagnostics(material, omega_carrier)

    effective_time = max(final_time - SOURCE_DELAY, 0.0)
    expected_e_radius = modes.elastic.phase_velocity * effective_time
    formal_t_radius = modes.thermal.phase_velocity * effective_time
    thermal_decay_exponent = modes.thermal.attenuation * formal_t_radius

    @printf("N                                  = %d\n", N)
    @printf("elements per coordinate            = %d\n", K1D)
    @printf("quadrilateral elements             = %d\n", md.K)
    @printf("nominal nodal spacing              = %.4e m\n", DOMAIN_LENGTH / (K1D * N))
    @printf("Figure 6 conductivity gamma        = %.6e W/(m K)\n", conductivity)
    @printf("thermal relaxation time tau        = %.6e s\n", material.tau)
    @printf("V_I                                = %.2f m/s\n", properties.vI)
    @printf("V_S                                = %.2f m/s\n", properties.vS)
    @printf("V_A                                = %.2f m/s\n", properties.vA)
    @printf("V_Einf                             = %.2f m/s\n", properties.vEinf)
    @printf("V_Tinf                             = %.2f m/s\n", properties.vTinf)
    @printf("source wavelet convention          = %s\n", String(wavelet_mode))
    @printf("actual carrier frequency           = %.6e Hz\n", carrier_frequency_hz)
    @printf("finite-frequency E phase speed     = %.2f m/s\n", modes.elastic.phase_velocity)
    @printf("finite-frequency E attenuation     = %.4e 1/m\n", modes.elastic.attenuation)
    @printf("finite-frequency T phase speed     = %.2f m/s\n", modes.thermal.phase_velocity)
    @printf("finite-frequency T attenuation     = %.4e 1/m\n", modes.thermal.attenuation)
    @printf("T-wave attenuation length          = %.4e m\n", modes.thermal.attenuation_length)
    @printf("source-peak E radius prediction    = %.3f mm\n", 1e3 * expected_e_radius)
    @printf("formal source-peak T radius        = %.3f mm\n", 1e3 * formal_t_radius)
    @printf("T attenuation exponent at radius   = %.2f\n", thermal_decay_exponent)
    println("The Figure 6 T mode is expected to be diffusive and concentrated near the source.")

    # ---------------------------------------------------------------------
    # Initial state and centred regularised heat source
    # ---------------------------------------------------------------------
    u0 = zeros(eltype(x), size(x, 1), size(x, 2), NSTATE)
    source_shape = normalized_gaussian_source(
        rd,
        md;
        x0=0.0,
        z0=0.0,
        sigma=source_sigma,
    )

    source_filter_dt = parse_env_alias(
        Float64,
        "FIG6_SOURCE_FILTER_DT",
        "DG_SOURCE_FILTER_DT",
        min(material.tau / 20, inv(200 * SOURCE_FREQUENCY)),
    )
    source_filter_dt > 0 || error("FIG6_SOURCE_FILTER_DT must be positive")

    source = CarcioneHeatSource(
        source_shape;
        tau=material.tau,
        tmax=final_time,
        amplitude=source_amplitude,
        f0=SOURCE_FREQUENCY,
        t0=SOURCE_DELAY,
        filter_dt=source_filter_dt,
        wavelet_mode=wavelet_mode,
    )

    cache = RHSCacheTE(u0, rd)
    params = (; rd, md, material, cache, source)

    source_q = volume_values(rd.Vq, source_shape)
    source_integral = sum(md.wJq .* source_q)
    @printf("discrete source integral           = %.16e\n", source_integral)
    isapprox(source_integral, 1.0; atol=1e-11, rtol=1e-11) ||
        error("normalised source integral is $source_integral rather than one")

    write_source_history(
        joinpath(outdir, "source_history.csv"),
        source,
        final_time,
        nsaves,
    )

    # ---------------------------------------------------------------------
    # Explicit time step, progress callback, and solve
    # ---------------------------------------------------------------------
    dt, dtinfo = estimate_dt(
        material,
        rd,
        md;
        CFL=cfl_target,
        relaxation_fraction=relaxation_fraction,
    )
    reported_cfl = compute_cfl(dt, material, rd, md)

    @printf("published Fourier-solver dt        = %.4e s\n", PUBLISHED_DT)
    @printf("DG element-width estimate          = %.4e m\n", dtinfo.h)
    @printf("DG coupled maximum speed           = %.2f m/s\n", dtinfo.cmax)
    @printf("DG wave-limited dt                 = %.4e s\n", dtinfo.dt_wave)
    @printf("DG relaxation-limited dt           = %.4e s\n", dtinfo.dt_relax)
    @printf("selected DG dt                     = %.4e s (%s limited)\n", dt, string(dtinfo.limiting))
    @printf("reported DG CFL                    = %.4f\n", reported_cfl)
    @printf("estimated time steps               = %d\n", ceil(Int, final_time / dt))

    tspan = (0.0, final_time)
    save_times = collect(range(tspan[1], tspan[2]; length=nsaves))
    prob = ODEProblem(rhs_thermoelastic_br1!, u0, tspan, params)

    progress_cb = make_fixed_step_progress_callback(
        tspan,
        dt;
        cfl=reported_cfl,
        description="Carcione Figure 6 - 2D",
        max_updates=progress_updates,
    )

    println("Starting SSPRK54 solve...")
    sol = solve(
        prob,
        SSPRK54();
        dt=dt,
        adaptive=false,
        saveat=save_times,
        save_everystep=false,
        save_start=true,
        save_end=true,
        callback=progress_cb,
        progress=false,
    )
    println("Solve complete with ", length(sol.t), " stored snapshots")

    for (index, state) in enumerate(sol.u)
        all(isfinite, state) ||
            error("nonfinite state at snapshot $index, t=$(sol.t[index])")
    end

    # ---------------------------------------------------------------------
    # Interpolation and optional VTK output
    # ---------------------------------------------------------------------
    xplot = rd.Vp * x
    zplot = rd.Vp * y
    ufinal = interpolate_state(rd, sol.u[end])

    # ---------------------------------------------------------------------
    # Longitudinal-versus-shear diagnostics
    # ---------------------------------------------------------------------
    # At the final snapshot, the S-wave radius and the high-frequency T-wave
    # radius are almost identical for the reference material. Radius alone
    # therefore cannot identify an S wave. Use polarization and vorticity.
    candidate_radius = properties.vS * effective_time
    diagnostic_half_width = max(7.5e-4, 4 * source_sigma)

    global_mode_diag = longitudinal_shear_diagnostics(
        sol.u[end],
        rd,
        md,
        material;
        rmin=3 * source_sigma,
        rmax=0.49 * DOMAIN_LENGTH,
    )
    candidate_mode_diag = longitudinal_shear_diagnostics(
        sol.u[end],
        rd,
        md,
        material;
        rmin=max(3 * source_sigma, candidate_radius - diagnostic_half_width),
        rmax=min(0.49 * DOMAIN_LENGTH, candidate_radius + diagnostic_half_width),
    )

    @printf(
        "Global tangential/radial velocity ratio = %.3e\n",
        global_mode_diag.tangential_to_radial,
    )
    @printf(
        "Global vorticity/dilatation ratio        = %.3e\n",
        global_mode_diag.vorticity_to_dilatation,
    )
    @printf(
        "Candidate-ring tangential/radial ratio   = %.3e\n",
        candidate_mode_diag.tangential_to_radial,
    )
    @printf(
        "Candidate-ring vorticity/dilatation      = %.3e\n",
        candidate_mode_diag.vorticity_to_dilatation,
    )

    write_mode_diagnostics(
        joinpath(outdir, "longitudinal_shear_diagnostic.txt"),
        global_mode_diag,
        candidate_mode_diag,
        candidate_radius,
    )

    if write_vtk
        vtk_ids = unique(round.(
            Int,
            range(1, length(sol.u); length=min(nvtk, length(sol.u))),
        ))
        println("Writing VTK snapshots...")
        for (counter, index) in enumerate(vtk_ids)
            uplot = interpolate_state(rd, sol.u[index])
            filename = joinpath(outdir, @sprintf("carcione_fig6_%04d.vtu", counter))
            export_quad_subcells_vtu(filename, xplot, zplot, uplot, material)
        end
    end

    # ---------------------------------------------------------------------
    # Figure 6-style snapshots
    # ---------------------------------------------------------------------
    vz = vec(@view ufinal[:, :, IVZ])
    theta = vec(@view ufinal[:, :, ITHETA])
    vz_norm, vz_scale = normalized_field(vz)
    theta_norm, theta_scale = normalized_field(theta)

    # The published grayscale limits are not reported.  Controlled clipping
    # reveals the weaker E-wave while the unscaled maxima remain in the titles.
    vz_display = clamp.(vz_norm ./ display_clip, -1.0, 1.0)
    theta_display = clamp.(theta_norm ./ display_clip, -1.0, 1.0)

    x_mm = 1e3 .* vec(xplot)
    z_mm = 1e3 .* vec(zplot)

    p_vz = scatter(
        x_mm,
        z_mm;
        marker_z=vz_display,
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
        title=@sprintf(
            "Figure 6: v_z/max|v_z| at %.3f us\nmax|v_z|=%.3e m/s",
            1e6 * sol.t[end],
            vz_scale,
        ),
    )

    p_theta = scatter(
        x_mm,
        z_mm;
        marker_z=theta_display,
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
        title=@sprintf(
            "Figure 6: theta/max|theta| at %.3f us\nmax|theta|=%.3e K",
            1e6 * sol.t[end],
            theta_scale,
        ),
    )

    # Figure 6 has a propagating E-wave and a strongly attenuated, diffusive
    # T mode.  Overlay only the finite-frequency E-wave prediction; drawing a
    # high-frequency T-wave circle would be physically misleading here.
    angles = range(0, 2pi; length=500)
    for plot_object in (p_vz, p_theta)
        plot!(
            plot_object,
            1e3 * expected_e_radius .* cos.(angles),
            1e3 * expected_e_radius .* sin.(angles);
            label=false,
            linewidth=1,
            linestyle=:dash,
        )
    end

    savefig(
        plot(p_vz, p_theta; layout=(1, 2), size=(1350, 570)),
        joinpath(outdir, "carcione_figure6_DG.png"),
    )

    # Plot radial and tangential velocity using the SAME amplitude scale.
    # Never normalize v_t by its own maximum: doing so makes tiny numerical
    # shear contamination look as large as the physical longitudinal wave.
    vx_plot = vec(@view ufinal[:, :, IVX])
    vz_plot = vec(@view ufinal[:, :, IVZ])
    x_plot_vec = vec(xplot)
    z_plot_vec = vec(zplot)
    r_plot = @. sqrt(x_plot_vec^2 + z_plot_vec^2)
    inv_r_plot = similar(r_plot)
    r_plot_tol = sqrt(eps(Float64)) * max(1.0, maximum(r_plot))
    @inbounds for i in eachindex(r_plot)
        inv_r_plot[i] = r_plot[i] > r_plot_tol ? inv(r_plot[i]) : 0.0
    end
    vr_plot = @. (x_plot_vec * vx_plot + z_plot_vec * vz_plot) * inv_r_plot
    vt_plot = @. (-z_plot_vec * vx_plot + x_plot_vec * vz_plot) * inv_r_plot
    vr_scale = maximum(abs, vr_plot)
    vr_scale = vr_scale > eps(Float64) ? vr_scale : 1.0

    p_vr = scatter(
        x_mm,
        z_mm;
        marker_z=clamp.(vr_plot ./ vr_scale, -1.0, 1.0),
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
        title="radial velocity / max|v_r|",
    )
    p_vt = scatter(
        x_mm,
        z_mm;
        marker_z=clamp.(vt_plot ./ vr_scale, -1.0, 1.0),
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
        title=@sprintf(
            "tangential velocity / max|v_r|\n||v_t||/||v_r||=%.2e",
            global_mode_diag.tangential_to_radial,
        ),
    )
    savefig(
        plot(p_vr, p_vt; layout=(1, 2), size=(1350, 570)),
        joinpath(outdir, "longitudinal_shear_decomposition.png"),
    )

    # ---------------------------------------------------------------------
    # Radial diagnostics
    # ---------------------------------------------------------------------
    r_v, profile_v = radial_rms_profile(
        vec(xplot),
        vec(zplot),
        vz;
        nbins=radial_bins,
    )
    r_t, profile_t = radial_rms_profile(
        vec(xplot),
        vec(zplot),
        theta;
        nbins=radial_bins,
    )

    open(joinpath(outdir, "radial_profiles.csv"), "w") do io
        println(io, "radius_m,vz_rms,theta_rms")
        for index in eachindex(r_v)
            @printf(
                io,
                "%.16e,%.16e,%.16e\n",
                r_v[index],
                profile_v[index],
                profile_t[index],
            )
        end
    end

    finite_v = filter(isfinite, profile_v)
    finite_t = filter(isfinite, profile_t)
    profile_v_scale = isempty(finite_v) ? 1.0 : maximum(finite_v)
    profile_t_scale = isempty(finite_t) ? 1.0 : maximum(finite_t)
    profile_v_scale = profile_v_scale > 0 ? profile_v_scale : 1.0
    profile_t_scale = profile_t_scale > 0 ? profile_t_scale : 1.0

    p_profile = plot(
        1e3 .* r_v,
        profile_v ./ profile_v_scale;
        linewidth=2,
        label="RMS v_z",
        xlabel="radius (mm)",
        ylabel="normalised radial RMS",
        title="Figure 6 radial wavefront diagnostic",
    )
    plot!(
        p_profile,
        1e3 .* r_t,
        profile_t ./ profile_t_scale;
        linewidth=2,
        label="RMS theta",
    )
    vline!(
        p_profile,
        [1e3 * expected_e_radius];
        linestyle=:dash,
        label="finite-frequency E prediction",
    )
    if isfinite(modes.thermal.attenuation_length)
        thermal_near_field = 3 * modes.thermal.attenuation_length
        vspan!(
            p_profile,
            [0.0, 1e3 * thermal_near_field];
            label="3 T attenuation lengths",
            alpha=0.12,
        )
    end
    savefig(p_profile, joinpath(outdir, "radial_profiles.png"))

    # ---------------------------------------------------------------------
    # Energy history and semi-discrete energy audit
    # ---------------------------------------------------------------------
    energy = write_energy_history(
        joinpath(outdir, "energy.csv"),
        sol,
        rd,
        md,
        material,
    )
    max_balance_residual = write_energy_balance(
        joinpath(outdir, "energy_balance.csv"),
        sol,
        params,
        rd,
        md,
        material,
        source,
    )
    @printf("maximum relative energy residual  = %.3e\n", max_balance_residual)

    total_energy = getproperty.(energy, :total)
    energy_scale = isempty(total_energy) ? 1.0 : maximum(total_energy)
    energy_scale = energy_scale > 0 ? energy_scale : 1.0

    p_energy = plot(
        1e6 .* sol.t,
        total_energy ./ energy_scale;
        linewidth=2,
        label="total",
        xlabel="time (us)",
        ylabel="energy/max energy",
        title="Discrete thermoelastic energy",
    )
    plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :elastic) ./ energy_scale; label="elastic")
    plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :kinetic) ./ energy_scale; label="kinetic")
    plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :temperature) ./ energy_scale; label="temperature")
    plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :heat_flux) ./ energy_scale; label="heat flux")
    savefig(p_energy, joinpath(outdir, "energy_history.png"))

    # ---------------------------------------------------------------------
    # Run metadata
    # ---------------------------------------------------------------------
    open(joinpath(outdir, "run_summary.txt"), "w") do io
        println(io, "Carcione et al. (2019), Figure 6 standalone DG experiment")
        println(io, "DOI: 10.1190/GEO2018-0448.1")
        println(io)
        println(io, "DG scheme: weak elastic stress / strong velocity")
        println(io, "Mechanical flux: Chan--Shukla scalar penalty flux")
        println(io, "Thermal scheme: weak temperature / strong Cattaneo heat flux, BR1 traces")
        println(io)
        @printf(io, "N = %d\n", N)
        @printf(io, "K1D = %d\n", K1D)
        @printf(io, "domain_length_m = %.16e\n", DOMAIN_LENGTH)
        @printf(io, "final_time_s = %.16e\n", final_time)
        @printf(io, "conductivity_W_per_mK = %.16e\n", conductivity)
        @printf(io, "tau_s = %.16e\n", material.tau)
        @printf(io, "source_frequency_parameter = %.16e\n", SOURCE_FREQUENCY)
        @printf(io, "actual_carrier_frequency_Hz = %.16e\n", carrier_frequency_hz)
        @printf(io, "source_delay_s = %.16e\n", SOURCE_DELAY)
        @printf(io, "source_sigma_m = %.16e\n", source_sigma)
        @printf(io, "source_amplitude = %.16e\n", source_amplitude)
        @printf(io, "source_filter_dt_s = %.16e\n", source.filter_dt)
        println(io, "source_wavelet_mode = ", wavelet_mode)
        @printf(io, "source_integral = %.16e\n", source_integral)
        @printf(io, "display_clip_fraction = %.16e\n", display_clip)
        @printf(io, "selected_dt_s = %.16e\n", dt)
        @printf(io, "target_CFL = %.16e\n", cfl_target)
        @printf(io, "reported_CFL = %.16e\n", reported_cfl)
        @printf(io, "finite_frequency_E_phase_speed_m_per_s = %.16e\n", modes.elastic.phase_velocity)
        @printf(io, "finite_frequency_E_attenuation_per_m = %.16e\n", modes.elastic.attenuation)
        @printf(io, "finite_frequency_T_phase_speed_m_per_s = %.16e\n", modes.thermal.phase_velocity)
        @printf(io, "finite_frequency_T_attenuation_per_m = %.16e\n", modes.thermal.attenuation)
        @printf(io, "finite_frequency_T_attenuation_length_m = %.16e\n", modes.thermal.attenuation_length)
        @printf(io, "source_peak_predicted_E_radius_m = %.16e\n", expected_e_radius)
        @printf(io, "formal_source_peak_T_radius_m = %.16e\n", formal_t_radius)
        @printf(io, "T_attenuation_exponent_at_formal_radius = %.16e\n", thermal_decay_exponent)
        @printf(io, "max_relative_energy_residual = %.16e\n", max_balance_residual)
        @printf(io, "global_tangential_to_radial_L2 = %.16e\n", global_mode_diag.tangential_to_radial)
        @printf(io, "global_tangential_kinetic_fraction = %.16e\n", global_mode_diag.tangential_kinetic_fraction)
        @printf(io, "global_vorticity_to_dilatation_L2 = %.16e\n", global_mode_diag.vorticity_to_dilatation)
        @printf(io, "candidate_ring_tangential_to_radial_L2 = %.16e\n", candidate_mode_diag.tangential_to_radial)
        @printf(io, "candidate_ring_vorticity_to_dilatation_L2 = %.16e\n", candidate_mode_diag.vorticity_to_dilatation)
        println(io)
        println(io, "The Figure 6 thermal mode is highly attenuated and should appear as a")
        println(io, "diffusive near-source feature, not as a clean high-frequency T-wave ring.")
        println(io, "The paper specifies the source waveform and location but not its absolute")
        println(io, "amplitude or spatial regularisation. These remain explicit run parameters.")
    end

    println("Saved Figure 6 comparison image: ", joinpath(outdir, "carcione_figure6_DG.png"))
    println("Saved longitudinal/shear plot: ", joinpath(outdir, "longitudinal_shear_decomposition.png"))
    println("Saved mode diagnostic: ", joinpath(outdir, "longitudinal_shear_diagnostic.txt"))
    println("Saved radial diagnostic: ", joinpath(outdir, "radial_profiles.png"))
    println("Saved energy history: ", joinpath(outdir, "energy.csv"))
    println("Saved energy balance: ", joinpath(outdir, "energy_balance.csv"))
    println("Saved source history: ", joinpath(outdir, "source_history.csv"))
    println("All outputs saved in: ", outdir)
    return outdir
end

end # module CarcioneFigure6Driver

Base.invokelatest(CarcioneFigure6Driver.main)

