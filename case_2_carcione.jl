# Standalone driver for Carcione et al. (2019), Figure 7.


module CarcioneFigure7Driver

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
const FIGURE7_CONDUCTIVITY = 4.5e6
const BODY_TEXT_CONDUCTIVITY = 1.13e6

parse_env(::Type{T}, name::AbstractString, default) where {T} =
    parse(T, get(ENV, name, string(default)))

function parse_bool_env(name::AbstractString, default::Bool)
    value = lowercase(strip(get(ENV, name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$name must be true/false, yes/no, on/off, or 1/0"))
end

function make_output_dir(; root=joinpath(@__DIR__, "output_figure7"), keep_last=8)
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

function figure7_material(
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
    radius = sqrt.(x .^ 2 .+ z .^ 2)
    rmax = maximum(radius)
    edges = collect(range(0.0, rmax; length=nbins + 1))
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])
    sums = zeros(Float64, nbins)
    counts = zeros(Int, nbins)

    for i in eachindex(radius)
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
    DG.NSTATE == 8 || error("Figure 7 driver expects the eight-field 2D state")

    final_time = parse_env(Float64, "FIG7_FINAL_TIME", PUBLISHED_FINAL_TIME)
    N = parse_env(Int, "FIG7_N", 4)
    # At 3.5 MHz the limiting T-wave wavelength is about 0.433 mm.
    # K1D=58 with N=4 gives only about one element per T wavelength; use
    # roughly two elements per wavelength.  The odd value 117 also places
    # the centered source inside the central element rather than at a
    # four-element vertex.
    default_K1D = N == 4 ? 117 : begin
        candidate = ceil(Int, 2 * (DOMAIN_LENGTH / REFERENCE_DX) / N)
        isodd(candidate) ? candidate : candidate + 1
    end
    K1D = parse_env(Int, "FIG7_K1D", default_K1D)
    cfl_target = parse_env(Float64, "FIG7_CFL", 0.20)
    penalty_scale = parse_env(Float64, "FIG7_PENALTY_SCALE", 1.0)
    br1_penalty_scale = parse_env(Float64, "FIG7_BR1_PENALTY_SCALE", 0.0)
    conductivity = parse_env(Float64, "FIG7_CONDUCTIVITY", FIGURE7_CONDUCTIVITY)
    # The source is centered in the paper, but its spatial regularization is
    # not reported.  Use a narrow Gaussian so that the short T-wave spectrum
    # is not strongly suppressed.
    source_sigma = parse_env(Float64, "FIG7_SOURCE_SIGMA", 0.5 * REFERENCE_DX)
    source_amplitude = parse_env(Float64, "FIG7_SOURCE_AMPLITUDE", 1.0)
    nsaves = parse_env(Int, "FIG7_NSAVES", 41)
    nvtk = parse_env(Int, "FIG7_NVTK", 7)
    write_vtk = parse_bool_env("FIG7_WRITE_VTK", true)
    display_clip = parse_env(Float64, "FIG7_DISPLAY_CLIP", 0.20)
    wavelet_name = lowercase(strip(get(ENV, "FIG7_WAVELET_MODE", "cyclic_hz")))
    wavelet_mode = if wavelet_name == "printed"
        :printed
    elseif wavelet_name in ("cyclic_hz", "hz", "frequency")
        :cyclic_hz
    else
        error("FIG7_WAVELET_MODE must be printed or cyclic_hz")
    end

    final_time > 0 || error("FIG7_FINAL_TIME must be positive")
    N >= 1 || error("FIG7_N must be at least 1")
    K1D >= 1 || error("FIG7_K1D must be at least 1")
    cfl_target > 0 || error("FIG7_CFL must be positive")
    conductivity > 0 || error("FIG7_CONDUCTIVITY must be positive")
    source_sigma > 0 || error("FIG7_SOURCE_SIGMA must be positive")
    nsaves >= 2 || error("FIG7_NSAVES must be at least 2")
    0 < display_clip <= 1 || error("FIG7_DISPLAY_CLIP must lie in (0,1]")

    outdir = make_output_dir()
    println("Carcione Figure 7 standalone DG experiment")
    println("Output directory: ", outdir)

    rd = RefElemData(Quad(), SBP(), N)
    VXY_ref, EToV = uniform_mesh(Quad(), K1D)
    half_length = DOMAIN_LENGTH / 2
    VXY = (half_length .* VXY_ref[1], half_length .* VXY_ref[2])
    md = MeshData(VXY, EToV, rd)
    md = make_periodic(md)
    (; x, y) = md

    material = figure7_material(
        conductivity;
        penalty_scale=penalty_scale,
        br1_penalty_scale=br1_penalty_scale,
    )
    properties = carcione_reference_properties(material)

    effective_time = max(final_time - SOURCE_DELAY, 0.0)
    expected_e_radius = properties.vEinf * effective_time
    expected_t_radius = properties.vTinf * effective_time

    @printf("N                                  = %d\n", N)
    @printf("elements per coordinate            = %d\n", K1D)
    @printf("quadrilateral elements             = %d\n", md.K)
    @printf("nominal nodal spacing              = %.4e m\n", DOMAIN_LENGTH / (K1D * N))
    @printf("Figure 7 conductivity gamma        = %.6e W/(m K)\n", conductivity)
    @printf("body-text conflicting value        = %.6e W/(m K)\n", BODY_TEXT_CONDUCTIVITY)
    @printf("thermal relaxation time tau        = %.6e s\n", material.tau)
    @printf("V_I                                = %.2f m/s\n", properties.vI)
    @printf("V_S                                = %.2f m/s\n", properties.vS)
    @printf("V_A                                = %.2f m/s\n", properties.vA)
    @printf("V_Einf                             = %.2f m/s\n", properties.vEinf)
    @printf("V_Tinf                             = %.2f m/s\n", properties.vTinf)
    @printf("source-peak E radius prediction    = %.3f mm\n", 1e3 * expected_e_radius)
    @printf("source-peak T radius prediction    = %.3f mm\n", 1e3 * expected_t_radius)
    @printf("source wavelet convention          = %s\n", String(wavelet_mode))

    u0 = zeros(eltype(x), size(x, 1), size(x, 2), NSTATE)
    source_shape = normalized_gaussian_source(
        rd,
        md;
        x0=0.0,
        z0=0.0,
        sigma=source_sigma,
    )

    source_filter_dt = parse_env(
        Float64,
        "FIG7_SOURCE_FILTER_DT",
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
        wavelet_mode=wavelet_mode,
    )

    cache = RHSCacheTE(u0, rd)
    params = (; rd, md, material, cache, source)

    source_integral = sum(md.wJq .* (rd.Vq * source_shape))
    @printf("discrete source integral           = %.16e\n", source_integral)
    isapprox(source_integral, 1.0; atol=1e-11, rtol=1e-11) ||
        error("normalized source integral is $source_integral rather than one")
    write_source_history(joinpath(outdir, "source_history.csv"), source, final_time, nsaves)

    dt, dtinfo = estimate_dt(
        material,
        rd,
        md;
        CFL=cfl_target,
        relaxation_fraction=0.5,
    )
    reported_cfl = compute_cfl(dt, material, rd, md)
    @printf("published Fourier-solver dt        = %.4e s\n", PUBLISHED_DT)
    @printf("DG wave-limited dt                 = %.4e s\n", dtinfo.dt_wave)
    @printf("DG relaxation-limited dt           = %.4e s\n", dtinfo.dt_relax)
    @printf("selected DG dt                     = %.4e s\n", dt)
    @printf("reported DG CFL                    = %.4f\n", reported_cfl)
    @printf("estimated time steps               = %d\n", ceil(Int, final_time / dt))

    tspan = (0.0, final_time)
    save_times = collect(range(tspan[1], tspan[2]; length=nsaves))
    prob = ODEProblem(rhs_thermoelastic_br1!, u0, tspan, params)

    # `reported_cfl` was already computed with the current DG API:
    #     compute_cfl(dt, material, rd, md)
    progress_cb = make_fixed_step_progress_callback(
        tspan,
        dt;
        cfl=reported_cfl,
        description="Carcione Figure 7 - 2D",
        max_updates=200,
    )

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


    for (i, U) in enumerate(sol.u)
        all(isfinite, U) || error("nonfinite state at snapshot $i, t=$(sol.t[i])")
    end

    xplot = rd.Vp * x
    zplot = rd.Vp * y
    ufinal = interpolate_state(rd, sol.u[end])

    if write_vtk
        vtk_ids = unique(round.(Int, range(1, length(sol.u); length=min(nvtk, length(sol.u)))))
        for (counter, index) in enumerate(vtk_ids)
            uplot = interpolate_state(rd, sol.u[index])
            filename = joinpath(outdir, @sprintf("carcione_fig7_%04d.vtu", counter))
            export_quad_subcells_vtu(filename, xplot, zplot, uplot, material)
        end
    end

    vz = vec(@view ufinal[:, :, IVZ])
    theta = vec(@view ufinal[:, :, ITHETA])
    vz_norm, vz_scale = normalized_field(vz)
    theta_norm, theta_scale = normalized_field(theta)

    # The paper does not report the grayscale limits.  Keep the true maxima
    # in the title and create a controlled contrast-enhanced display so the
    # weaker E-wave is not hidden by the stronger T-wave.
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
        title=@sprintf("Figure 7: v_z/max|v_z| at %.3f us\nmax|v_z|=%.3e m/s", 1e6 * sol.t[end], vz_scale),
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
        title=@sprintf("Figure 7: theta/max|theta| at %.3f us\nmax|theta|=%.3e K", 1e6 * sol.t[end], theta_scale),
    )

    angles = range(0, 2pi; length=500)
    for p in (p_vz, p_theta)
        plot!(p, 1e3 * expected_e_radius .* cos.(angles), 1e3 * expected_e_radius .* sin.(angles); label=false, linewidth=1, linestyle=:dash)
        plot!(p, 1e3 * expected_t_radius .* cos.(angles), 1e3 * expected_t_radius .* sin.(angles); label=false, linewidth=1, linestyle=:dot)
    end
    savefig(plot(p_vz, p_theta; layout=(1, 2), size=(1350, 570)), joinpath(outdir, "carcione_figure7_DG.png"))

    r_v, profile_v = radial_rms_profile(vec(xplot), vec(zplot), vz)
    r_t, profile_t = radial_rms_profile(vec(xplot), vec(zplot), theta)
    open(joinpath(outdir, "radial_profiles.csv"), "w") do io
        println(io, "radius_m,vz_rms,theta_rms")
        for i in eachindex(r_v)
            @printf(io, "%.16e,%.16e,%.16e\n", r_v[i], profile_v[i], profile_t[i])
        end
    end

    finite_v = filter(isfinite, profile_v)
    finite_t = filter(isfinite, profile_t)
    pv_scale = isempty(finite_v) ? 1.0 : maximum(finite_v)
    pt_scale = isempty(finite_t) ? 1.0 : maximum(finite_t)
    pv_scale = pv_scale > 0 ? pv_scale : 1.0
    pt_scale = pt_scale > 0 ? pt_scale : 1.0
    p_profile = plot(
        1e3 .* r_v,
        profile_v ./ pv_scale;
        linewidth=2,
        label="RMS v_z",
        xlabel="radius (mm)",
        ylabel="normalized radial RMS",
        title="Figure 7 radial wavefront diagnostic",
    )
    plot!(p_profile, 1e3 .* r_t, profile_t ./ pt_scale; linewidth=2, label="RMS theta")
    vline!(p_profile, [1e3 * expected_e_radius]; linestyle=:dash, label="V_Einf prediction")
    vline!(p_profile, [1e3 * expected_t_radius]; linestyle=:dot, label="V_Tinf prediction")
    savefig(p_profile, joinpath(outdir, "radial_profiles.png"))

    energy = write_energy_history(joinpath(outdir, "energy.csv"), sol, rd, md, material)
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
    Escale = isempty(total_energy) ? 1.0 : maximum(total_energy)
    Escale = Escale > 0 ? Escale : 1.0
    p_energy = plot(
        1e6 .* sol.t,
        total_energy ./ Escale;
        linewidth=2,
        label="total",
        xlabel="time (us)",
        ylabel="energy/max energy",
        title="Discrete thermoelastic energy",
    )
    plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :elastic) ./ Escale; label="elastic")
    plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :kinetic) ./ Escale; label="kinetic")
    plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :temperature) ./ Escale; label="temperature")
    plot!(p_energy, 1e6 .* sol.t, getproperty.(energy, :heat_flux) ./ Escale; label="heat flux")
    savefig(p_energy, joinpath(outdir, "energy_history.png"))

    open(joinpath(outdir, "run_summary.txt"), "w") do io
        println(io, "Carcione et al. (2019), Figure 7 standalone DG experiment")
        println(io, "DOI: 10.1190/GEO2018-0448.1")
        println(io)
        println(io, "The default conductivity follows the Figure 7 caption and Figure 2: 4.5e6 W/(m K).")
        println(io, "The paper body contains a conflicting 1.13e6 value; use FIG7_CONDUCTIVITY to test it.")
        println(io)
        @printf(io, "N = %d\n", N)
        @printf(io, "K1D = %d\n", K1D)
        @printf(io, "domain_length_m = %.16e\n", DOMAIN_LENGTH)
        @printf(io, "final_time_s = %.16e\n", final_time)
        @printf(io, "conductivity_W_per_mK = %.16e\n", conductivity)
        @printf(io, "tau_s = %.16e\n", material.tau)
        @printf(io, "source_frequency_Hz = %.16e\n", SOURCE_FREQUENCY)
        @printf(io, "source_delay_s = %.16e\n", SOURCE_DELAY)
        @printf(io, "source_sigma_m = %.16e\n", source_sigma)
        @printf(io, "source_amplitude = %.16e\n", source_amplitude)
        @printf(io, "display_clip_fraction = %.16e\n", display_clip)
        println(io, "source_wavelet_mode = ", wavelet_mode)
        @printf(io, "source_integral = %.16e\n", source_integral)
        @printf(io, "selected_dt_s = %.16e\n", dt)
        @printf(io, "target_CFL = %.16e\n", cfl_target)
        @printf(io, "reported_CFL = %.16e\n", reported_cfl)
        @printf(io, "V_Einf_m_per_s = %.16e\n", properties.vEinf)
        @printf(io, "V_Tinf_m_per_s = %.16e\n", properties.vTinf)
        @printf(io, "source_peak_predicted_E_radius_m = %.16e\n", expected_e_radius)
        @printf(io, "source_peak_predicted_T_radius_m = %.16e\n", expected_t_radius)
        @printf(io, "max_relative_energy_residual = %.16e\n", max_balance_residual)
    end

    println("Saved Figure 7 comparison image: ", joinpath(outdir, "carcione_figure7_DG.png"))
    println("Saved radial diagnostic: ", joinpath(outdir, "radial_profiles.png"))
    println("All outputs saved in: ", outdir)
    return outdir
end

end # module CarcioneFigure7Driver

Base.invokelatest(CarcioneFigure7Driver.main)


