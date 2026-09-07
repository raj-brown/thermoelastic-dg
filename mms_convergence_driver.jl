# Standalone manufactured-solution convergence driver.
#
# Place this file beside main.jl in the repository root:
#
#   repo/
#     main.jl
#     mms_convergence_driver.jl   <-- this file
#     src/DG.jl
#
# Run with:
#
#   julia --project=. mms_convergence_driver.jl
#
# This script does not include, modify, or call main.jl. It loads src/DG.jl
# directly and writes all results under mms_output/.
#
# Required DG.jl API (the proposed eight-field scheme):
#   NSTATE == 8
#   state [s_xx,s_zz,s_xz,v_x,v_z,theta,q_x,q_z]
#   make_material, RHSCacheTE, rhs_thermoelastic_br1!, estimate_dt
#   CarcioneHeatSource (used only to construct an identically zero source)

module MMSConvergenceDriver

ENV["GKSwstype"] = "100"

using Dates
using LinearAlgebra
using OrdinaryDiffEq
using OrdinaryDiffEqSSPRK
using Plots
using Printf
using StartUpDG

const REPO_ROOT = @__DIR__
const DG_FILE = joinpath(REPO_ROOT, "src", "DG.jl")
isfile(DG_FILE) || error("Cannot find src/DG.jl relative to $(REPO_ROOT)")

include(DG_FILE)
using .DG

gr()

# -----------------------------------------------------------------------------
# Required API check
# -----------------------------------------------------------------------------
const REQUIRED_DG_SYMBOLS = (
    :NSTATE,
    :STATE_FIELD_NAMES,
    :ISXX,
    :ISZZ,
    :ISXZ,
    :IVX,
    :IVZ,
    :ITHETA,
    :IQX,
    :IQZ,
    :make_material,
    :RHSCacheTE,
    :rhs_thermoelastic_br1!,
    :estimate_dt,
    :maximum_characteristic_speed,
    :CarcioneHeatSource,
)

for symbol in REQUIRED_DG_SYMBOLS
    isdefined(DG, symbol) || error(
        "src/DG.jl does not define DG.$symbol. " *
        "This driver requires the revised eight-field weak-stress/strong-velocity " *
        "thermoelastic-Cattaneo implementation."
    )
end

DG.NSTATE == 8 || error(
    "This driver expects eight fields " *
    "[s_xx,s_zz,s_xz,v_x,v_z,theta,q_x,q_z], but DG.NSTATE=$(DG.NSTATE)."
)

# -----------------------------------------------------------------------------
# Environment input helpers
# -----------------------------------------------------------------------------
parse_env(::Type{T}, name::AbstractString, default) where {T} =
    parse(T, get(ENV, name, string(default)))

function parse_int_list(name::AbstractString, default::AbstractString)
    text = get(ENV, name, default)
    values = parse.(Int, strip.(split(text, ',')))
    isempty(values) && error("$name must contain at least one integer")
    all(>(0), values) || error("all entries in $name must be positive")
    return sort(unique(values))
end

function make_output_dir(; root=joinpath(REPO_ROOT, "mms_output"))
    mkpath(root)
    stamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS_sss")
    outdir = joinpath(root, stamp)
    mkpath(outdir)
    return outdir
end

# -----------------------------------------------------------------------------
# Smooth periodic manufactured solution on Omega = [-1,1]^2
# -----------------------------------------------------------------------------
struct MMSData{T}
    kx::T
    kz::T
    omega::T
    phase_x::T
    phase_z::T
    amplitudes::NTuple{8,T}
end

function MMSData(;
    kx::Real=pi,
    kz::Real=pi,
    omega::Real=1.3,
    phase_x::Real=0.23,
    phase_z::Real=-0.31,
    amplitudes=(0.8, -0.6, 0.5, 0.7, -0.9, 0.4, 0.3, -0.2),
)
    length(amplitudes) == 8 || error("MMS amplitudes must contain eight entries")
    T = promote_type(
        Float64,
        typeof(kx),
        typeof(kz),
        typeof(omega),
        typeof(phase_x),
        typeof(phase_z),
        map(typeof, amplitudes)...,
    )
    amps = ntuple(i -> convert(T, amplitudes[i]), 8)
    return MMSData(
        convert(T, kx),
        convert(T, kz),
        convert(T, omega),
        convert(T, phase_x),
        convert(T, phase_z),
        amps,
    )
end

"""
Fill `u` with the exact manufactured state

  s_xx  = a1 sin(X) cos(Z) cos(omega*t)
  s_zz  = a2 cos(X) sin(Z) sin(omega*t)
  s_xz  = a3 sin(X) sin(Z) cos(omega*t)
  v_x   = a4 cos(X) cos(Z) sin(omega*t)
  v_z   = a5 sin(X) cos(Z) cos(omega*t)
  theta = a6 cos(X) sin(Z) cos(omega*t)
  q_x   = a7 sin(X) sin(Z) sin(omega*t)
  q_z   = a8 cos(X) cos(Z) cos(omega*t),

where X=kx*x+phase_x and Z=kz*z+phase_z.
"""
function exact_state!(u, x, z, time, mms::MMSData)
    size(x) == size(z) || throw(DimensionMismatch("x and z must have equal size"))
    size(u, 1) == size(x, 1) || throw(DimensionMismatch("u and x point counts differ"))
    size(u, 2) == size(x, 2) || throw(DimensionMismatch("u and x element counts differ"))
    size(u, 3) == DG.NSTATE || throw(DimensionMismatch("u must have $(DG.NSTATE) fields"))

    a1, a2, a3, a4, a5, a6, a7, a8 = mms.amplitudes
    st, ct = sincos(mms.omega * time)

    @inbounds for k in axes(x, 2), i in axes(x, 1)
        sx, cx = sincos(mms.kx * x[i, k] + mms.phase_x)
        sz, cz = sincos(mms.kz * z[i, k] + mms.phase_z)

        u[i, k, DG.ISXX]   = a1 * sx * cz * ct
        u[i, k, DG.ISZZ]   = a2 * cx * sz * st
        u[i, k, DG.ISXZ]   = a3 * sx * sz * ct
        u[i, k, DG.IVX]    = a4 * cx * cz * st
        u[i, k, DG.IVZ]    = a5 * sx * cz * ct
        u[i, k, DG.ITHETA] = a6 * cx * sz * ct
        u[i, k, DG.IQX]    = a7 * sx * sz * st
        u[i, k, DG.IQZ]    = a8 * cx * cz * ct
    end

    return u
end

# -----------------------------------------------------------------------------
# O(1) anisotropic material for MMS verification
# -----------------------------------------------------------------------------
function make_mms_material(; penalty_scale=1.0, br1_penalty_scale=0.25)
    C = [
        4.0   1.0    0.25
        1.0   3.0   -0.15
        0.25 -0.15   1.20
    ]

    b = [0.30, -0.20, 0.10]

    K = [
        0.80 0.12
        0.12 0.60
    ]

    return DG.make_material(
        rho=1.30,
        C=C,
        b=b,
        heat_capacity=1.70,
        Tref=1.20,
        K=K,
        tau=0.40,
        penalty_scale=penalty_scale,
        br1_penalty_scale=br1_penalty_scale,
    )
end

# -----------------------------------------------------------------------------
# Exact forcing for the complete first-order thermoelastic-Cattaneo system
# -----------------------------------------------------------------------------
"""
Add forcing so that `exact_state!` satisfies

  s_t     = C E(v)                              + f_s,
  v_t     = rho^(-1) div(s - b theta)           + f_v,
  theta_t = -(Tref b^T E(v) + div(q))/c         + f_theta,
  q_t     = -(q + K grad(theta))/tau            + f_q.

The forcing is added in state-rate form, so it can be applied after the
homogeneous DG operator without changing any numerical flux.
"""
function add_mms_forcing!(du, x, z, time, material, mms::MMSData)
    C11 = material.C[1, 1]
    C12 = material.C[1, 2]
    C13 = material.C[1, 3]
    C22 = material.C[2, 2]
    C23 = material.C[2, 3]
    C33 = material.C[3, 3]

    b1, b2, b3 = material.b
    K11 = material.K[1, 1]
    K12 = material.K[1, 2]
    K21 = material.K[2, 1]
    K22 = material.K[2, 2]

    rho = material.rho
    heat_capacity = material.heat_capacity
    Tref = material.Tref
    tau = material.tau

    kx = mms.kx
    kz = mms.kz
    omega = mms.omega
    a1, a2, a3, a4, a5, a6, a7, a8 = mms.amplitudes
    st, ct = sincos(omega * time)

    @inbounds for k in axes(x, 2), i in axes(x, 1)
        sx, cx = sincos(kx * x[i, k] + mms.phase_x)
        sz, cz = sincos(kz * z[i, k] + mms.phase_z)

        # Exact state values used by the Cattaneo source terms.
        qx = a7 * sx * sz * st
        qz = a8 * cx * cz * ct

        # Exact time derivatives.
        sxx_t = -a1 * omega * sx * cz * st
        szz_t =  a2 * omega * cx * sz * ct
        sxz_t = -a3 * omega * sx * sz * st
        vx_t =    a4 * omega * cx * cz * ct
        vz_t =   -a5 * omega * sx * cz * st
        theta_t = -a6 * omega * cx * sz * st
        qx_t =    a7 * omega * sx * sz * ct
        qz_t =   -a8 * omega * cx * cz * st

        # Velocity derivatives and engineering-shear strain rate.
        vx_x = -a4 * kx * sx * cz * st
        vx_z = -a4 * kz * cx * sz * st
        vz_x =  a5 * kx * cx * cz * ct
        vz_z = -a5 * kz * sx * sz * ct

        epsxx = vx_x
        epszz = vz_z
        gammaxz = vx_z + vz_x

        # Elastic-stress and temperature derivatives.
        sxx_x = a1 * kx * cx * cz * ct
        szz_z = a2 * kz * cx * cz * st
        sxz_x = a3 * kx * cx * sz * ct
        sxz_z = a3 * kz * sx * cz * ct
        theta_x = -a6 * kx * sx * sz * ct
        theta_z =  a6 * kz * cx * cz * ct

        # Total stress Sigma = s - b*theta.
        div_sigma_x =
            (sxx_x - b1 * theta_x) +
            (sxz_z - b3 * theta_z)

        div_sigma_z =
            (sxz_x - b3 * theta_x) +
            (szz_z - b2 * theta_z)

        # Heat-flux divergence.
        qx_x =  a7 * kx * cx * sz * st
        qz_z = -a8 * kz * cx * sz * ct
        div_q = qx_x + qz_z

        b_dot_eps = b1 * epsxx + b2 * epszz + b3 * gammaxz

        # C E(v) in engineering-shear Voigt ordering.
        Ceps_xx = C11 * epsxx + C12 * epszz + C13 * gammaxz
        Ceps_zz = C12 * epsxx + C22 * epszz + C23 * gammaxz
        Ceps_xz = C13 * epsxx + C23 * epszz + C33 * gammaxz

        du[i, k, DG.ISXX] += sxx_t - Ceps_xx
        du[i, k, DG.ISZZ] += szz_t - Ceps_zz
        du[i, k, DG.ISXZ] += sxz_t - Ceps_xz

        du[i, k, DG.IVX] += vx_t - div_sigma_x / rho
        du[i, k, DG.IVZ] += vz_t - div_sigma_z / rho

        du[i, k, DG.ITHETA] +=
            theta_t + (Tref * b_dot_eps + div_q) / heat_capacity

        du[i, k, DG.IQX] +=
            qx_t + (qx + K11 * theta_x + K12 * theta_z) / tau

        du[i, k, DG.IQZ] +=
            qz_t + (qz + K21 * theta_x + K22 * theta_z) / tau
    end

    return nothing
end

function rhs_mms!(du, u, parameters, time)
    DG.rhs_thermoelastic_br1!(du, u, parameters, time)
    add_mms_forcing!(
        du,
        parameters.md.x,
        parameters.md.y,
        time,
        parameters.material,
        parameters.mms,
    )
    return nothing
end

# -----------------------------------------------------------------------------
# Error calculation with an independent overintegrated Gauss rule
# -----------------------------------------------------------------------------
#
# Do not measure an SBP/GLL solution with the same collocated Lobatto rule used
# by the solver.  For SBP quadrilaterals rd.Vq may be I and md.wJq contains the
# Lobatto weights.  That produces a discrete nodal norm and can exaggerate
# odd/even cancellation.  The helpers below construct an independent tensor
# Gauss-Legendre rule and interpolate the degree-N DG polynomial to it.

function gauss_legendre_rule(n::Int)
    n >= 1 || throw(ArgumentError("Gauss rule must contain at least one point"))
    if n == 1
        return [0.0], [2.0]
    end

    beta = [j / sqrt(4j^2 - 1) for j in 1:(n - 1)]
    F = eigen(SymTridiagonal(zeros(n), beta))
    nodes = collect(F.values)
    weights = 2.0 .* vec(F.vectors[1, :]).^2
    return nodes, weights
end

function unique_sorted_nodes(values)
    values = sort(collect(float.(values)))
    isempty(values) && return values

    scale = max(maximum(abs, values), 1.0)
    tolerance = 500 * eps(eltype(values)) * scale
    nodes = [values[1]]

    for value in @view values[2:end]
        if abs(value - nodes[end]) > tolerance
            push!(nodes, value)
        end
    end
    return nodes
end

function lagrange_interpolation_matrix(nodes, evaluation_points)
    n = length(nodes)
    m = length(evaluation_points)
    n >= 1 || throw(ArgumentError("interpolation node set is empty"))

    barycentric_weights = ones(promote_type(eltype(nodes), eltype(evaluation_points)), n)
    for j in 1:n
        for ell in 1:n
            if ell != j
                barycentric_weights[j] /= nodes[j] - nodes[ell]
            end
        end
    end

    V = zeros(eltype(barycentric_weights), m, n)
    scale = max(maximum(abs, nodes), 1.0)
    tolerance = 500 * eps(eltype(V)) * scale

    for i in 1:m
        x = evaluation_points[i]
        matched = findfirst(j -> abs(x - nodes[j]) <= tolerance, 1:n)

        if matched !== nothing
            V[i, matched] = one(eltype(V))
            continue
        end

        denominator = zero(eltype(V))
        for j in 1:n
            denominator += barycentric_weights[j] / (x - nodes[j])
        end
        for j in 1:n
            V[i, j] = (barycentric_weights[j] / (x - nodes[j])) / denominator
        end
    end
    return V
end

function reference_nodes_2d(rd)
    if hasproperty(rd, :rst)
        rst = getproperty(rd, :rst)
        return vec(rst[1]), vec(rst[2])
    elseif hasproperty(rd, :r) && hasproperty(rd, :s)
        return vec(getproperty(rd, :r)), vec(getproperty(rd, :s))
    else
        error("Cannot obtain the two-dimensional reference nodes from RefElemData")
    end
end

function overintegrated_error_rule(rd, md)
    default_n1d = max(2 * rd.N + 3, 8)
    n1d = parse(Int, get(ENV, "MMS_ERROR_QUAD_1D", string(default_n1d)))
    n1d >= rd.N + 2 || error(
        "MMS_ERROR_QUAD_1D must be at least N+2=$(rd.N + 2); got $n1d",
    )

    gauss_nodes, gauss_weights = gauss_legendre_rule(n1d)
    rnodes_flat, snodes_flat = reference_nodes_2d(rd)
    rnodes = unique_sorted_nodes(rnodes_flat)
    snodes = unique_sorted_nodes(snodes_flat)

    length(rnodes) == rd.N + 1 || error(
        "expected $(rd.N + 1) unique r nodes, found $(length(rnodes))",
    )
    length(snodes) == rd.N + 1 || error(
        "expected $(rd.N + 1) unique s nodes, found $(length(snodes))",
    )

    Lr = lagrange_interpolation_matrix(rnodes, gauss_nodes)
    Ls = lagrange_interpolation_matrix(snodes, gauss_nodes)

    Np = length(rnodes_flat)
    Nq = n1d^2
    Verror = zeros(promote_type(eltype(Lr), eltype(Ls)), Nq, Np)
    rq = zeros(eltype(Verror), Nq)
    sq = similar(rq)
    wq = similar(rq)

    node_r_index = [argmin(abs.(rnodes .- r)) for r in rnodes_flat]
    node_s_index = [argmin(abs.(snodes .- s)) for s in snodes_flat]

    row = 0
    for js in 1:n1d, ir in 1:n1d
        row += 1
        rq[row] = gauss_nodes[ir]
        sq[row] = gauss_nodes[js]
        wq[row] = gauss_weights[ir] * gauss_weights[js]

        for node in 1:Np
            Verror[row, node] =
                Lr[ir, node_r_index[node]] * Ls[js, node_s_index[node]]
        end
    end

    xq = Verror * md.x
    zq = Verror * md.y
    Jq = Verror * md.J
    wJq = reshape(wq, :, 1) .* abs.(Jq)

    return (; Verror, xq, zq, wJq, n1d)
end

function mms_errors(u, rd, md, material, mms, time)
    error_rule = overintegrated_error_rule(rd, md)
    Verror = error_rule.Verror
    xq = error_rule.xq
    zq = error_rule.zq
    wJq = error_rule.wJq

    Nq = size(Verror, 1)
    K = size(u, 2)
    size(wJq) == (Nq, K) || throw(DimensionMismatch(
        "overintegrated weights have size $(size(wJq)); expected $((Nq, K))",
    ))

    Tq = promote_type(eltype(u), eltype(wJq))
    numerical = Array{Tq}(undef, Nq, K, DG.NSTATE)
    exact = similar(numerical)

    @views for field in 1:DG.NSTATE
        mul!(numerical[:, :, field], Verror, u[:, :, field])
    end
    exact_state!(exact, xq, zq, time, mms)

    component_abs = zeros(Tq, DG.NSTATE)
    component_exact = similar(component_abs)

    @inbounds for field in 1:DG.NSTATE
        error_squared = zero(Tq)
        exact_squared = zero(Tq)

        for k in axes(wJq, 2), i in axes(wJq, 1)
            weight = wJq[i, k]
            difference = numerical[i, k, field] - exact[i, k, field]
            error_squared += weight * difference^2
            exact_squared += weight * exact[i, k, field]^2
        end

        component_abs[field] = sqrt(max(error_squared, zero(Tq)))
        component_exact[field] = sqrt(max(exact_squared, zero(Tq)))
    end

    component_rel = similar(component_abs)
    @inbounds for field in 1:DG.NSTATE
        component_rel[field] = component_abs[field] /
            max(component_exact[field], eps(Tq))
    end

    combined_abs = norm(component_abs)
    combined_exact = norm(component_exact)
    combined_rel = combined_abs / max(combined_exact, eps(Tq))

    S = material.S
    Kinv = material.Kinv
    energy_error_squared = zero(Tq)
    energy_exact_squared = zero(Tq)

    @inbounds for k in axes(wJq, 2), i in axes(wJq, 1)
        weight = wJq[i, k]

        es1 = numerical[i, k, DG.ISXX] - exact[i, k, DG.ISXX]
        es2 = numerical[i, k, DG.ISZZ] - exact[i, k, DG.ISZZ]
        es3 = numerical[i, k, DG.ISXZ] - exact[i, k, DG.ISXZ]
        ev1 = numerical[i, k, DG.IVX] - exact[i, k, DG.IVX]
        ev2 = numerical[i, k, DG.IVZ] - exact[i, k, DG.IVZ]
        et = numerical[i, k, DG.ITHETA] - exact[i, k, DG.ITHETA]
        eq1 = numerical[i, k, DG.IQX] - exact[i, k, DG.IQX]
        eq2 = numerical[i, k, DG.IQZ] - exact[i, k, DG.IQZ]

        stress_error =
            es1 * (S[1, 1] * es1 + S[1, 2] * es2 + S[1, 3] * es3) +
            es2 * (S[2, 1] * es1 + S[2, 2] * es2 + S[2, 3] * es3) +
            es3 * (S[3, 1] * es1 + S[3, 2] * es2 + S[3, 3] * es3)

        flux_error =
            eq1 * (Kinv[1, 1] * eq1 + Kinv[1, 2] * eq2) +
            eq2 * (Kinv[2, 1] * eq1 + Kinv[2, 2] * eq2)

        energy_error_squared += weight * (
            stress_error +
            material.rho * (ev1^2 + ev2^2) +
            material.heat_capacity / material.Tref * et^2 +
            material.tau / material.Tref * flux_error
        )

        xs1 = exact[i, k, DG.ISXX]
        xs2 = exact[i, k, DG.ISZZ]
        xs3 = exact[i, k, DG.ISXZ]
        xv1 = exact[i, k, DG.IVX]
        xv2 = exact[i, k, DG.IVZ]
        xt = exact[i, k, DG.ITHETA]
        xq1 = exact[i, k, DG.IQX]
        xq2 = exact[i, k, DG.IQZ]

        stress_exact =
            xs1 * (S[1, 1] * xs1 + S[1, 2] * xs2 + S[1, 3] * xs3) +
            xs2 * (S[2, 1] * xs1 + S[2, 2] * xs2 + S[2, 3] * xs3) +
            xs3 * (S[3, 1] * xs1 + S[3, 2] * xs2 + S[3, 3] * xs3)

        flux_exact =
            xq1 * (Kinv[1, 1] * xq1 + Kinv[1, 2] * xq2) +
            xq2 * (Kinv[2, 1] * xq1 + Kinv[2, 2] * xq2)

        energy_exact_squared += weight * (
            stress_exact +
            material.rho * (xv1^2 + xv2^2) +
            material.heat_capacity / material.Tref * xt^2 +
            material.tau / material.Tref * flux_exact
        )
    end

    energy_abs = sqrt(max(energy_error_squared, zero(Tq)))
    energy_exact = sqrt(max(energy_exact_squared, zero(Tq)))
    energy_rel = energy_abs / max(energy_exact, eps(Tq))

    return (;
        component_abs,
        component_rel,
        combined_abs,
        combined_rel,
        energy_abs,
        energy_rel,
        error_quadrature_points_1d=error_rule.n1d,
    )
end

# -----------------------------------------------------------------------------
# One N/K solve
# -----------------------------------------------------------------------------
function run_case(
    N::Int,
    K1D::Int,
    final_time,
    cfl,
    material,
    mms;
    time_scheme::AbstractString,
    tolerance,
    fixed_dt_fraction,
)
    rd = RefElemData(Quad(), SBP(), N)
    VXY, EToV = uniform_mesh(Quad(), K1D)
    md = make_periodic(MeshData(VXY, EToV, rd))

    u0 = zeros(eltype(md.x), size(md.x, 1), size(md.x, 2), DG.NSTATE)
    exact_state!(u0, md.x, md.y, 0.0, mms)

    # The production RHS in some repository versions requires a source object.
    # This object is identically zero, so the manufactured forcing remains the
    # only forcing in the calculation.
    zero_spatial = zeros(eltype(md.x), size(md.x))
    zero_source = DG.CarcioneHeatSource(
        zero_spatial;
        amplitude=0.0,
        f0=1.0,
        t0=0.0,
        tau=material.tau,
        tmax=final_time,
        filter_dt=min(material.tau / 20, final_time / 100),
    )

    cache = DG.RHSCacheTE(u0, rd)
    parameters = (;
        rd,
        md,
        material,
        cache,
        source=zero_source,
        mms,
    )

    dt_stable, dtinfo = DG.estimate_dt(
        material,
        rd,
        md;
        CFL=cfl,
        relaxation_fraction=0.25,
    )

    problem = ODEProblem(rhs_mms!, u0, (0.0, final_time), parameters)

    local solution
    elapsed = if lowercase(time_scheme) == "vern9"
        @elapsed solution = solve(
            problem,
            Vern9();
            adaptive=true,
            abstol=tolerance,
            reltol=tolerance,
            dt=min(dt_stable / 4, final_time / 100),
            dtmax=dt_stable,
            saveat=[final_time],
            save_start=false,
            save_end=true,
            save_everystep=false,
            dense=false,
            progress=false,
        )
    elseif lowercase(time_scheme) == "ssprk54"
        dt_fixed = fixed_dt_fraction * dt_stable
        @elapsed solution = solve(
            problem,
            SSPRK54();
            adaptive=false,
            dt=dt_fixed,
            tstops=[final_time],
            saveat=[final_time],
            save_start=false,
            save_end=true,
            save_everystep=false,
            dense=false,
            progress=false,
        )
    else
        error("MMS_TIME_SCHEME must be 'vern9' or 'ssprk54'")
    end

    ufinal = solution.u[end]
    all(isfinite, ufinal) || error("non-finite solution for N=$N, K1D=$K1D")

    errors = mms_errors(ufinal, rd, md, material, mms, final_time)
    dofs = DG.NSTATE * size(u0, 1) * size(u0, 2)
    nsteps = hasproperty(solution, :destats) ? solution.destats.naccept : -1

    return (;
        N,
        K1D,
        h=dtinfo.h,
        dofs,
        dt_stable,
        nsteps,
        elapsed,
        errors,
        rd,
        md,
        ufinal,
    )
end

# -----------------------------------------------------------------------------
# Output writers
# -----------------------------------------------------------------------------
function write_csv(path, results)
    open(path, "w") do io
        header = [
            "N",
            "K1D",
            "h",
            "dofs",
            "dt_stability_cap",
            "accepted_steps",
            "wall_seconds",
            "combined_L2_absolute",
            "combined_L2_relative",
            "combined_L2_rate",
            "energy_absolute",
            "energy_relative",
            "energy_rate",
        ]

        for name in DG.STATE_FIELD_NAMES
            push!(header, "$(name)_L2_absolute")
            push!(header, "$(name)_L2_relative")
        end

        println(io, join(header, ','))

        for row in results
            values = Any[
                row.N,
                row.K1D,
                row.h,
                row.dofs,
                row.dt_stable,
                row.nsteps,
                row.elapsed,
                row.combined_abs,
                row.combined_rel,
                row.rate_l2,
                row.energy_abs,
                row.energy_rel,
                row.rate_energy,
            ]

            for field in 1:DG.NSTATE
                push!(values, row.component_abs[field])
                push!(values, row.component_rel[field])
            end

            println(io, join(values, ','))
        end
    end
end

function write_latex_table(path, results)
    open(path, "w") do io
        println(io, raw"\begin{table}[htbp]")
        println(io, raw"\centering")
        println(io, raw"\caption{Manufactured-solution convergence of the complete thermoelastic DG discretization.}")
        println(io, raw"\label{tab:mms-convergence}")
        println(io, raw"\begin{tabular}{rrrrrrrr}")
        println(io, raw"\toprule")
        println(io, raw"$N$ & $K_{1D}$ & $h$ & DOFs & rel. $L^2$ & rate & rel. energy & rate \\")
        println(io, raw"\midrule")

        for row in results
            l2rate = isfinite(row.rate_l2) ? @sprintf("%.2f", row.rate_l2) : "--"
            erate = isfinite(row.rate_energy) ? @sprintf("%.2f", row.rate_energy) : "--"
            @printf(
                io,
                "%d & %d & %.3e & %d & %.3e & %s & %.3e & %s \\\\\n",
                row.N,
                row.K1D,
                row.h,
                row.dofs,
                row.combined_rel,
                l2rate,
                row.energy_rel,
                erate,
            )
        end

        println(io, raw"\bottomrule")
        println(io, raw"\end{tabular}")
        println(io, raw"\end{table}")
    end
end

function make_plots(outdir, results, degrees)
    p_l2 = plot(
        xscale=:log10,
        yscale=:log10,
        xlabel="element width h",
        ylabel="relative combined L2 error",
        title="Manufactured-solution h-convergence",
        legend=:bottomright,
        linewidth=2,
    )

    p_energy = plot(
        xscale=:log10,
        yscale=:log10,
        xlabel="element width h",
        ylabel="relative energy-norm error",
        title="Energy-norm h-convergence",
        legend=:bottomright,
        linewidth=2,
    )

    p_order = plot(
        xscale=:log10,
        xlabel="element width h",
        ylabel="observed order",
        title="Observed energy-norm order",
        legend=:best,
        linewidth=2,
    )

    for N in degrees
        subset = sort(filter(row -> row.N == N, results); by=row -> row.h, rev=true)
        hs = getproperty.(subset, :h)
        l2errors = getproperty.(subset, :combined_rel)
        energyerrors = getproperty.(subset, :energy_rel)

        plot!(p_l2, hs, l2errors; marker=:circle, label="N=$N")
        plot!(p_energy, hs, energyerrors; marker=:circle, label="N=$N")

        if length(subset) > 1
            plot!(
                p_order,
                hs[2:end],
                getproperty.(subset[2:end], :rate_energy);
                marker=:circle,
                label="N=$N",
            )
            hline!(p_order, [N + 1]; linestyle=:dash, label="target $(N + 1)")
        end
    end

    savefig(p_l2, joinpath(outdir, "mms_h_convergence_L2.png"))
    savefig(p_energy, joinpath(outdir, "mms_h_convergence_energy.png"))
    savefig(p_order, joinpath(outdir, "mms_observed_orders.png"))
end

# -----------------------------------------------------------------------------
# Main refinement study
# -----------------------------------------------------------------------------
function main()
    degrees = parse_int_list("MMS_DEGREES", "1,2,3,4")
    meshes = parse_int_list("MMS_MESHES", "4,8,16,32")
    final_time = parse_env(Float64, "MMS_FINAL_TIME", 0.20)
    cfl = parse_env(Float64, "MMS_CFL", 0.05)
    penalty_scale = parse_env(Float64, "MMS_PENALTY_SCALE", 0.25)
    br1_penalty_scale = parse_env(Float64, "MMS_BR1_PENALTY_SCALE", 0.25)
    time_scheme = lowercase(get(ENV, "MMS_TIME_SCHEME", "vern9"))
    tolerance = parse_env(Float64, "MMS_TOL", 1.0e-13)
    fixed_dt_fraction = parse_env(Float64, "MMS_FIXED_DT_FRACTION", 0.25)

    final_time > 0 || error("MMS_FINAL_TIME must be positive")
    cfl > 0 || error("MMS_CFL must be positive")
    penalty_scale >= 0 || error("MMS_PENALTY_SCALE must be nonnegative")
    br1_penalty_scale >= 0 || error("MMS_BR1_PENALTY_SCALE must be nonnegative")
    tolerance > 0 || error("MMS_TOL must be positive")
    0 < fixed_dt_fraction <= 1 || error("MMS_FIXED_DT_FRACTION must be in (0,1]")
    time_scheme in ("vern9", "ssprk54") ||
        error("MMS_TIME_SCHEME must be 'vern9' or 'ssprk54'")

    outdir = make_output_dir()
    material = make_mms_material(
        penalty_scale=penalty_scale,
        br1_penalty_scale=br1_penalty_scale,
    )
    mms = MMSData()

    println("Manufactured-solution convergence study")
    println("Existing main.jl is not used or modified.")
    println("Output directory: ", outdir)
    println("Degrees: ", join(degrees, ", "))
    println("Meshes: ", join(meshes, ", "))
    @printf("Final time: %.6e\n", final_time)
    @printf("CFL: %.4f\n", cfl)
    @printf("Time scheme: %s\n", time_scheme)
    @printf("Tolerance: %.3e\n", tolerance)
    @printf("Mechanical penalty scale: %.4f\n", penalty_scale)
    @printf("BR1 temperature-jump scale: %.4f\n", br1_penalty_scale)
    @printf(
        "Maximum coupled characteristic speed: %.6e\n\n",
        DG.maximum_characteristic_speed(material),
    )

    results = NamedTuple[]

    println("  N   K1D          h        DOFs      rel L2   rate    rel energy   rate   seconds")
    println("----------------------------------------------------------------------------------")

    for N in degrees
        previous = nothing

        for K1D in meshes
            case = run_case(
                N,
                K1D,
                final_time,
                cfl,
                material,
                mms;
                time_scheme=time_scheme,
                tolerance=tolerance,
                fixed_dt_fraction=fixed_dt_fraction,
            )

            rate_l2 = previous === nothing ? NaN :
                log(previous.combined_abs / case.errors.combined_abs) /
                log(previous.h / case.h)

            rate_energy = previous === nothing ? NaN :
                log(previous.energy_abs / case.errors.energy_abs) /
                log(previous.h / case.h)

            row = (;
                N=case.N,
                K1D=case.K1D,
                h=case.h,
                dofs=case.dofs,
                dt_stable=case.dt_stable,
                nsteps=case.nsteps,
                elapsed=case.elapsed,
                combined_abs=case.errors.combined_abs,
                combined_rel=case.errors.combined_rel,
                rate_l2,
                energy_abs=case.errors.energy_abs,
                energy_rel=case.errors.energy_rel,
                rate_energy,
                component_abs=copy(case.errors.component_abs),
                component_rel=copy(case.errors.component_rel),
            )
            push!(results, row)

            l2rate = isfinite(rate_l2) ? @sprintf("%.2f", rate_l2) : "--"
            erate = isfinite(rate_energy) ? @sprintf("%.2f", rate_energy) : "--"

            @printf(
                "%3d %6d %10.3e %11d %11.3e %6s %13.3e %6s %9.2f\n",
                N,
                K1D,
                case.h,
                case.dofs,
                case.errors.combined_rel,
                l2rate,
                case.errors.energy_rel,
                erate,
                case.elapsed,
            )

            previous = (;
                h=case.h,
                combined_abs=case.errors.combined_abs,
                energy_abs=case.errors.energy_abs,
            )
        end

        println()
    end

    csv_path = joinpath(outdir, "mms_convergence.csv")
    tex_path = joinpath(outdir, "mms_convergence_table.tex")
    summary_path = joinpath(outdir, "mms_summary.txt")

    write_csv(csv_path, results)
    write_latex_table(tex_path, results)
    make_plots(outdir, results, degrees)

    open(summary_path, "w") do io
        println(io, "Manufactured-solution convergence study")
        println(io, "Domain: [-1,1]^2, periodic")
        println(io, "State: [s_xx,s_zz,s_xz,v_x,v_z,theta,q_x,q_z]")
        println(io, "Time scheme: ", time_scheme)
        @printf(io, "final time = %.16e\n", final_time)
        @printf(io, "CFL = %.16e\n", cfl)
        @printf(io, "tolerance = %.16e\n", tolerance)
        @printf(io, "mechanical penalty scale = %.16e\n", penalty_scale)
        @printf(io, "BR1 penalty scale = %.16e\n", br1_penalty_scale)
        println(io, "degrees = ", join(degrees, ','))
        println(io, "meshes = ", join(meshes, ','))
        println(io)
        println(io, "Expected asymptotic relative L2 order for a smooth solution: N+1.")
        println(io, "Errors are evaluated with an independent overintegrated Gauss rule.")
        println(io, "Observed rates are formed from absolute errors, not mesh-dependent relative denominators.")
        println(io, "Use the finest two or three meshes when reporting an observed rate.")
        println(io)

        for N in degrees
            subset = sort(filter(row -> row.N == N, results); by=row -> row.h, rev=true)
            lastrow = subset[end]
            @printf(
                io,
                "N=%d: finest rel L2=%.6e, rate=%.4f; rel energy=%.6e, rate=%.4f\n",
                N,
                lastrow.combined_rel,
                lastrow.rate_l2,
                lastrow.energy_rel,
                lastrow.rate_energy,
            )
        end
    end

    println()
    println("Saved CSV: ", csv_path)
    println("Saved LaTeX table: ", tex_path)
    println("Saved plots and summary in: ", outdir)
    println("Expected asymptotic spatial order: approximately N+1.")

    return nothing
end

end # module MMSConvergenceDriver


Base.invokelatest(MMSConvergenceDriver.main)

