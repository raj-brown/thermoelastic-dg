module DG

# Thermoelastic model and experiment:
#   J. M. Carcione et al., Geophysics 84 (2019), T1-T11,
#   DOI 10.1190/geo2018-0448.1.
using LinearAlgebra
using WriteVTK

export ThermoelasticMaterial
export CarcioneHeatSource
export RHSCacheTE
export make_material
export carcione_reference_material
export carcione_reference_properties
export carcione_pulse
export paper_heat_source
export first_order_heat_supply
export normalized_gaussian_source
export rhs_thermoelastic_br1!
export maximum_characteristic_speed
export estimate_dt
export compute_cfl
export energy_components
export energy_balance
export interpolate_state
export export_quad_subcells_vtu
export ManufacturedSolution
export manufactured_material
export manufactured_state!
export add_manufactured_forcing!
export rhs_thermoelastic_mms!
export manufactured_errors
export NSTATE, STATE_FIELD_NAMES
export ISXX, ISZZ, ISXZ, IVX, IVZ, ITHETA, IQX, IQZ

# -----------------------------------------------------------------------------
# State ordering
# -----------------------------------------------------------------------------
# The stored stresses are ELASTIC stresses s = C:epsilon.  The physical total
# thermoelastic stress is Sigma = s - b*theta and is reconstructed whenever the
# momentum equation, traction flux, or output routines need it.
const ISXX   = 1
const ISZZ   = 2
const ISXZ   = 3
const IVX    = 4
const IVZ    = 5
const ITHETA = 6
const IQX    = 7
const IQZ    = 8
const NSTATE = 8

const STATE_FIELD_NAMES = (
    "s_xx", "s_zz", "s_xz", "v_x", "v_z", "theta", "q_x", "q_z"
)

# -----------------------------------------------------------------------------
# Material data
# -----------------------------------------------------------------------------
struct ThermoelasticMaterial{T,MC<:AbstractMatrix{T},VB<:AbstractVector{T},MK<:AbstractMatrix{T}}
    rho::T
    C::MC
    S::MC
    b::VB
    heat_capacity::T
    Tref::T
    K::MK
    Kinv::MK
    tau::T
    alpha_v::T
    alpha_sigma::T
    eta_br1::T
end

"""
    _normal_matrix(nx, nz)

Return the two-dimensional engineering-shear normal operator

    A_n = [nx  0;
           0  nz;
           nz nx].
"""
function _normal_matrix(nx::T, nz::T) where {T}
    return T[nx 0; 0 nz; nz nx]
end

"""
    _maximum_mechanical_speed(C, rho; nsamples=128)

Sample the anisotropic acoustic tensor A_n' C A_n on the unit circle and
return a conservative maximum elastic speed.
"""
function _maximum_mechanical_speed(C::AbstractMatrix, rho::Real; nsamples::Int=128)
    cmax = 0.0
    for j in 0:(nsamples - 1)
        angle = 2pi * j / nsamples
        nx = cos(angle)
        nz = sin(angle)
        An = _normal_matrix(nx, nz)
        Qn = Symmetric(An' * C * An)
        cmax = max(cmax, sqrt(maximum(eigvals(Qn)) / rho))
    end
    return cmax
end

"""
    make_material(; rho, C, b, heat_capacity, Tref, K, tau,
                    penalty_scale=1, br1_penalty_scale=0,
                    alpha_v=nothing, alpha_sigma=nothing, eta_br1=nothing)

Construct a homogeneous two-dimensional thermoelastic material.  `C` is the
3-by-3 stiffness matrix in the engineering-shear Voigt ordering
`[xx, zz, xz]`, `b` is the three-component thermal-stress vector, and `K` is
the 2-by-2 conductivity tensor.

The Chan--Shukla flux uses scalar penalties.  In dimensional variables their
natural scales are a scalar impedance `rho*c_ref` and its reciprocal.  The
`penalty_scale` multiplies both.  `br1_penalty_scale=0` selects pure BR1.
"""
function make_material(;
    rho,
    C,
    b,
    heat_capacity,
    Tref,
    K,
    tau,
    penalty_scale=1.0,
    br1_penalty_scale=0.0,
    alpha_v=nothing,
    alpha_sigma=nothing,
    eta_br1=nothing,
)
    T = promote_type(
        Float64,
        typeof(rho),
        eltype(C),
        eltype(b),
        typeof(heat_capacity),
        typeof(Tref),
        eltype(K),
        typeof(tau),
    )

    rhoT = convert(T, rho)
    cT = convert(T, heat_capacity)
    TrefT = convert(T, Tref)
    tauT = convert(T, tau)
    Cmat = Matrix{T}(C)
    bvec = Vector{T}(b)
    Kmat = Matrix{T}(K)

    size(Cmat) == (3, 3) || throw(ArgumentError("C must be 3-by-3"))
    length(bvec) == 3 || throw(ArgumentError("b must have length 3"))
    size(Kmat) == (2, 2) || throw(ArgumentError("K must be 2-by-2"))
    rhoT > 0 || throw(ArgumentError("rho must be positive"))
    cT > 0 || throw(ArgumentError("heat_capacity must be positive"))
    TrefT > 0 || throw(ArgumentError("Tref must be positive"))
    tauT > 0 || throw(ArgumentError("tau must be positive"))
    isposdef(Symmetric(Cmat)) || throw(ArgumentError("C must be symmetric positive definite"))
    isposdef(Symmetric(Kmat)) || throw(ArgumentError("K must be symmetric positive definite"))

    Smat = inv(Cmat)
    Kinv = inv(Kmat)

    cref = _maximum_mechanical_speed(Cmat, rhoT)
    Zref = rhoT * cref

    alpha_vT = isnothing(alpha_v) ? convert(T, penalty_scale) * Zref : convert(T, alpha_v)
    alpha_sigmaT = isnothing(alpha_sigma) ? convert(T, penalty_scale) / Zref : convert(T, alpha_sigma)

    kmax = maximum(eigvals(Symmetric(Kmat)))
    Zthermal = sqrt(cT * kmax / tauT)
    etaT = isnothing(eta_br1) ? convert(T, br1_penalty_scale) * Zthermal : convert(T, eta_br1)

    alpha_vT >= 0 || throw(ArgumentError("alpha_v must be nonnegative"))
    alpha_sigmaT >= 0 || throw(ArgumentError("alpha_sigma must be nonnegative"))
    etaT >= 0 || throw(ArgumentError("eta_br1 must be nonnegative"))

    return ThermoelasticMaterial(
        rhoT,
        Cmat,
        Smat,
        bvec,
        cT,
        TrefT,
        Kmat,
        Kinv,
        tauT,
        alpha_vT,
        alpha_sigmaT,
        etaT,
    )
end

"""
    carcione_reference_material(; penalty_scale=1, br1_penalty_scale=0)

Return the homogeneous isotropic material used for the coupled heat-source
experiment in Figure 6 of Carcione et al. (2019).  The Lamé parameters are
computed from the reported isothermal P-wave and S-wave speeds.  The reported
thermal-stress coefficient beta = 7.92e4 Pa/K is used directly.
"""
function carcione_reference_material(;
    penalty_scale::Real=1.0,
    br1_penalty_scale::Real=0.0,
)
    rho = 2650.0
    vI = 2457.0
    vS = 1505.0
    heat_capacity = 117.0
    conductivity = 10.5
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

"""
    carcione_reference_properties(material)

Return characteristic reference values for the isotropic Carcione material.
For the supplied material these reproduce approximately
`V_A = 3480 m/s`, `V_Einf = 3980 m/s`, and `V_Tinf = 1517 m/s`.
"""
function carcione_reference_properties(material::ThermoelasticMaterial)
    rho = material.rho
    C11 = material.C[1, 1]
    C33 = material.C[3, 3]
    beta = material.b[1]
    c = material.heat_capacity
    Tref = material.Tref
    gamma = material.K[1, 1]
    tau = material.tau

    vI = sqrt(C11 / rho)
    vS = sqrt(C33 / rho)
    coupling_speed_sq = beta^2 * Tref / (rho * c)
    vA = sqrt(vI^2 + coupling_speed_sq)
    Minf = gamma / (c * tau)
    discriminant = max((vA^2 + Minf)^2 - 4Minf * vI^2, 0.0)
    vEinf = sqrt(0.5 * (vA^2 + Minf + sqrt(discriminant)))
    vTinf = sqrt(0.5 * (vA^2 + Minf - sqrt(discriminant)))

    return (;
        vI,
        vS,
        vA,
        vEinf,
        vTinf,
        beta,
        tau,
        thermal_diffusivity=gamma / c,
    )
end

# -----------------------------------------------------------------------------
# Carcione heat-source pulse and spatial regularization
# -----------------------------------------------------------------------------
"""
    carcione_pulse(t, f0, t0=3/(2*f0); mode=:printed)

Return the Gaussian-modulated source pulse.  Two conventions are retained
because equation (35) in Carcione et al. is printed without `2*pi`, while
`f0` is described in MHz:

* `mode=:printed` uses the equation exactly as typeset,
  `cos(f0*(t-t0))*exp(-2*(f0*(t-t0))^2)`.
* `mode=:cyclic_hz` interprets `f0` as cycles/s,
  `cos(2*pi*f0*(t-t0))*exp(-2*(f0*(t-t0))^2)`.

The second convention produces a pulse whose carrier is actually 3.5 MHz and
is useful when reproducing the appearance of Figure 7.  Both conventions are
available so the ambiguity is explicit rather than hidden.
"""
@inline function carcione_pulse(
    t,
    f0,
    t0=3 / (2 * f0);
    mode::Symbol=:printed,
)
    xi = (t - t0) * f0
    if mode === :printed
        return cos(xi) * exp(-2 * xi^2)
    elseif mode === :cyclic_hz
        return cospi(2 * xi) * exp(-2 * xi^2)
    else
        throw(ArgumentError(
            "unknown Carcione pulse mode $(repr(mode)); " *
            "use :printed or :cyclic_hz",
        ))
    end
end

struct CarcioneHeatSource{T,A<:AbstractMatrix{T},V<:AbstractVector{T}}
    spatial::A
    amplitude::T
    f0::T
    t0::T
    tau::T
    tmax::T
    filter_dt::T
    wavelet_mode::Symbol
    supply_history::V
end

"""
    CarcioneHeatSource(spatial; amplitude=1, f0=3.5e6,
                       t0=3/(2*f0), tau, tmax, filter_dt=nothing,
                       wavelet_mode=:printed)

Construct the centered heat source.  `wavelet_mode` is forwarded to
`carcione_pulse`; use `:printed` for the literal equation (35) or
`:cyclic_hz` when `f0` is interpreted as cycles per second.

`amplitude*carcione_pulse(t,f0,t0)` is Carcione's source `q(t)` in

    K*Delta(theta) = c*(theta_t + tau*theta_tt)
                     + Tref*b:E(v) + tau*Tref*b:E(v_t) + q.

The energy-stable first-order Cattaneo formulation uses a heat-supply rate
`r(t)` in

    c*theta_t + Tref*b:E(v) + div(q_flux) = r.

Exact elimination of the physical heat flux gives the source combination
`r + tau*r_t`.  Therefore, to reproduce Carcione's second-order forcing,
this constructor precomputes the causal scalar filter

    tau*r_t + r = -q(t),    r(0)=0.

The filter is integrated with RK4 on a uniform subgrid.  The default filter
step resolves both the relaxation time and the source pulse.
"""
function CarcioneHeatSource(
    spatial::AbstractMatrix;
    amplitude::Real=1.0,
    f0::Real=3.5e6,
    t0::Real=3 / (2 * f0),
    tau::Real,
    tmax::Real,
    filter_dt=nothing,
    wavelet_mode::Symbol=:printed,
)
    wavelet_mode in (:printed, :cyclic_hz) || throw(ArgumentError(
        "wavelet_mode must be :printed or :cyclic_hz",
    ))

    T = promote_type(
        Float64,
        eltype(spatial),
        typeof(amplitude),
        typeof(f0),
        typeof(t0),
        typeof(tau),
        typeof(tmax),
    )

    amplitudeT = convert(T, amplitude)
    f0T = convert(T, f0)
    t0T = convert(T, t0)
    tauT = convert(T, tau)
    tmaxT = convert(T, tmax)

    f0T > 0 || throw(ArgumentError("f0 must be positive"))
    tauT > 0 || throw(ArgumentError("tau must be positive"))
    tmaxT > 0 || throw(ArgumentError("tmax must be positive"))

    requested_dt = if isnothing(filter_dt)
        min(tauT / 20, inv(200 * f0T))
    else
        convert(T, filter_dt)
    end
    requested_dt > 0 || throw(ArgumentError("filter_dt must be positive"))

    nsteps = max(1, ceil(Int, tmaxT / requested_dt))
    dtT = tmaxT / nsteps
    history = zeros(T, nsteps + 1)

    qpaper(t) = amplitudeT * carcione_pulse(
        t,
        f0T,
        t0T;
        mode=wavelet_mode,
    )
    source_rhs(t, r) = (-qpaper(t) - r) / tauT

    @inbounds for n in 1:nsteps
        tn = (n - 1) * dtT
        rn = history[n]
        k1 = source_rhs(tn, rn)
        k2 = source_rhs(tn + dtT / 2, rn + dtT * k1 / 2)
        k3 = source_rhs(tn + dtT / 2, rn + dtT * k2 / 2)
        k4 = source_rhs(tn + dtT, rn + dtT * k3)
        history[n + 1] = rn + (dtT / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    end

    return CarcioneHeatSource(
        Matrix{T}(spatial),
        amplitudeT,
        f0T,
        t0T,
        tauT,
        tmaxT,
        dtT,
        wavelet_mode,
        history,
    )
end

"""Return Carcione's second-order heat-source history `q(t)`."""
@inline paper_heat_source(source::CarcioneHeatSource, t) =
    source.amplitude * carcione_pulse(
        t,
        source.f0,
        source.t0;
        mode=source.wavelet_mode,
    )

"""
    first_order_heat_supply(source, t)

Return the linearly interpolated first-order supply `r(t)` satisfying
`tau*r_t + r = -q(t)`.  The precomputed interval is `[0, source.tmax]`.
"""
@inline function first_order_heat_supply(source::CarcioneHeatSource, t)
    t <= zero(t) && return source.supply_history[1]
    t >= source.tmax && return source.supply_history[end]

    coordinate = t / source.filter_dt
    left = clamp(floor(Int, coordinate) + 1, 1, length(source.supply_history) - 1)
    fraction = coordinate - (left - 1)
    return muladd(
        fraction,
        source.supply_history[left + 1] - source.supply_history[left],
        source.supply_history[left],
    )
end

"""
    normalized_gaussian_source(rd, md; x0=0, z0=0, sigma)

Construct a smooth approximation to a two-dimensional point source and
normalize its physical-domain integral to one using StartUpDG volume
quadrature.  This avoids the four-copy ambiguity of placing a DG point source
at a shared element vertex.
"""
function normalized_gaussian_source(rd, md; x0::Real=0.0, z0::Real=0.0, sigma::Real)
    sigma > 0 || throw(ArgumentError("sigma must be positive"))
    (; x, y, wJq) = md
    (; Vq) = rd

    source = @. exp(-((x - x0)^2 + (y - z0)^2) / (2sigma^2))
    source_q = Vq * source
    integral = sum(wJq .* source_q)

    isfinite(integral) && integral > 0 ||
        throw(ArgumentError("source normalization failed; integral = $integral"))

    source ./= integral
    return source
end

# -----------------------------------------------------------------------------
# Cache
# -----------------------------------------------------------------------------
# Volume scratch slots
const V_SIGXX       = 1
const V_SIGZZ       = 2
const V_SIGXZ       = 3
const V_VX_XJ       = 4
const V_VX_ZJ       = 5
const V_VZ_XJ       = 6
const V_VZ_ZJ       = 7
const V_EPSXX       = 8
const V_EPSZZ       = 9
const V_GAMMAXZ     = 10
const V_DR          = 11
const V_DS          = 12
const V_DXJ         = 13
const V_DZJ         = 14
const V_DIVSIGXJ    = 15
const V_DIVSIGZJ    = 16
const V_GRADTHETAX  = 17
const V_GRADTHETAZ  = 18
const V_DIVQ         = 19
const V_LIFTVOL     = 20
const NVOLUME_CACHE = 20

# Face scratch slots
const F_SIGXXM       = 1
const F_SIGXXP       = 2
const F_SIGZZM       = 3
const F_SIGZZP       = 4
const F_SIGXZM       = 5
const F_SIGXZP       = 6
const F_VXM          = 7
const F_VXP          = 8
const F_VZM          = 9
const F_VZP          = 10
const F_THETAM       = 11
const F_THETAP       = 12
const F_QXM          = 13
const F_QXP          = 14
const F_QZM          = 15
const F_QZP          = 16
const F_DSIGXX       = 17
const F_DSIGZZ       = 18
const F_DSIGXZ       = 19
const F_DVX          = 20
const F_DVZ          = 21
const F_DTHETA       = 22
const F_DQX          = 23
const F_DQZ          = 24
const F_DTX          = 25
const F_DTZ          = 26
const F_DVHATX       = 27
const F_DVHATZ       = 28
const F_STRESSCORRXX = 29
const F_STRESSCORRZZ = 30
const F_STRESSCORRXZ = 31
const F_TRACTIONCORRX = 32
const F_TRACTIONCORRZ = 33
const F_QFLUXCORR    = 34
const F_THETACORRX   = 35
const F_THETACORRZ   = 36
const F_LIFTBUF      = 37
const NFACE_CACHE    = 37

struct RHSCacheTE{T,AV<:Array{T,3},AF<:Array{T,3}}
    volume::AV
    face::AF
end

function RHSCacheTE(u::AbstractArray, rd)
    T = eltype(u)
    Np, K = size(u, 1), size(u, 2)
    Nf = size(rd.Vf, 1)
    volume = zeros(T, Np, K, NVOLUME_CACHE)
    face = zeros(T, Nf, K, NFACE_CACHE)
    return RHSCacheTE(volume, face)
end

@inline function _gather_plus!(up, um, mapP)
    @inbounds for i in eachindex(up)
        up[i] = um[mapP[i]]
    end
    return nothing
end

@inline function _derivative_numerators!(dxJ, dzJ, dr, ds, field, Dr, Ds, rxJ, sxJ, ryJ, syJ)
    mul!(dr, Dr, field)
    mul!(ds, Ds, field)
    @. dxJ = rxJ * dr + sxJ * ds
    @. dzJ = ryJ * dr + syJ * ds
    return nothing
end

# -----------------------------------------------------------------------------
# Semi-discrete weak--strong DG operator
# -----------------------------------------------------------------------------
"""
    rhs_thermoelastic_br1!(du, u, parameters, time)

Eight-field two-dimensional thermoelastic DG operator:

    u = [s_xx, s_zz, s_xz, v_x, v_z, theta, q_x, q_z].

The elastic stress equation is the weak equation evaluated through its
strong-equivalent SBP form.  Velocity is in strong form.  The mechanical
numerical flux is the scalar Chan--Shukla penalty flux, with total stress
`Sigma = s - b*theta`.  Temperature is weak and heat flux is strong, using
BR1 traces.  `eta_br1 = 0` gives pure BR1.

`parameters` must provide `rd`, `md`, `material`, and `cache`.  A physical
heat source is optional; omit `source` or set it to `nothing` for a source-free
problem.  The mesh is assumed periodic in this experiment, so every face has
a valid neighbor in `mapP`.
"""
function rhs_thermoelastic_br1!(du, u, parameters, time)
    (; rd, md, material, cache) = parameters
    source = hasproperty(parameters, :source) ? getproperty(parameters, :source) : nothing
    (; Vf, Dr, Ds, LIFT) = rd
    (; rxJ, sxJ, ryJ, syJ, nx, ny, J, Jf, mapP) = md

    m = material
    c = cache
    fill!(du, zero(eltype(du)))

    C11 = m.C[1, 1]
    C12 = m.C[1, 2]
    C13 = m.C[1, 3]
    C22 = m.C[2, 2]
    C23 = m.C[2, 3]
    C33 = m.C[3, 3]

    b1, b2, b3 = m.b
    K11 = m.K[1, 1]
    K12 = m.K[1, 2]
    K21 = m.K[2, 1]
    K22 = m.K[2, 2]

    @views begin
        # State views
        sxx = u[:, :, ISXX]
        szz = u[:, :, ISZZ]
        sxz = u[:, :, ISXZ]
        vx = u[:, :, IVX]
        vz = u[:, :, IVZ]
        theta = u[:, :, ITHETA]
        qx = u[:, :, IQX]
        qz = u[:, :, IQZ]

        dsxx_dt = du[:, :, ISXX]
        dszz_dt = du[:, :, ISZZ]
        dsxz_dt = du[:, :, ISXZ]
        dvx_dt = du[:, :, IVX]
        dvz_dt = du[:, :, IVZ]
        dtheta_dt = du[:, :, ITHETA]
        dqx_dt = du[:, :, IQX]
        dqz_dt = du[:, :, IQZ]

        # Volume cache views
        sigxx = c.volume[:, :, V_SIGXX]
        sigzz = c.volume[:, :, V_SIGZZ]
        sigxz = c.volume[:, :, V_SIGXZ]
        vx_xJ = c.volume[:, :, V_VX_XJ]
        vx_zJ = c.volume[:, :, V_VX_ZJ]
        vz_xJ = c.volume[:, :, V_VZ_XJ]
        vz_zJ = c.volume[:, :, V_VZ_ZJ]
        epsxx = c.volume[:, :, V_EPSXX]
        epszz = c.volume[:, :, V_EPSZZ]
        gammaxz = c.volume[:, :, V_GAMMAXZ]
        dr = c.volume[:, :, V_DR]
        ds = c.volume[:, :, V_DS]
        dxJ = c.volume[:, :, V_DXJ]
        dzJ = c.volume[:, :, V_DZJ]
        divsigxJ = c.volume[:, :, V_DIVSIGXJ]
        divsigzJ = c.volume[:, :, V_DIVSIGZJ]
        gradtheta_x = c.volume[:, :, V_GRADTHETAX]
        gradtheta_z = c.volume[:, :, V_GRADTHETAZ]
        divq = c.volume[:, :, V_DIVQ]
        liftvol = c.volume[:, :, V_LIFTVOL]

        # Face cache views
        sigxxm = c.face[:, :, F_SIGXXM]
        sigxxp = c.face[:, :, F_SIGXXP]
        sigzzm = c.face[:, :, F_SIGZZM]
        sigzzp = c.face[:, :, F_SIGZZP]
        sigxzm = c.face[:, :, F_SIGXZM]
        sigxzp = c.face[:, :, F_SIGXZP]
        vxm = c.face[:, :, F_VXM]
        vxp = c.face[:, :, F_VXP]
        vzm = c.face[:, :, F_VZM]
        vzp = c.face[:, :, F_VZP]
        thetam = c.face[:, :, F_THETAM]
        thetap = c.face[:, :, F_THETAP]
        qxm = c.face[:, :, F_QXM]
        qxp = c.face[:, :, F_QXP]
        qzm = c.face[:, :, F_QZM]
        qzp = c.face[:, :, F_QZP]
        dsigxx = c.face[:, :, F_DSIGXX]
        dsigzz = c.face[:, :, F_DSIGZZ]
        dsigxz = c.face[:, :, F_DSIGXZ]
        dvx = c.face[:, :, F_DVX]
        dvz = c.face[:, :, F_DVZ]
        dtheta = c.face[:, :, F_DTHETA]
        dqx = c.face[:, :, F_DQX]
        dqz = c.face[:, :, F_DQZ]
        dtx = c.face[:, :, F_DTX]
        dtz = c.face[:, :, F_DTZ]
        dvhatx = c.face[:, :, F_DVHATX]
        dvhatz = c.face[:, :, F_DVHATZ]
        stresscorrxx = c.face[:, :, F_STRESSCORRXX]
        stresscorrzz = c.face[:, :, F_STRESSCORRZZ]
        stresscorrxz = c.face[:, :, F_STRESSCORRXZ]
        tractioncorrx = c.face[:, :, F_TRACTIONCORRX]
        tractioncorrz = c.face[:, :, F_TRACTIONCORRZ]
        qfluxcorr = c.face[:, :, F_QFLUXCORR]
        thetacorrx = c.face[:, :, F_THETACORRX]
        thetacorrz = c.face[:, :, F_THETACORRZ]
        liftbuf = c.face[:, :, F_LIFTBUF]

        # -----------------------------------------------------------------
        # Total stress Sigma = s - b*theta
        # -----------------------------------------------------------------
        @. sigxx = sxx - b1 * theta
        @. sigzz = szz - b2 * theta
        @. sigxz = sxz - b3 * theta

        # -----------------------------------------------------------------
        # Face traces and plus states
        # -----------------------------------------------------------------
        mul!(sigxxm, Vf, sigxx)
        mul!(sigzzm, Vf, sigzz)
        mul!(sigxzm, Vf, sigxz)
        mul!(vxm, Vf, vx)
        mul!(vzm, Vf, vz)
        mul!(thetam, Vf, theta)
        mul!(qxm, Vf, qx)
        mul!(qzm, Vf, qz)

        _gather_plus!(sigxxp, sigxxm, mapP)
        _gather_plus!(sigzzp, sigzzm, mapP)
        _gather_plus!(sigxzp, sigxzm, mapP)
        _gather_plus!(vxp, vxm, mapP)
        _gather_plus!(vzp, vzm, mapP)
        _gather_plus!(thetap, thetam, mapP)
        _gather_plus!(qxp, qxm, mapP)
        _gather_plus!(qzp, qzm, mapP)

        # Jump convention: plus minus minus, matching the theory.
        @. dsigxx = sigxxp - sigxxm
        @. dsigzz = sigzzp - sigzzm
        @. dsigxz = sigxzp - sigxzm
        @. dvx = vxp - vxm
        @. dvz = vzp - vzm
        @. dtheta = thetap - thetam
        @. dqx = qxp - qxm
        @. dqz = qzp - qzm

        # A_n' jump(Sigma): total-traction jump using the local outward normal.
        @. dtx = nx * dsigxx + ny * dsigxz
        @. dtz = nx * dsigxz + ny * dsigzz

        # Chan--Shukla v-hat correction:
        # vhat - vminus = 1/2 jump(v) + alpha_sigma/2 A_n' jump(Sigma).
        @. dvhatx = 0.5 * dvx + 0.5 * m.alpha_sigma * dtx
        @. dvhatz = 0.5 * dvz + 0.5 * m.alpha_sigma * dtz

        # -----------------------------------------------------------------
        # Weak elastic-stress equation, evaluated in strong-equivalent form
        # S s_t = E(v) + LIFT(A_n(vhat-vminus)).
        # -----------------------------------------------------------------
        _derivative_numerators!(vx_xJ, vx_zJ, dr, ds, vx, Dr, Ds, rxJ, sxJ, ryJ, syJ)
        _derivative_numerators!(vz_xJ, vz_zJ, dr, ds, vz, Dr, Ds, rxJ, sxJ, ryJ, syJ)

        @. stresscorrxx = nx * dvhatx
        @. stresscorrzz = ny * dvhatz
        @. stresscorrxz = ny * dvhatx + nx * dvhatz

        @. liftbuf = stresscorrxx * Jf
        mul!(liftvol, LIFT, liftbuf)
        @. epsxx = (vx_xJ + liftvol) / J

        @. liftbuf = stresscorrzz * Jf
        mul!(liftvol, LIFT, liftbuf)
        @. epszz = (vz_zJ + liftvol) / J

        @. liftbuf = stresscorrxz * Jf
        mul!(liftvol, LIFT, liftbuf)
        @. gammaxz = (vx_zJ + vz_xJ + liftvol) / J

        @. dsxx_dt = C11 * epsxx + C12 * epszz + C13 * gammaxz
        @. dszz_dt = C12 * epsxx + C22 * epszz + C23 * gammaxz
        @. dsxz_dt = C13 * epsxx + C23 * epszz + C33 * gammaxz

        # -----------------------------------------------------------------
        # Strong velocity equation
        # rho v_t = D(Sigma) + LIFT(that - tminus).
        # -----------------------------------------------------------------
        _derivative_numerators!(dxJ, dzJ, dr, ds, sigxx, Dr, Ds, rxJ, sxJ, ryJ, syJ)
        divsigxJ .= dxJ

        _derivative_numerators!(dxJ, dzJ, dr, ds, sigxz, Dr, Ds, rxJ, sxJ, ryJ, syJ)
        @. divsigxJ += dzJ
        divsigzJ .= dxJ

        _derivative_numerators!(dxJ, dzJ, dr, ds, sigzz, Dr, Ds, rxJ, sxJ, ryJ, syJ)
        @. divsigzJ += dzJ

        # A_n' A_n jump(v), written without forming A_n.
        @. tractioncorrx =
            0.5 * dtx +
            0.5 * m.alpha_v * (
                nx * (nx * dvx) +
                ny * (ny * dvx + nx * dvz)
            )

        @. tractioncorrz =
            0.5 * dtz +
            0.5 * m.alpha_v * (
                ny * (ny * dvz) +
                nx * (ny * dvx + nx * dvz)
            )

        @. liftbuf = tractioncorrx * Jf
        mul!(liftvol, LIFT, liftbuf)
        @. dvx_dt = (divsigxJ + liftvol) / (m.rho * J)

        @. liftbuf = tractioncorrz * Jf
        mul!(liftvol, LIFT, liftbuf)
        @. dvz_dt = (divsigzJ + liftvol) / (m.rho * J)

        # -----------------------------------------------------------------
        # Weak temperature equation
        # c theta_t = r - Tref*b'E_h(v) - div_h(q).
        # Using the same epsxx/epszz/gammaxz arrays enforces the exact
        # semi-discrete thermoelastic cancellation proved in the manuscript.
        # -----------------------------------------------------------------
        _derivative_numerators!(dxJ, dzJ, dr, ds, qx, Dr, Ds, rxJ, sxJ, ryJ, syJ)
        divq .= dxJ
        _derivative_numerators!(dxJ, dzJ, dr, ds, qz, Dr, Ds, rxJ, sxJ, ryJ, syJ)
        @. divq += dzJ

        # qhat_n - qminus_n = 1/2 jump(q).n - eta_br1 jump(theta).
        @. qfluxcorr = 0.5 * (nx * dqx + ny * dqz) - m.eta_br1 * dtheta
        @. liftbuf = qfluxcorr * Jf
        mul!(liftvol, LIFT, liftbuf)
        @. divq = (divq + liftvol) / J

        if source === nothing
            @. dtheta_dt = (
                -m.Tref * (b1 * epsxx + b2 * epszz + b3 * gammaxz) -
                divq
            ) / m.heat_capacity
        else
            source_supply = first_order_heat_supply(source, time)
            @. dtheta_dt = (
                source_supply * source.spatial -
                m.Tref * (b1 * epsxx + b2 * epszz + b3 * gammaxz) -
                divq
            ) / m.heat_capacity
        end

        # -----------------------------------------------------------------
        # Strong Cattaneo heat-flux equation with the BR1 theta trace
        # tau q_t + q + K grad_h(theta) = 0.
        # -----------------------------------------------------------------
        _derivative_numerators!(dxJ, dzJ, dr, ds, theta, Dr, Ds, rxJ, sxJ, ryJ, syJ)

        @. thetacorrx = 0.5 * dtheta * nx
        @. liftbuf = thetacorrx * Jf
        mul!(liftvol, LIFT, liftbuf)
        @. gradtheta_x = (dxJ + liftvol) / J

        @. thetacorrz = 0.5 * dtheta * ny
        @. liftbuf = thetacorrz * Jf
        mul!(liftvol, LIFT, liftbuf)
        @. gradtheta_z = (dzJ + liftvol) / J

        @. dqx_dt = -(qx + K11 * gradtheta_x + K12 * gradtheta_z) / m.tau
        @. dqz_dt = -(qz + K21 * gradtheta_x + K22 * gradtheta_z) / m.tau
    end

    return nothing
end

# -----------------------------------------------------------------------------
# Periodic manufactured solution and convergence diagnostics
# -----------------------------------------------------------------------------

"""
    ManufacturedSolution(; kx=pi, kz=pi, omega=1.3, phase_x=0.23, phase_z=-0.31, amplitudes=(...))

Parameters for a smooth periodic manufactured solution on `[-1,1]^2`.
The eight amplitudes correspond to

    [s_xx, s_zz, s_xz, v_x, v_z, theta, q_x, q_z].

The trigonometric fields are deliberately independent, so the manufactured
forcing exercises every mechanical, thermoelastic, conductivity, and
Cattaneo-relaxation term in `rhs_thermoelastic_br1!`.
"""
struct ManufacturedSolution{T}
    kx::T
    kz::T
    omega::T
    phase_x::T
    phase_z::T
    amplitudes::NTuple{NSTATE,T}
end

function ManufacturedSolution(;
    kx::Real=pi,
    kz::Real=pi,
    omega::Real=1.3,
    phase_x::Real=0.23,
    phase_z::Real=-0.31,
    amplitudes=(0.8, -0.6, 0.5, 0.7, -0.9, 0.4, 0.3, -0.2),
)
    length(amplitudes) == NSTATE ||
        throw(ArgumentError("amplitudes must contain $NSTATE entries"))

    T = promote_type(
        Float64,
        typeof(kx),
        typeof(kz),
        typeof(omega),
        typeof(phase_x),
        typeof(phase_z),
        map(typeof, amplitudes)...,
    )
    amps = ntuple(i -> convert(T, amplitudes[i]), NSTATE)
    return ManufacturedSolution(
        convert(T, kx),
        convert(T, kz),
        convert(T, omega),
        convert(T, phase_x),
        convert(T, phase_z),
        amps,
    )
end

"""
    manufactured_material(; penalty_scale=1, br1_penalty_scale=0.25)

Return an O(1), homogeneous anisotropic material for the manufactured-solution
study.  Both the stiffness matrix and conductivity tensor are symmetric
positive definite, and the thermal-stress vector has a nonzero shear entry.
This makes the test exercise the anisotropic code path without the very small
physical relaxation time of the Carcione benchmark.
"""
function manufactured_material(;
    penalty_scale::Real=1.0,
    br1_penalty_scale::Real=0.25,
)
    C = [
        4.0   1.0    0.25
        1.0   3.0   -0.15
        0.25 -0.15   1.20
    ]
    b = [0.30, -0.20, 0.10]
    K = [0.80 0.12; 0.12 0.60]

    return make_material(
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

"""
    manufactured_state!(u, x, z, time, mms)

Fill `u` with the exact periodic manufactured state.  The coordinates may be
solution nodes, quadrature nodes, or plotting nodes, provided `x` and `z` have
the same two-dimensional shape as the first two dimensions of `u`.
"""
function manufactured_state!(
    u::AbstractArray,
    x::AbstractMatrix,
    z::AbstractMatrix,
    time,
    mms::ManufacturedSolution,
)
    size(x) == size(z) || throw(DimensionMismatch("x and z must have the same size"))
    size(u, 1) == size(x, 1) || throw(DimensionMismatch("u and x have different point counts"))
    size(u, 2) == size(x, 2) || throw(DimensionMismatch("u and x have different element counts"))
    size(u, 3) == NSTATE || throw(DimensionMismatch("u must contain $NSTATE fields"))

    a1, a2, a3, a4, a5, a6, a7, a8 = mms.amplitudes
    st, ct = sincos(mms.omega * time)

    @inbounds for k in axes(x, 2), i in axes(x, 1)
        sx, cx = sincos(mms.kx * x[i, k] + mms.phase_x)
        sz, cz = sincos(mms.kz * z[i, k] + mms.phase_z)

        u[i, k, ISXX]   =  a1 * sx * cz * ct
        u[i, k, ISZZ]   =  a2 * cx * sz * st
        u[i, k, ISXZ]   =  a3 * sx * sz * ct
        u[i, k, IVX]    =  a4 * cx * cz * st
        u[i, k, IVZ]    =  a5 * sx * cz * ct
        u[i, k, ITHETA] =  a6 * cx * sz * ct
        u[i, k, IQX]    =  a7 * sx * sz * st
        u[i, k, IQZ]    =  a8 * cx * cz * ct
    end

    return u
end

"""
    add_manufactured_forcing!(du, x, z, time, material, mms)

Add the exact source terms needed to make `manufactured_state!` solve the
complete first-order thermoelastic-Cattaneo system.  The additions are made
in state-rate form:

    s_t     = C E(v)                         + f_s,
    v_t     = rho^(-1) div(s - b theta)      + f_v,
    theta_t = -(Tref b' E(v) + div(q))/c     + f_theta,
    q_t     = -(q + K grad(theta))/tau       + f_q.

Consequently, this routine can be added directly after the homogeneous DG
operator, without changing the numerical fluxes.
"""
function add_manufactured_forcing!(
    du::AbstractArray,
    x::AbstractMatrix,
    z::AbstractMatrix,
    time,
    material::ThermoelasticMaterial,
    mms::ManufacturedSolution,
)
    size(du, 1) == size(x, 1) || throw(DimensionMismatch("du and x have different point counts"))
    size(du, 2) == size(x, 2) || throw(DimensionMismatch("du and x have different element counts"))
    size(du, 3) == NSTATE || throw(DimensionMismatch("du must contain $NSTATE fields"))

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

        # Exact state values needed by the Cattaneo forcing.
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

        # Derivatives of elastic stress and temperature.
        sxx_x =  a1 * kx * cx * cz * ct
        szz_z =  a2 * kz * cx * cz * st
        sxz_x =  a3 * kx * cx * sz * ct
        sxz_z =  a3 * kz * sx * cz * ct
        theta_x = -a6 * kx * sx * sz * ct
        theta_z =  a6 * kz * cx * cz * ct

        # Total stress Sigma = s - b theta.
        divsig_x = (sxx_x - b1 * theta_x) + (sxz_z - b3 * theta_z)
        divsig_z = (sxz_x - b3 * theta_x) + (szz_z - b2 * theta_z)

        # Divergence of the exact heat flux.
        qx_x =  a7 * kx * cx * sz * st
        qz_z = -a8 * kz * cx * sz * ct
        divq = qx_x + qz_z

        b_dot_eps = b1 * epsxx + b2 * epszz + b3 * gammaxz

        # Exact constitutive rate C E(v).
        constitutive_xx = C11 * epsxx + C12 * epszz + C13 * gammaxz
        constitutive_zz = C12 * epsxx + C22 * epszz + C23 * gammaxz
        constitutive_xz = C13 * epsxx + C23 * epszz + C33 * gammaxz

        # Add state-rate forcing.
        du[i, k, ISXX] += sxx_t - constitutive_xx
        du[i, k, ISZZ] += szz_t - constitutive_zz
        du[i, k, ISXZ] += sxz_t - constitutive_xz

        du[i, k, IVX] += vx_t - divsig_x / rho
        du[i, k, IVZ] += vz_t - divsig_z / rho

        du[i, k, ITHETA] +=
            theta_t + (Tref * b_dot_eps + divq) / heat_capacity

        du[i, k, IQX] +=
            qx_t + (qx + K11 * theta_x + K12 * theta_z) / tau
        du[i, k, IQZ] +=
            qz_t + (qz + K21 * theta_x + K22 * theta_z) / tau
    end

    return nothing
end

"""
    rhs_thermoelastic_mms!(du, u, parameters, time)

Manufactured-solution wrapper around `rhs_thermoelastic_br1!`.  The parameter
container must provide `mms` in addition to the standard DG parameters.  No
physical heat source is required; omit `source` or set it to `nothing`.
"""
function rhs_thermoelastic_mms!(du, u, parameters, time)
    rhs_thermoelastic_br1!(du, u, parameters, time)
    add_manufactured_forcing!(
        du,
        parameters.md.x,
        parameters.md.y,
        time,
        parameters.material,
        parameters.mms,
    )
    return nothing
end

"""
    manufactured_errors(u, rd, md, material, mms, time)

Compute componentwise absolute and relative L2 errors, an unweighted combined
L2 error, and the physically weighted thermoelastic energy-norm error.  All
integrals use the volume quadrature stored in `rd` and `md`.
"""
function manufactured_errors(
    u::AbstractArray,
    rd,
    md,
    material::ThermoelasticMaterial,
    mms::ManufacturedSolution,
    time,
)
    (; Vq) = rd
    (; wJq) = md

    nq = size(Vq, 1)
    nelements = size(u, 2)
    T = promote_type(eltype(u), eltype(wJq), typeof(time))

    numerical = Array{T}(undef, nq, nelements, NSTATE)
    exact = similar(numerical)

    @views for field in 1:NSTATE
        mul!(numerical[:, :, field], Vq, u[:, :, field])
    end

    xq = Vq * md.x
    zq = Vq * md.y
    manufactured_state!(exact, xq, zq, time, mms)

    component_absolute = zeros(T, NSTATE)
    component_exact = zeros(T, NSTATE)

    @inbounds for field in 1:NSTATE
        error_squared = zero(T)
        exact_squared = zero(T)
        for k in axes(wJq, 2), i in axes(wJq, 1)
            weight = wJq[i, k]
            difference = numerical[i, k, field] - exact[i, k, field]
            error_squared += weight * difference^2
            exact_squared += weight * exact[i, k, field]^2
        end
        component_absolute[field] = sqrt(max(error_squared, zero(T)))
        component_exact[field] = sqrt(max(exact_squared, zero(T)))
    end

    component_relative = similar(component_absolute)
    @inbounds for field in 1:NSTATE
        component_relative[field] = component_absolute[field] /
            max(component_exact[field], eps(T))
    end

    combined_absolute = norm(component_absolute)
    combined_exact = norm(component_exact)
    combined_relative = combined_absolute / max(combined_exact, eps(T))

    S = material.S
    Kinv = material.Kinv
    energy_error_squared = zero(T)
    energy_exact_squared = zero(T)

    @inbounds for k in axes(wJq, 2), i in axes(wJq, 1)
        weight = wJq[i, k]

        es1 = numerical[i, k, ISXX] - exact[i, k, ISXX]
        es2 = numerical[i, k, ISZZ] - exact[i, k, ISZZ]
        es3 = numerical[i, k, ISXZ] - exact[i, k, ISXZ]
        ev1 = numerical[i, k, IVX] - exact[i, k, IVX]
        ev2 = numerical[i, k, IVZ] - exact[i, k, IVZ]
        etheta = numerical[i, k, ITHETA] - exact[i, k, ITHETA]
        eq1 = numerical[i, k, IQX] - exact[i, k, IQX]
        eq2 = numerical[i, k, IQZ] - exact[i, k, IQZ]

        stress_error =
            es1 * (S[1, 1] * es1 + S[1, 2] * es2 + S[1, 3] * es3) +
            es2 * (S[2, 1] * es1 + S[2, 2] * es2 + S[2, 3] * es3) +
            es3 * (S[3, 1] * es1 + S[3, 2] * es2 + S[3, 3] * es3)
        flux_error =
            eq1 * (Kinv[1, 1] * eq1 + Kinv[1, 2] * eq2) +
            eq2 * (Kinv[2, 1] * eq1 + Kinv[2, 2] * eq2)

        energy_error_density =
            stress_error +
            material.rho * (ev1^2 + ev2^2) +
            material.heat_capacity / material.Tref * etheta^2 +
            material.tau / material.Tref * flux_error

        xs1 = exact[i, k, ISXX]
        xs2 = exact[i, k, ISZZ]
        xs3 = exact[i, k, ISXZ]
        xv1 = exact[i, k, IVX]
        xv2 = exact[i, k, IVZ]
        xtheta = exact[i, k, ITHETA]
        xq1 = exact[i, k, IQX]
        xq2 = exact[i, k, IQZ]

        stress_exact =
            xs1 * (S[1, 1] * xs1 + S[1, 2] * xs2 + S[1, 3] * xs3) +
            xs2 * (S[2, 1] * xs1 + S[2, 2] * xs2 + S[2, 3] * xs3) +
            xs3 * (S[3, 1] * xs1 + S[3, 2] * xs2 + S[3, 3] * xs3)
        flux_exact =
            xq1 * (Kinv[1, 1] * xq1 + Kinv[1, 2] * xq2) +
            xq2 * (Kinv[2, 1] * xq1 + Kinv[2, 2] * xq2)

        energy_exact_density =
            stress_exact +
            material.rho * (xv1^2 + xv2^2) +
            material.heat_capacity / material.Tref * xtheta^2 +
            material.tau / material.Tref * flux_exact

        energy_error_squared += weight * energy_error_density
        energy_exact_squared += weight * energy_exact_density
    end

    energy_absolute = sqrt(max(energy_error_squared, zero(T)))
    energy_exact = sqrt(max(energy_exact_squared, zero(T)))
    energy_relative = energy_absolute / max(energy_exact, eps(T))

    return (;
        component_absolute,
        component_relative,
        combined_absolute,
        combined_relative,
        energy_absolute,
        energy_relative,
    )
end

# -----------------------------------------------------------------------------
# Characteristic speed and time-step estimate
# -----------------------------------------------------------------------------
function _principal_symbol(material::ThermoelasticMaterial, nx::Real, nz::Real)
    T = eltype(material.C)
    An = _normal_matrix(convert(T, nx), convert(T, nz))
    B = zeros(T, NSTATE, NSTATE)

    # Ordering [s(3), v(2), theta, q(2)].
    B[1:3, 4:5] .= material.C * An
    B[4:5, 1:3] .= An' / material.rho
    B[4:5, 6] .= -(An' * material.b) / material.rho
    B[6, 4:5] .= -(material.Tref / material.heat_capacity) .* vec(material.b' * An)
    B[6, 7] = -nx / material.heat_capacity
    B[6, 8] = -nz / material.heat_capacity
    B[7:8, 6] .= -(material.K * T[nx, nz]) / material.tau

    return B
end

"""
    maximum_characteristic_speed(material; nsamples=128)

Compute the largest absolute eigenvalue of the complete coupled principal
symbol, sampled over propagation directions.  For the isotropic Carcione
material this returns approximately 3980 m/s.
"""
function maximum_characteristic_speed(material::ThermoelasticMaterial; nsamples::Int=128)
    nsamples >= 4 || throw(ArgumentError("nsamples must be at least 4"))
    cmax = 0.0
    max_imaginary = 0.0

    for j in 0:(nsamples - 1)
        angle = 2pi * j / nsamples
        B = _principal_symbol(material, cos(angle), sin(angle))
        values = eigvals(B)
        cmax = max(cmax, maximum(abs.(values)))
        max_imaginary = max(max_imaginary, maximum(abs, imag.(values)))
    end

    if max_imaginary > 1e-8 * max(1.0, cmax)
        @warn "The sampled principal symbol has non-negligible complex eigenvalues" max_imaginary
    end

    return cmax
end

"""
    estimate_dt(material, rd, md; CFL=0.25, relaxation_fraction=0.5)

Return a conservative explicit time step for the square affine mesh used in
the Carcione experiment.  The wave restriction uses `(N+1)^2`, while the
local Cattaneo relaxation restriction is a fraction of `tau`.
"""
function estimate_dt(
    material::ThermoelasticMaterial,
    rd,
    md;
    CFL::Real=0.25,
    relaxation_fraction::Real=0.5,
    speed_samples::Int=128,
)
    CFL > 0 || throw(ArgumentError("CFL must be positive"))
    relaxation_fraction > 0 || throw(ArgumentError("relaxation_fraction must be positive"))

    # For the affine square cells in main.jl, J = (h/2)^2.
    h = 2 * minimum(sqrt.(abs.(md.J)))
    cmax = maximum_characteristic_speed(material; nsamples=speed_samples)
    degree_factor = (rd.N + 1)^2

    dt_wave = CFL * h / (degree_factor * cmax)
    dt_relax = relaxation_fraction * material.tau
    dt = min(dt_wave, dt_relax)

    info = (;
        h,
        cmax,
        degree_factor,
        dt_wave,
        dt_relax,
        limiting=(dt_wave <= dt_relax ? :wave : :relaxation),
    )
    return dt, info
end

function compute_cfl(dt::Real, material::ThermoelasticMaterial, rd, md)
    h = 2 * minimum(sqrt.(abs.(md.J)))
    cmax = maximum_characteristic_speed(material)
    return (rd.N + 1)^2 * cmax * dt / h
end

# -----------------------------------------------------------------------------
# Diagnostics and output
# -----------------------------------------------------------------------------
"""
    energy_components(u, rd, md, material)

Evaluate the semi-discrete energy components with StartUpDG volume
quadrature.  The result includes the Cattaneo relaxation dissipation rate.
"""
function energy_components(u, rd, md, material::ThermoelasticMaterial)
    (; Vq) = rd
    (; wJq) = md

    @views begin
        sxxq = Vq * u[:, :, ISXX]
        szzq = Vq * u[:, :, ISZZ]
        sxzq = Vq * u[:, :, ISXZ]
        vxq = Vq * u[:, :, IVX]
        vzq = Vq * u[:, :, IVZ]
        thetaq = Vq * u[:, :, ITHETA]
        qxq = Vq * u[:, :, IQX]
        qzq = Vq * u[:, :, IQZ]
    end

    S = material.S
    epsxxq = @. S[1, 1] * sxxq + S[1, 2] * szzq + S[1, 3] * sxzq
    epszzq = @. S[2, 1] * sxxq + S[2, 2] * szzq + S[2, 3] * sxzq
    gammaxzq = @. S[3, 1] * sxxq + S[3, 2] * szzq + S[3, 3] * sxzq

    Kinv = material.Kinv
    Kiqx = @. Kinv[1, 1] * qxq + Kinv[1, 2] * qzq
    Kiqz = @. Kinv[2, 1] * qxq + Kinv[2, 2] * qzq

    elastic_density = @. 0.5 * (sxxq * epsxxq + szzq * epszzq + sxzq * gammaxzq)
    kinetic_density = @. 0.5 * material.rho * (vxq^2 + vzq^2)
    temperature_density = @. 0.5 * material.heat_capacity / material.Tref * thetaq^2
    heat_flux_density = @. 0.5 * material.tau / material.Tref * (qxq * Kiqx + qzq * Kiqz)
    relaxation_density = @. (qxq * Kiqx + qzq * Kiqz) / material.Tref

    elastic = sum(wJq .* elastic_density)
    kinetic = sum(wJq .* kinetic_density)
    temperature = sum(wJq .* temperature_density)
    heat_flux = sum(wJq .* heat_flux_density)
    relaxation = sum(wJq .* relaxation_density)
    total = elastic + kinetic + temperature + heat_flux

    return (; elastic, kinetic, temperature, heat_flux, total, relaxation)
end


"""
    energy_balance(u, du, rd, md, material;
                   source_supply=0, source_spatial=nothing)

Evaluate the instantaneous semi-discrete energy identity

    dE/dt + D_rel + D_mech + D_BR1 = P_source.

`du` must be the result of `rhs_thermoelastic_br1!(du,u,...)` at the same
state and time.  Interior-face integrals are assembled from both element
sides and divided by two.  The returned `residual` should be near roundoff
for a periodic constant-coefficient SBP/DGSEM mesh, up to time-integration
and floating-point errors in the supplied state.
"""
function energy_balance(
    u,
    du,
    rd,
    md,
    material::ThermoelasticMaterial;
    source_supply::Real=0.0,
    source_spatial=nothing,
)
    (; Vq, Vf, wf) = rd
    (; wJq, Jf, nx, ny, mapP) = md

    @views begin
        sxxq = Vq * u[:, :, ISXX]
        szzq = Vq * u[:, :, ISZZ]
        sxzq = Vq * u[:, :, ISXZ]
        vxq = Vq * u[:, :, IVX]
        vzq = Vq * u[:, :, IVZ]
        thetaq = Vq * u[:, :, ITHETA]
        qxq = Vq * u[:, :, IQX]
        qzq = Vq * u[:, :, IQZ]

        dsxxq = Vq * du[:, :, ISXX]
        dszzq = Vq * du[:, :, ISZZ]
        dsxzq = Vq * du[:, :, ISXZ]
        dvxq = Vq * du[:, :, IVX]
        dvzq = Vq * du[:, :, IVZ]
        dthetaq = Vq * du[:, :, ITHETA]
        dqxq = Vq * du[:, :, IQX]
        dqzq = Vq * du[:, :, IQZ]
    end

    S = material.S
    Kinv = material.Kinv

    epsxxq = @. S[1, 1] * sxxq + S[1, 2] * szzq + S[1, 3] * sxzq
    epszzq = @. S[2, 1] * sxxq + S[2, 2] * szzq + S[2, 3] * sxzq
    gammaxzq = @. S[3, 1] * sxxq + S[3, 2] * szzq + S[3, 3] * sxzq

    Kiqx = @. Kinv[1, 1] * qxq + Kinv[1, 2] * qzq
    Kiqz = @. Kinv[2, 1] * qxq + Kinv[2, 2] * qzq

    energy_rate_density = @. (
        epsxxq * dsxxq +
        epszzq * dszzq +
        gammaxzq * dsxzq +
        material.rho * (vxq * dvxq + vzq * dvzq) +
        (material.heat_capacity / material.Tref) * thetaq * dthetaq +
        (material.tau / material.Tref) * (Kiqx * dqxq + Kiqz * dqzq)
    )

    denergy_dt = sum(wJq .* energy_rate_density)
    relaxation = sum(wJq .* (@. (qxq * Kiqx + qzq * Kiqz) / material.Tref))

    # Reconstruct total stress and face jumps.
    @views begin
        sigxx = @. u[:, :, ISXX] - material.b[1] * u[:, :, ITHETA]
        sigzz = @. u[:, :, ISZZ] - material.b[2] * u[:, :, ITHETA]
        sigxz = @. u[:, :, ISXZ] - material.b[3] * u[:, :, ITHETA]

        sigxxf = Vf * sigxx
        sigzzf = Vf * sigzz
        sigxzf = Vf * sigxz
        vxf = Vf * u[:, :, IVX]
        vzf = Vf * u[:, :, IVZ]
        thetaf = Vf * u[:, :, ITHETA]
    end

    dsigxx = sigxxf[mapP] - sigxxf
    dsigzz = sigzzf[mapP] - sigzzf
    dsigxz = sigxzf[mapP] - sigxzf
    dvx = vxf[mapP] - vxf
    dvz = vzf[mapP] - vzf
    dtheta = thetaf[mapP] - thetaf

    dtx = @. nx * dsigxx + ny * dsigxz
    dtz = @. nx * dsigxz + ny * dsigzz

    # A_n[[v]] in engineering-shear Voigt ordering.
    anv1 = @. nx * dvx
    anv2 = @. ny * dvz
    anv3 = @. ny * dvx + nx * dvz

    mechanical_density = @. (
        0.5 * material.alpha_sigma * (dtx^2 + dtz^2) +
        0.5 * material.alpha_v * (anv1^2 + anv2^2 + anv3^2)
    )

    thermal_density = @. (
        (material.eta_br1 / material.Tref) * dtheta^2
    )

    face_weights = reshape(wf, :, 1) .* Jf
    # Every periodic interior face is represented once from each neighboring
    # element, hence the factor 1/2.
    mechanical_face = 0.5 * sum(face_weights .* mechanical_density)
    thermal_face = 0.5 * sum(face_weights .* thermal_density)

    source_power = zero(denergy_dt)
    if source_spatial !== nothing && !iszero(source_supply)
        sourceq = Vq * source_spatial
        source_power = source_supply / material.Tref * sum(wJq .* thetaq .* sourceq)
    end

    residual = denergy_dt + relaxation + mechanical_face + thermal_face - source_power

    return (;
        denergy_dt,
        relaxation,
        mechanical_face,
        thermal_face,
        source_power,
        residual,
    )
end

function interpolate_state(rd, u)
    nplot = size(rd.Vp, 1)
    K = size(u, 2)
    result = Array{eltype(u)}(undef, nplot, K, NSTATE)
    @views for field in 1:NSTATE
        mul!(result[:, :, field], rd.Vp, u[:, :, field])
    end
    return result
end

"""
    export_quad_subcells_vtu(filename, x, z, uplot, material)

Export the eight primary fields together with the reconstructed total stress.
`x`, `z`, and `uplot` must already be interpolated with `rd.Vp`.
"""
function export_quad_subcells_vtu(
    filename::AbstractString,
    x::AbstractMatrix,
    z::AbstractMatrix,
    uplot::AbstractArray,
    material::ThermoelasticMaterial,
)
    size(x) == size(z) || throw(DimensionMismatch("x and z must have the same size"))
    ndims(uplot) == 3 || throw(DimensionMismatch("uplot must be a three-dimensional array"))

    Np, K = size(x)
    size(uplot, 1) == Np || throw(DimensionMismatch("uplot point count does not match x"))
    size(uplot, 2) == K || throw(DimensionMismatch("uplot element count does not match x"))
    size(uplot, 3) == NSTATE || throw(DimensionMismatch("uplot must contain $NSTATE fields"))

    n1 = round(Int, sqrt(Np))
    n1 * n1 == Np || throw(ArgumentError("plotting points per element must be a perfect square"))

    npoints = Np * K
    points = Matrix{Float64}(undef, 3, npoints)
    points[1, :] .= vec(x)
    points[2, :] .= vec(z)
    points[3, :] .= 0.0

    cells = MeshCell[]
    for k in 1:K
        offset = (k - 1) * Np
        ids = reshape(collect(1:Np), n1, n1)
        for j in 1:(n1 - 1), i in 1:(n1 - 1)
            a = offset + ids[i, j]
            b = offset + ids[i + 1, j]
            c = offset + ids[i + 1, j + 1]
            d = offset + ids[i, j + 1]
            push!(cells, MeshCell(VTKCellTypes.VTK_TRIANGLE, [a, b, c]))
            push!(cells, MeshCell(VTKCellTypes.VTK_TRIANGLE, [a, c, d]))
        end
    end

    @views begin
        sxx = uplot[:, :, ISXX]
        szz = uplot[:, :, ISZZ]
        sxz = uplot[:, :, ISXZ]
        theta = uplot[:, :, ITHETA]

        sigxx = @. sxx - material.b[1] * theta
        sigzz = @. szz - material.b[2] * theta
        sigxz = @. sxz - material.b[3] * theta

        base = endswith(lowercase(filename), ".vtu") ? filename[1:(end - 4)] : filename
        vtk_grid(base, points, cells) do vtk
            vtk["elastic_s_xx", VTKPointData()] = vec(sxx)
            vtk["elastic_s_zz", VTKPointData()] = vec(szz)
            vtk["elastic_s_xz", VTKPointData()] = vec(sxz)
            vtk["total_sigma_xx", VTKPointData()] = vec(sigxx)
            vtk["total_sigma_zz", VTKPointData()] = vec(sigzz)
            vtk["total_sigma_xz", VTKPointData()] = vec(sigxz)
            vtk["v_x", VTKPointData()] = vec(uplot[:, :, IVX])
            vtk["v_z", VTKPointData()] = vec(uplot[:, :, IVZ])
            vtk["theta", VTKPointData()] = vec(theta)
            vtk["q_x", VTKPointData()] = vec(uplot[:, :, IQX])
            vtk["q_z", VTKPointData()] = vec(uplot[:, :, IQZ])
        end

        return base * ".vtu"
    end
end

end # module DG
