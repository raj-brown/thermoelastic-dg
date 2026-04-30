module DG

using LinearAlgebra
using WriteVTK
using ProgressMeter, Printf


export RHSCacheTE
export rhs_thermoelastic_ldg!
export rhs_thermoelastic_cldg!
export rick
export export_quad_subcells_vtu
export make_progress_callback
export compute_cfl
export compute_dt
using OrdinaryDiffEq



# State:
# u[:,:,1] = σxx
# u[:,:,2] = σyy
# u[:,:,3] = σxy
# u[:,:,4] = vx
# u[:,:,5] = vy
# u[:,:,6] = T
# u[:,:,7] = ψ

struct RHSCacheTE
    # elastic derivatives
    s11r; s11s; s11x; s11y
    s12r; s12s; s12x; s12y
    s22r; s22s; s22x; s22y
    v1r;  v1s;  v1x;  v1y
    v2r;  v2s;  v2x;  v2y

    # thermal LDG
    Tr; Ts; Tx; Ty
    qxr; qxs; qxx
    qyr; qys; qyy
    lapT

    # Pi derivatives
    Pix; Piy
    Pir; Pis; Pix_x
    Pjr; Pjs; Piy_y
    divPi

    # face elastic
    s11m; s11p; ds11
    s12m; s12p; ds12
    s22m; s22p; ds22
    v1m;  v1p;  dv1
    v2m;  v2p;  dv2

    # face thermal
    Tm; Tp; dT
    qxm; qxp; dqx
    qym; qyp; dqy

    # elastic flux/work
    divsx; divsy; vxy
    nSx; nSy
    nv11x; nv11y; nvxy
    fpenalty_s11; fpenalty_s22; fpenalty_s12
    fpenalty_v1;  fpenalty_v2
    fluxS11; fluxS22; fluxS12
    fluxv1; fluxv2

    # thermal flux/work
    fluxT
    fluxq

    # buffers
    liftbuf
    tmp1
    tmp2
    tmp3
end


function RHSCacheTE(u, rd)
    Np, K = size(u, 1), size(u, 2)
    NfpK = size(rd.Vf, 1)

    vol  = zeros(Np, K)
    face = zeros(NfpK, K)

    return RHSCacheTE(
        # elastic derivatives
        copy(vol), copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol), copy(vol),

        # thermal LDG
        copy(vol), copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol),
        copy(vol),

        # Pi
        copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol),
        copy(vol),

        # face elastic
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),

        # face thermal
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),

        # elastic flux/work
        copy(vol), copy(vol), copy(vol),
        copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face),

        # thermal flux/work
        copy(face), copy(face),

        # buffers
        copy(face), copy(vol), copy(vol), copy(vol)
    )
end


function rhs_thermoelastic_ldg!(du, u, parameters, time)

    # ========================================================
    # O(1) nondimensional verification parameters
    # ========================================================
    ρ = 1.0

    μ = 1.0
    λ = 1.0

    c11 = λ + 2μ
    c22 = c11
    c12 = λ
    c33 = μ

    invrho = 1.0 / ρ

    cT = 1.0
    γ  = 1.0
    T0 = 1.0
    β  = 0.05
    τ  = 1.0

    q_source = 0.0

    (; rd, md, pt_src, cache) = parameters
    (; Vf, Dr, Ds, LIFT) = rd
    (; rxJ, sxJ, ryJ, syJ, nxJ, nyJ, nx, ny, J, Jf, mapP) = md

    c = cache

    hmin = minimum(sqrt.(J))
    αT = 5.0 * (rd.N + 1)^2 / hmin

    s11 = @view u[:, :, 1]
    s22 = @view u[:, :, 2]
    s12 = @view u[:, :, 3]
    v1  = @view u[:, :, 4]
    v2  = @view u[:, :, 5]
    T   = @view u[:, :, 6]
    ψ   = @view u[:, :, 7]

    # ========================================================
    # Elastic derivatives
    # ========================================================
    mul!(c.s11r, Dr, s11)
    mul!(c.s11s, Ds, s11)
    @. c.s11x = rxJ * c.s11r + sxJ * c.s11s
    @. c.s11y = ryJ * c.s11r + syJ * c.s11s

    mul!(c.s12r, Dr, s12)
    mul!(c.s12s, Ds, s12)
    @. c.s12x = rxJ * c.s12r + sxJ * c.s12s
    @. c.s12y = ryJ * c.s12r + syJ * c.s12s

    mul!(c.s22r, Dr, s22)
    mul!(c.s22s, Ds, s22)
    @. c.s22x = rxJ * c.s22r + sxJ * c.s22s
    @. c.s22y = ryJ * c.s22r + syJ * c.s22s

    mul!(c.v1r, Dr, v1)
    mul!(c.v1s, Ds, v1)
    @. c.v1x = rxJ * c.v1r + sxJ * c.v1s
    @. c.v1y = ryJ * c.v1r + syJ * c.v1s

    mul!(c.v2r, Dr, v2)
    mul!(c.v2s, Ds, v2)
    @. c.v2x = rxJ * c.v2r + sxJ * c.v2s
    @. c.v2y = ryJ * c.v2r + syJ * c.v2s

    # ========================================================
    # Elastic face values
    # ========================================================
    mul!(c.s11m, Vf, s11)
    mul!(c.s12m, Vf, s12)
    mul!(c.s22m, Vf, s22)
    mul!(c.v1m,  Vf, v1)
    mul!(c.v2m,  Vf, v2)

    @inbounds for i in eachindex(c.s11p)
        p = mapP[i]
        c.s11p[i] = c.s11m[p]
        c.s12p[i] = c.s12m[p]
        c.s22p[i] = c.s22m[p]
        c.v1p[i]  = c.v1m[p]
        c.v2p[i]  = c.v2m[p]
    end

    @. c.ds11 = c.s11p - c.s11m
    @. c.ds12 = c.s12p - c.s12m
    @. c.ds22 = c.s22p - c.s22m
    @. c.dv1  = c.v1p  - c.v1m
    @. c.dv2  = c.v2p  - c.v2m

    # ========================================================
    # Elastic numerical fluxes
    # ========================================================
    @. c.nSx = nx * c.ds11 + ny * c.ds12
    @. c.nSy = nx * c.ds12 + ny * c.ds22

    @. c.nv11x = c.dv1 * nx
    @. c.nv11y = c.dv2 * ny
    @. c.nvxy  = c.dv2 * nx + c.dv1 * ny

    @. c.fpenalty_s11 = c.nSx * nxJ
    @. c.fpenalty_s22 = c.nSy * nyJ
    @. c.fpenalty_s12 = c.nSy * nxJ + c.nSx * nyJ

    @. c.fpenalty_v1 = nxJ * c.nv11x + nyJ * c.nvxy
    @. c.fpenalty_v2 = nxJ * c.nvxy  + nyJ * c.nv11y

    @. c.fluxS11 = c.nv11x + c.fpenalty_s11
    @. c.fluxS22 = c.nv11y + c.fpenalty_s22
    @. c.fluxS12 = c.nvxy  + c.fpenalty_s12

    @. c.fluxv1 = c.nSx + c.fpenalty_v1
    @. c.fluxv2 = c.nSy + c.fpenalty_v2

    # ========================================================
    # Elastic RHS before constitutive law
    # ========================================================
    @. c.liftbuf = c.fluxS11 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 1] = (c.v1x + 0.5 * c.tmp1) / J

    @. c.liftbuf = c.fluxS22 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 2] = (c.v2y + 0.5 * c.tmp1) / J

    @. c.vxy = c.v1y + c.v2x
    @. c.liftbuf = c.fluxS12 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 3] = (c.vxy + 0.5 * c.tmp1) / J

    @. c.divsx = c.s11x + c.s12y
    @. c.divsy = c.s12x + c.s22y

    @. c.liftbuf = c.fluxv1 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 4] = invrho * (c.divsx + 0.5 * c.tmp1) / J

    @. c.liftbuf = c.fluxv2 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 5] = invrho * (c.divsy + 0.5 * c.tmp1) / J

    # ========================================================
    # Π = acceleration and div(Π)
    # ========================================================
    c.Pix .= @view du[:, :, 4]
    c.Piy .= @view du[:, :, 5]

    mul!(c.Pir, Dr, c.Pix)
    mul!(c.Pis, Ds, c.Pix)
    @. c.Pix_x = rxJ * c.Pir + sxJ * c.Pis

    mul!(c.Pjr, Dr, c.Piy)
    mul!(c.Pjs, Ds, c.Piy)
    @. c.Piy_y = ryJ * c.Pjr + syJ * c.Pjs

    @. c.divPi = c.Pix_x + c.Piy_y

    # ========================================================
    # LDG thermal Laplacian
    #
    # q = ∇T
    # ΔT = ∇·q
    #
    # Flux choice:
    # T_hat = T_plus
    # q_hat·n = q_minus·n - αT [T]
    # ========================================================

    # ---------------------------
    # Stage 1: q = ∇T
    # ---------------------------
    mul!(c.Tr, Dr, T)
    mul!(c.Ts, Ds, T)

    @. c.Tx = rxJ * c.Tr + sxJ * c.Ts
    @. c.Ty = ryJ * c.Tr + syJ * c.Ts

    mul!(c.Tm, Vf, T)

    @inbounds for i in eachindex(c.Tp)
        p = mapP[i]
        c.Tp[i] = c.Tm[p]
    end

    @. c.dT = c.Tp - c.Tm

    # Strong-form correction for T_hat = T_plus.
    @. c.liftbuf = c.dT * nxJ * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. c.Tx -= c.tmp1 / J

    @. c.liftbuf = c.dT * nyJ * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. c.Ty += c.tmp1 / J

    # ---------------------------
    # Stage 2: ΔT = ∇·q
    # ---------------------------
    mul!(c.qxr, Dr, c.Tx)
    mul!(c.qxs, Ds, c.Tx)
    @. c.qxx = rxJ * c.qxr + sxJ * c.qxs

    mul!(c.qyr, Dr, c.Ty)
    mul!(c.qys, Ds, c.Ty)
    @. c.qyy = ryJ * c.qyr + syJ * c.qys

    mul!(c.qxm, Vf, c.Tx)
    mul!(c.qym, Vf, c.Ty)

    @inbounds for i in eachindex(c.qxp)
        p = mapP[i]
        c.qxp[i] = c.qxm[p]
        c.qyp[i] = c.qym[p]
    end

    @. c.dqx = c.qxp - c.qxm
    @. c.dqy = c.qyp - c.qym

    # q_hat = q_minus - αT [T] n
    # Therefore (q_hat - q_minus)·n = -αT [T].
    @. c.fluxq = -αT * c.dT

    @. c.liftbuf = c.fluxq * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)

    # Strong-form correction.
    @. c.lapT = c.qxx + c.qyy - c.tmp1 / J

    # ========================================================
    # Final thermoelastic RHS
    # ========================================================
    c.tmp1 .= @view du[:, :, 1]   # e_xx_dot
    c.tmp2 .= @view du[:, :, 2]   # e_yy_dot
    c.tmp3 .= @view du[:, :, 3]   # 2e_xy_dot

    @. du[:, :, 1] = c11 * c.tmp1 + c12 * c.tmp2 - β * ψ
    @. du[:, :, 2] = c12 * c.tmp1 + c22 * c.tmp2 - β * ψ
    @. du[:, :, 3] = c33 * c.tmp3

    # Velocity equations are already assigned in du[:,:,4:5].

    # Temperature equation
    @. du[:, :, 6] = ψ

    # Relaxation equation.
    # Safer verification version: no τ divΠ term.
    # @. du[:, :, 7] =
    #     (γ * c.lapT -
    #      q_source -
    #      T0 * β * (c.tmp1 + c.tmp2)) / cT -
    #     ψ / τ

    # Full PDE version, enable after basic stability check:
    @. du[:, :, 7] =
        (γ * c.lapT -
         q_source -
         T0 * β * ((c.tmp1 + c.tmp2) + τ * c.divPi)) / cT -
        ψ / τ

    @. du[:, :, 5] += pt_src * rick(time, 10.0)

    return nothing
end



function rhs_thermoelastic_cldg!(du, u, parameters, time)

    # ========================================================
    # O(1) nondimensional verification parameters
    # ========================================================
    ρ = 1.0

    μ = 1.0
    λ = 1.0

    c11 = λ + 2μ
    c22 = c11
    c12 = λ
    c33 = μ

    invrho = 1.0 / ρ

    cT = 1.0
    γ  = 1.0
    T0 = 1.0
    β  = 0.05
    τ  = 1.0

    q_source = 0.0

    (; rd, md, pt_src, cache) = parameters
    (; Vf, Dr, Ds, LIFT) = rd
    (; rxJ, sxJ, ryJ, syJ, nxJ, nyJ, nx, ny, J, Jf, mapP) = md

    c = cache

    hmin = minimum(sqrt.(J))
    η = 10.0 * (rd.N + 1)^2 / hmin

    s11 = @view u[:, :, 1]
    s22 = @view u[:, :, 2]
    s12 = @view u[:, :, 3]
    v1  = @view u[:, :, 4]
    v2  = @view u[:, :, 5]
    T   = @view u[:, :, 6]
    ψ   = @view u[:, :, 7]

    # ========================================================
    # Elastic derivatives
    # ========================================================
    mul!(c.s11r, Dr, s11)
    mul!(c.s11s, Ds, s11)
    @. c.s11x = rxJ * c.s11r + sxJ * c.s11s
    @. c.s11y = ryJ * c.s11r + syJ * c.s11s

    mul!(c.s12r, Dr, s12)
    mul!(c.s12s, Ds, s12)
    @. c.s12x = rxJ * c.s12r + sxJ * c.s12s
    @. c.s12y = ryJ * c.s12r + syJ * c.s12s

    mul!(c.s22r, Dr, s22)
    mul!(c.s22s, Ds, s22)
    @. c.s22x = rxJ * c.s22r + sxJ * c.s22s
    @. c.s22y = ryJ * c.s22r + syJ * c.s22s

    mul!(c.v1r, Dr, v1)
    mul!(c.v1s, Ds, v1)
    @. c.v1x = rxJ * c.v1r + sxJ * c.v1s
    @. c.v1y = ryJ * c.v1r + syJ * c.v1s

    mul!(c.v2r, Dr, v2)
    mul!(c.v2s, Ds, v2)
    @. c.v2x = rxJ * c.v2r + sxJ * c.v2s
    @. c.v2y = ryJ * c.v2r + syJ * c.v2s

    # ========================================================
    # Elastic face values
    # ========================================================
    mul!(c.s11m, Vf, s11)
    mul!(c.s12m, Vf, s12)
    mul!(c.s22m, Vf, s22)
    mul!(c.v1m,  Vf, v1)
    mul!(c.v2m,  Vf, v2)

    @inbounds for i in eachindex(c.s11p)
        p = mapP[i]
        c.s11p[i] = c.s11m[p]
        c.s12p[i] = c.s12m[p]
        c.s22p[i] = c.s22m[p]
        c.v1p[i]  = c.v1m[p]
        c.v2p[i]  = c.v2m[p]
    end

    @. c.ds11 = c.s11p - c.s11m
    @. c.ds12 = c.s12p - c.s12m
    @. c.ds22 = c.s22p - c.s22m
    @. c.dv1  = c.v1p  - c.v1m
    @. c.dv2  = c.v2p  - c.v2m

    # ========================================================
    # Elastic numerical fluxes
    # ========================================================
    @. c.nSx = nx * c.ds11 + ny * c.ds12
    @. c.nSy = nx * c.ds12 + ny * c.ds22

    @. c.nv11x = c.dv1 * nx
    @. c.nv11y = c.dv2 * ny
    @. c.nvxy  = c.dv2 * nx + c.dv1 * ny

    @. c.fpenalty_s11 = c.nSx * nxJ
    @. c.fpenalty_s22 = c.nSy * nyJ
    @. c.fpenalty_s12 = c.nSy * nxJ + c.nSx * nyJ

    @. c.fpenalty_v1 = nxJ * c.nv11x + nyJ * c.nvxy
    @. c.fpenalty_v2 = nxJ * c.nvxy  + nyJ * c.nv11y

    @. c.fluxS11 = c.nv11x + c.fpenalty_s11
    @. c.fluxS22 = c.nv11y + c.fpenalty_s22
    @. c.fluxS12 = c.nvxy  + c.fpenalty_s12

    @. c.fluxv1 = c.nSx + c.fpenalty_v1
    @. c.fluxv2 = c.nSy + c.fpenalty_v2

    # ========================================================
    # Elastic RHS before constitutive law
    # ========================================================
    @. c.liftbuf = c.fluxS11 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 1] = (c.v1x + 0.5 * c.tmp1) / J

    @. c.liftbuf = c.fluxS22 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 2] = (c.v2y + 0.5 * c.tmp1) / J

    @. c.vxy = c.v1y + c.v2x
    @. c.liftbuf = c.fluxS12 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 3] = (c.vxy + 0.5 * c.tmp1) / J

    @. c.divsx = c.s11x + c.s12y
    @. c.divsy = c.s12x + c.s22y

    @. c.liftbuf = c.fluxv1 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 4] = invrho * (c.divsx + 0.5 * c.tmp1) / J

    @. c.liftbuf = c.fluxv2 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 5] = invrho * (c.divsy + 0.5 * c.tmp1) / J

    # ========================================================
    # Π = acceleration and div(Π)
    # ========================================================
    c.Pix .= @view du[:, :, 4]
    c.Piy .= @view du[:, :, 5]

    mul!(c.Pir, Dr, c.Pix)
    mul!(c.Pis, Ds, c.Pix)
    @. c.Pix_x = rxJ * c.Pir + sxJ * c.Pis

    mul!(c.Pjr, Dr, c.Piy)
    mul!(c.Pjs, Ds, c.Piy)
    @. c.Piy_y = ryJ * c.Pjr + syJ * c.Pjs

    @. c.divPi = c.Pix_x + c.Piy_y

    # ========================================================
    # CLDG / SIPG thermal Laplacian
    #
    # Δ_h T = div(grad T) with compact numerical flux:
    #
    # q_hat · n = {grad T} · n - η [T]
    #
    # This is more robust than alternating LDG.
    # ========================================================

    # Volume gradient of T
    mul!(c.Tr, Dr, T)
    mul!(c.Ts, Ds, T)

    @. c.Tx = rxJ * c.Tr + sxJ * c.Ts
    @. c.Ty = ryJ * c.Tr + syJ * c.Ts

    # Volume divergence of grad T
    mul!(c.qxr, Dr, c.Tx)
    mul!(c.qxs, Ds, c.Tx)
    @. c.qxx = rxJ * c.qxr + sxJ * c.qxs

    mul!(c.qyr, Dr, c.Ty)
    mul!(c.qys, Ds, c.Ty)
    @. c.qyy = ryJ * c.qyr + syJ * c.qys

    # Face temperature
    mul!(c.Tm, Vf, T)

    @inbounds for i in eachindex(c.Tp)
        p = mapP[i]
        c.Tp[i] = c.Tm[p]
    end

    @. c.dT = c.Tp - c.Tm

    # Face gradients
    mul!(c.qxm, Vf, c.Tx)
    mul!(c.qym, Vf, c.Ty)

    @inbounds for i in eachindex(c.qxp)
        p = mapP[i]
        c.qxp[i] = c.qxm[p]
        c.qyp[i] = c.qym[p]
    end

    # Average normal derivative
    @. c.fluxT =
        0.5 * (
            nx * (c.qxm + c.qxp) +
            ny * (c.qym + c.qyp)
        ) -
        η * c.dT

    # Strong-form surface correction
    @. c.liftbuf = c.fluxT * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)

    @. c.lapT = c.qxx + c.qyy - c.tmp1 / J

    # ========================================================
    # Final thermoelastic RHS
    # ========================================================
    c.tmp1 .= @view du[:, :, 1]   # e_xx_dot
    c.tmp2 .= @view du[:, :, 2]   # e_yy_dot
    c.tmp3 .= @view du[:, :, 3]   # 2e_xy_dot

    @. du[:, :, 1] = c11 * c.tmp1 + c12 * c.tmp2 - β * ψ
    @. du[:, :, 2] = c12 * c.tmp1 + c22 * c.tmp2 - β * ψ
    @. du[:, :, 3] = c33 * c.tmp3

    # Temperature equation
    @. du[:, :, 6] = ψ

    # Relaxation equation: safer version without τ divΠ
    # @. du[:, :, 7] =
    #     (γ * c.lapT -
    #      q_source -
    #      T0 * β * (c.tmp1 + c.tmp2)) / cT -
    #     ψ / τ

    # Full model version:
    @. du[:, :, 7] =
        (γ * c.lapT -
         q_source -
         T0 * β * ((c.tmp1 + c.tmp2) + τ * c.divPi)) / cT -
        ψ / τ

    @. du[:, :, 5] += pt_src * rick(time, 10.0)

    return nothing
end

function compute_cfl(c, dt, md, rd)
    h = minimum(sqrt.(md.J))
    return (rd.N + 1)^2 * c * dt / h
end


function compute_dt(c11, ρ, τ, md; CFL_target = 0.02)
    # element size
    h = minimum(sqrt.(md.J))

    # wave speed
    c = sqrt(c11 / ρ)

    # time step limits
    dt_elastic = CFL_target * h / c
    dt_relax   = CFL_target * τ

    dt = min(dt_elastic, dt_relax)

    return dt, h, c
end



function rick(t, f0)
    tR = 1 / f0
    return 1e1 * (1 .- 2 * (π * f0 * (t .- tR)) .^ 2) .* exp.(-(π * f0 * (t .- tR)) .^ 2)
end


function export_quad_subcells_vtu(filename, x, y, u)
    @assert size(x) == size(y)
    @assert ndims(u) == 3

    Np, K = size(x)
    @assert size(u, 1) == Np
    @assert size(u, 2) == K

    q = round(Int, sqrt(Np))
    @assert q * q == Np "Np must be a perfect square."
    n1 = q

    npts = Np * K
    points = Matrix{Float64}(undef, 3, npts)

    field_names = ["sigma_xx", "sigma_yy", "sigma_xy", "v1", "v2", "T", "psi"]
    vals = [Vector{Float64}(undef, npts) for _ in 1:7]

    for k in 1:K
        off = (k - 1) * Np

        for i in 1:Np
            gid = off + i

            points[1, gid] = x[i, k]
            points[2, gid] = y[i, k]
            points[3, gid] = 0.0

            for fld in 1:7
                vals[fld][gid] = u[i, k, fld]
            end
        end
    end

    cells = MeshCell[]

    for k in 1:K
        off = (k - 1) * Np
        local_ids = reshape(collect(1:Np), n1, n1)

        for j in 1:n1-1, i in 1:n1-1
            a = off + local_ids[i, j]
            b = off + local_ids[i+1, j]
            c = off + local_ids[i+1, j+1]
            d = off + local_ids[i, j+1]

            # Use the same working indexing convention as your code.
            push!(cells, MeshCell(VTKCellTypes.VTK_TRIANGLE, [a, b, c]))
            push!(cells, MeshCell(VTKCellTypes.VTK_TRIANGLE, [a, c, d]))
        end
    end

    vtk_grid(filename, points, cells) do vtk
        for fld in 1:7
            vtk[field_names[fld], VTKPointData()] = vals[fld]
        end
    end

    return nothing
end


function make_progress_callback(tspan, dt, CFL)

    nsteps = ceil(Int, (tspan[2] - tspan[1]) / dt)
    prog = Progress(nsteps; desc="Solving", showspeed=true)

    step_counter = Ref(0)
    t_start = time()

    function progress_cb(integrator)
        t = integrator.t

        step_counter[] += 1
        elapsed = time() - t_start
        steps_per_sec = step_counter[] / max(elapsed, eps())
        eta = (nsteps - step_counter[]) / max(steps_per_sec, eps())

        ProgressMeter.next!(prog; showvalues=[
            (:t, @sprintf("%.3e", t)),
            (:CFL, @sprintf("%.3e", CFL)),
            (:steps_per_sec, @sprintf("%.1f", steps_per_sec)),
            (:ETA, @sprintf("%.1f", eta))
        ])
        return nothing
    end

    return DiscreteCallback((u,t,integrator)->true, progress_cb;
                            save_positions=(false,false))
end

end
