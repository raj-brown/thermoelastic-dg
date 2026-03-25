struct RHSCache
    # volume arrays
    s11r; s11s; s11x; s11y
    s12r; s12s; s12x; s12y
    s22r; s22s; s22x; s22y
    v1r;  v1s;  v1x;  v1y
    v2r;  v2s;  v2x;  v2y

    # face arrays
    s11m; s11p; ds11
    s12m; s12p; ds12
    s22m; s22p; ds22
    v1m;  v1p;  dv1
    v2m;  v2p;  dv2

    # flux + work
    divsx; divsy; vxy
    nSx; nSy
    nv11x; nv11y; nvxy

    fpenalty_s11; fpenalty_s22; fpenalty_s12
    fpenalty_v1;  fpenalty_v2

    fluxS11; fluxS22; fluxS12
    fluxv1;  fluxv2

    liftbuf
    tmp1; tmp2
end


function RHSCache(u, rd)
    Np, K = size(u,1), size(u,2)
    NfpK = size(rd.Vf,1)   # face nodes × elements

    vol  = zeros(Np, K)
    face = zeros(NfpK, K)

    RHSCache(
        copy(vol), copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol), copy(vol),
        copy(vol), copy(vol), copy(vol), copy(vol),

        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),
        copy(face), copy(face), copy(face),

        copy(vol), copy(vol), copy(vol),
        copy(face), copy(face),
        copy(face), copy(face), copy(face),

        copy(face), copy(face), copy(face),
        copy(face), copy(face),

        copy(face), copy(face), copy(face),
        copy(face), copy(face),

        copy(face),
        copy(vol), copy(vol)
    )
end


function rhs!(du, u, parameters, time)

    μ = 1.0
    λ = 2.0
    c11 = λ + 2μ
    c12 = λ
    c33 = μ
    c22 = c11
    invrho = 1.0

    (; rd, md, cache) = parameters
    (; Vf, Dr, Ds, LIFT) = rd
    (; rxJ, sxJ, ryJ, syJ, nxJ, nyJ, nx, ny, J, Jf, mapP) = md

    c = cache

    s11 = @view u[:, :, 1]
    s22 = @view u[:, :, 2]
    s12 = @view u[:, :, 3]
    v1  = @view u[:, :, 4]
    v2  = @view u[:, :, 5]

    # ===== Derivatives (no allocs)
    mul!(c.s11r, Dr, s11); mul!(c.s11s, Ds, s11)
    @. c.s11x = rxJ * c.s11r + sxJ * c.s11s
    @. c.s11y = ryJ * c.s11r + syJ * c.s11s

    mul!(c.s12r, Dr, s12); mul!(c.s12s, Ds, s12)
    @. c.s12x = rxJ * c.s12r + sxJ * c.s12s
    @. c.s12y = ryJ * c.s12r + syJ * c.s12s

    mul!(c.s22r, Dr, s22); mul!(c.s22s, Ds, s22)
    @. c.s22x = rxJ * c.s22r + sxJ * c.s22s
    @. c.s22y = ryJ * c.s22r + syJ * c.s22s

    mul!(c.v1r, Dr, v1); mul!(c.v1s, Ds, v1)
    @. c.v1x = rxJ * c.v1r + sxJ * c.v1s
    @. c.v1y = ryJ * c.v1r + syJ * c.v1s

    mul!(c.v2r, Dr, v2); mul!(c.v2s, Ds, v2)
    @. c.v2x = rxJ * c.v2r + sxJ * c.v2s
    @. c.v2y = ryJ * c.v2r + syJ * c.v2s

    # ===== Face values
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

    # ===== Physics
    @. c.divsx = c.s11x + c.s12y
    @. c.divsy = c.s12x + c.s22y
    @. c.vxy   = c.v1y + c.v2x

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
    @. c.fluxv1  = c.nSx   + c.fpenalty_v1
    @. c.fluxv2  = c.nSy   + c.fpenalty_v2

    # ===== LIFT
    @. c.liftbuf = c.fluxS11 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 1] = (c.v1x + 0.5 * c.tmp1) / J

    @. c.liftbuf = c.fluxS22 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 2] = (c.v2y + 0.5 * c.tmp1) / J

    @. c.liftbuf = c.fluxS12 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 3] = (c.vxy + 0.5 * c.tmp1) / J

    @. c.liftbuf = c.fluxv1 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 4] = (c.divsx + 0.5 * c.tmp1) / J

    @. c.liftbuf = c.fluxv2 * Jf
    mul!(c.tmp1, LIFT, c.liftbuf)
    @. du[:, :, 5] = (c.divsy + 0.5 * c.tmp1) / J

    # ===== Constitutive (no alloc)
    c.tmp1 .= @view du[:, :, 1]
    c.tmp2 .= @view du[:, :, 2]

    @. du[:, :, 1] = c11 * c.tmp1 + c12 * c.tmp2
    @. du[:, :, 2] = c12 * c.tmp1 + c22 * c.tmp2
    @. du[:, :, 3] = c33 * du[:, :, 3]
    @. du[:, :, 4] = du[:, :, 4] * invrho
    @. du[:, :, 5] = du[:, :, 5] * invrho

    return nothing
end