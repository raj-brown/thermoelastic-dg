using StaticArrays


function rhs!(du, u, parameters, time)


    μ = 1.0
    λ = 2.0
    c11 = λ + 2μ
    c12 = λ
    c33 = μ
    c22 = c11
    rho = 1.0

    Np, K = parameters.rd.Np, parameters.md.K

    (; rd, md) = parameters
    (; Vf, Dr, Ds, LIFT) = rd
    (; rxJ, sxJ, ryJ, syJ, nxJ, nyJ, nx, ny, J, Jf, mapP) = md
       


    s11 = view(u, :, :, 1)
    s22 = view(u, :, :, 2)
    s12 = view(u, :, :, 3)
    v1 = view(u, :, :, 4)
    v2 = view(u, :, :, 5)
    
    s11r = Dr * s11 
    s11s = Ds * s11
    
    s11x = rxJ .* s11r + sxJ .* s11s
    s11y = ryJ .* s11r + syJ .* s11s
    
    
    s12r = Dr * s12 
    s12s = Ds * s12
    
    s12x = rxJ .* s12r + sxJ .* s12s
    s12y = ryJ .* s12r + syJ .* s12s


    s22r = Dr * s22 
    s22s = Ds * s22
    s22x = rxJ .* s22r + sxJ .* s22s
    s22y = ryJ .* s22r + syJ .* s22s

    v1r = Dr * v1
    v1s = Ds * v1
    v1x = rxJ .* v1r + sxJ .* v1s
    v1y = ryJ .* v1r + syJ .* v1s

    v2r = Dr * v2
    v2s = Ds * v2
    v2x = rxJ .* v2r + sxJ .* v2s
    v2y = ryJ .* v2r + syJ .* v2s


    s11m = Vf * s11
    s11p = s11m[mapP]
    ds11 = s11p - s11m
    
    s12m = Vf * s12
    s12p = s12m[mapP]
    ds12 = s12p - s12m

    s22m = Vf * s22
    s22p = s22m[mapP]
    ds22 = s22p - s22m

    v1m = Vf * v1
    v1p = v1m[mapP]
    dv1 = v1p - v1m

    v2m = Vf * v2
    v2p = v2m[mapP]
    dv2 = v2p - v2m

    divsx = s11x + s12y
    divsy = s12x + s22y
    vxy   = v1y + v2x


    # Compute the velocity fluxes
    nSx = nx .* ds11 + ny.* ds12
    nSy = nx .* ds12 + ny.* ds22


    # Impose BC's


    # Stress fluxes
    nv11x = dv1 .* nx
    nv11y = dv2 .* ny
    nvxy = dv2.* nx + dv1 .* ny

    fcs11 = nv11x
    fcs22 = nv11y
    fcs12 = nvxy

    fcv1 = nSx
    fcv2 = nSy

    # Penalization term

    fpenalty_s11 = fcv1.*nxJ
    fpenalty_s22 = fcv2.*nyJ
    fpenalty_s12 = fcv2.*nxJ + fcv1.*nyJ

    fpenalty_v1 = nxJ.*fcs11 + nyJ.*fcs12
    fpenalty_v2 = nxJ.*fcs12 + nyJ.*fcs22
    
    tau = 1.0
    fluxS11 = @. fcs11 + tau * fpenalty_s11 
    fluxS22 = @. fcs22 + tau * fpenalty_s22
    fluxS12 = @. fcs12 + tau * fpenalty_s12 
    fluxv1 = @. fcv1 + tau * fpenalty_v1
    fluxv2 = @. fcv2 + tau * fpenalty_v2


    du[:, :, 1] .= (v1x + 0.5*LIFT *(fluxS11 .* Jf))./J
    du[:, :, 2] .= (v1y + 0.5*LIFT *(fluxS22 .* Jf))./J
    du[:, :, 3] .= (vxy + 0.5*LIFT *(fluxS12 .* Jf))./J
    du[:, :, 4] .= (divsx + 0.5*LIFT *(fluxv1 .* Jf))./J
    du[:, :, 5] .= (divsy + 0.5* LIFT *(fluxv2 .* Jf))./J

    du[:, :, 1] = c11.* du[:, :, 1] + c12.* du[:, :, 2]
    du[:, :, 2] = c12.* du[:, :, 1] + c22.* du[:, :, 2] 
    du[:, :, 3] = c33.* du[:, :, 3]

    du[:, :, 4] .= du[:, :, 4]./rho
    du[:, :, 5] .= du[:, :, 5]./rho
   
    nothing 
end


function point_force(t)
    println("Adding point force at time $t")
end



