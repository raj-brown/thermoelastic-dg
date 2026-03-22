using StaticArrays


function rhs!(dU, u, parameters, time)
    Np, K = parameters.rd.Np, parameters.md.K
    N = Np*K
    s11 = reshape(@view(u[1:N]), Np, K)
    s12 = reshape(@view(u[N+1:2N]), Np, K)
    s22 = reshape(@view(u[2N+1:3N]), Np, K)
    v1 = reshape(@view(u[3N+1:4N]), Np, K)
    v2 = reshape(@view(u[4N+1:5N]), Np, K)
   
    ds11 = reshape(@view(dU[1:N]), Np, K)
    ds12 = reshape(@view(dU[N+1:2N]), Np, K)
    ds22 = reshape(@view(dU[2N+1:3N]), Np, K)
    dv1 = reshape(@view(dU[3N+1:4N]), Np, K)
    dv2 = reshape(@view(dU[4N+1:5N]), Np, K)

    ds11 .= 0.0
    ds12 .= 0.0
    ds22 .= 0.0
    dv1 .= 0.0
    dv2 .= 0.0
    nothing 
end


function point_force(t)
    println("Adding point force at time $t")
end

function initial_condition(Np, K)
    s11 = zeros(Np, K)
    s12 = zeros(Np, K)
    s22 = zeros(Np, K)
    v1 = zeros(Np, K)
    v2 = zeros(Np, K)
    vcat(vec(s11), vec(s12), vec(s22), vec(v1), vec(v2))
end

