using StaticArrays


function rhs!(du_cache, u, parameters, time)
    Np, K = parameters.rd.Np, parameters.md.K
    N = Np*K

    (; rd, md) = parameters
    (; Vf, Dr, Ds, LIFT) = rd
    (; rxJ, sxJ, ryJ, syJ, nxJ, nyJ, nx, ny, J, Jf, mapP) = md
       


    s11 = reshape(@view(u[1:N]), Np, K)
    s12 = reshape(@view(u[N+1:2N]), Np, K)
    s22 = reshape(@view(u[2N+1:3N]), Np, K)
    v1 = reshape(@view(u[3N+1:4N]), Np, K)
    v2 = reshape(@view(u[4N+1:5N]), Np, K)

    uM = Vf * s11
    uP = uM[mapP]


   
   
   
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

