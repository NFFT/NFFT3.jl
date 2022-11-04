include("/home/arch/git/NFFT3.jl/src/NFFCT.jl")

function getphi(x::Vector{Float64}, k::Vector{Int64}, P::NFFCT{D}) where {D}
    p = 1
    for i = 1:D
        if (P.dcos[i])
            if k[i] ≠ 0
                p *= sqrt(2.0)*cos(pi*k[i]*x[i])
            end
        else
            p *= exp(-2.0*pi*im*k[i]*x[i])
        end
    end
    return p
end

function getMat(P::NFFCT{D}) where {D}
    a = length(findall(P.dcos))
    p = prod(P.NFFT_struct.N)
    n = p÷(2^a)
    X = copy(P.NFFT_struct.x)
    for i in range(1, D)
        if (P.dcos[i])
            X[i,:] .*= 2
        end
    end
    F = [getphi([X[j]], getk(P,l), P) for j = 1:P.NFFT_struct.M, l = 1:n]
    return F
end 

function getk(P::NFFCT{D}, i::Int64) where {D}
    k = zeros(Int64, D)
    pv = zeros(Int64, D)
    pv[D] = 1
    for i = D-1:-1:1
        pv[i] = pv[i+1] * P.NFFT_struct.N[i+1]
        if (P.dcos[i+1])
            pv[i] ÷= 2;
        end
    end
    for j=1:D
        if (P.dcos[j])
            k[j] = ((i - 1) ÷ pv[j]) % (P.NFFT_struct.N[j] ÷ 2)
        else
            k[j] = ((i - 1) ÷ pv[j]) % P.NFFT_struct.N[j] - P.NFFT_struct.N[j] ÷ 2
        end
    end
    return k
end

function test()
    FC = NFFCT((true,), (12,), 2)
    #FC.m = Int32(32)
    #C = NFCT((6,6), 20)
    x = [0.1,0.8]
    #x=rand(5)
    fhat = [1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im]
    #fhat=rand(6)+im*rand(6)
    #s = sqrt(2)
    #xh = 0.5 .* x
    #wu =[1,s,s,s,s,s,s,2,2,2,2,2,s,2,2,2,2,2,s,2,2,2,2,2,s,2,2,2,2,2,s,2,2,2,2,2]
    #fhhat = wu.*fhat
    FC.x = x
    #C.x = xh
    FC.fhat = fhat
    #C.fhat = fhhat

    #F = FC.NFFT_struct
    #println(F)
    #println(F.x)
    #println(unsafe_wrap(Vector{Float64}, Core.getfield(F, :x), F.M))
    #println(F.fhat)
    #println(unsafe_wrap(Vector{ComplexF64}, Core.getfield(F, :fhat), prod(F.N)))
    #println(F.N)
    #println(F.M)
    #println(F.n)
    #println(F.m)
    #println(F.f1)
    #println(F.f2)
    #println(F.init_done)
    #println(F.finalized)
    #F2 = NFFT{1}(F.N, F.M, F.n, F.m, F.f1, F.f2)
    #F2.x = F.x
    #F2.fhat = F.fhat
    #nfft_trafo_direct(F)
    #f1  = copy(F.f)

    #println(f1)

    #println(F.fhat)

    #nfft_trafo(F)
    #f2 = copy(F.f)

    #println(f2)

    #error_vector = f1 - f2

    #println(error_vector)

    #E_2 = norm(error_vector) / norm(f1)
    #E_infty = norm(error_vector, Inf) / norm(fhat, 1)

    #println(E_2)
    #println(E_infty)

    #nfct_trafo(C)
    #println(C.f)
    #println(FC.f)
    #println(sum(abs.(C.f-FC.f)))
    
    nffct_trafo(FC)

    println(sum(abs.(getMat(FC)*fhat-FC.f)))
end

function test2()
    FC = NFFCT((false,true,true), (30,20,12), 200)
    x=rand(3,200)
    f = rand(200)
    FC.x = x
    FC.f = f
    nffct_adjoint(FC)
    println(sum(abs.(Base.adjoint(getMat(FC)) * f - FC.fhat)))
end

test()