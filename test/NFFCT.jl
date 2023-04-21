using LinearAlgebra
using Test
using NFFT3

function getphi(dcos::NTuple{D,String}, x::Vector{Float64}, k::Vector{Int64})::ComplexF64 where {D}
    p = 1
    for (idx, s) in enumerate(dcos)
        if (NFFT3.BASES[s]==1)
            if k[idx] ≠ 0
                p *= sqrt(2.0)*cos(pi*k[idx]*x[idx])
            end
        elseif (NFFT3.BASES[s]==2)
            if k[idx] ≠ 0
                p *= sqrt(2.0)*cos(k[idx]*acos(2*x[idx]-1))
            end
        else
            p *= exp(-2.0*pi*im*k[idx]*x[idx])
        end
    end
    return p
end

function getMat(N::NTuple{D,Integer},M::Integer, x
    ::Array{Float64},dcos::NTuple{D,String}) where {D}
    a = sum(getindex.([NFFT3.BASES],dcos).>0)
    p = prod(N)
    n = p÷(2^a)
    X = copy(x)

    F = [getphi(dcos,X[:, j], getk(dcos,N,l)) for j = 1:M, l = 1:n]
    return F
end 

function getk(dcos::NTuple{D,String},N::NTuple{D,Integer}, i::Int64) where {D}
    k = zeros(Int64, D)
    pv = zeros(Int64, D)
    pv[D] = 1
    for i = D-1:-1:1
        pv[i] = pv[i+1] * N[i+1]
        if (NFFT3.BASES[dcos[i+1]]>0)
            pv[i] ÷= 2;
        end
    end
    for j=1:D
        if (NFFT3.BASES[dcos[j]]>0)
            k[j] = ((i - 1) ÷ pv[j]) % (N[j] ÷ 2)
        else
            k[j] = ((i - 1) ÷ pv[j]) % N[j] - N[j] ÷ 2
        end
    end
    return k
end

N    = (20,    8,    12,   4)
dcos = ("exp", "alg", "cos", "alg")
M = 100

X = rand(4, M)
a = sum(getindex.([NFFT3.BASES],dcos).>0)
p = prod(N)
l = p÷(2^a)
fhat = rand(prod(l)) + im * rand(prod(l))

p = NFFCT(dcos, N, M)
p.x = X
p.fhat = fhat

#NFFT3.nffct_trafo(p)
nffct_trafo(p)
f2 = p.f

F = getMat(N, M, X,dcos)

f1 = F * fhat

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

#NFFT3.nffct_adjoint(p)
nffct_adjoint(p)
f2 = p.fhat

f1 = F' * p.f

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)


