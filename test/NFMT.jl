using LinearAlgebra
using Test
using NFFT3

function getphi(basis_vect::NTuple{D,String}, x::Vector{Float64}, k::Vector{Int64})::ComplexF64 where {D}
    p = 1
    for (idx, s) in enumerate(basis_vect)
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
    ::Array{Float64},basis_vect::NTuple{D,String}) where {D}
    a = sum(getindex.([NFFT3.BASES],basis_vect).>0)
    p = prod(N)
    n = p÷(2^a)
    X = copy(x)

    F = [getphi(basis_vect,X[:, j], getk(basis_vect,N,l)) for j = 1:M, l = 1:n]
    return F
end 

function getk(basis_vect::NTuple{D,String},N::NTuple{D,Integer}, i::Int64) where {D}
    k = zeros(Int64, D)
    pv = zeros(Int64, D)
    pv[D] = 1
    for i = D-1:-1:1
        pv[i] = pv[i+1] * N[i+1]
        if (NFFT3.BASES[basis_vect[i+1]]>0)
            pv[i] ÷= 2;
        end
    end
    for j=1:D
        if (NFFT3.BASES[basis_vect[j]]>0)
            k[j] = ((i - 1) ÷ pv[j]) % (N[j] ÷ 2)
        else
            k[j] = ((i - 1) ÷ pv[j]) % N[j] - N[j] ÷ 2
        end
    end
    return k
end

N    = (20,    8,    12,   4)
basis_vect = ("exp", "alg", "cos", "alg")
M = 100

X = rand(4, M)
a = sum(getindex.([NFFT3.BASES],basis_vect).>0)
p = prod(N)
l = p÷(2^a)
fhat = rand(prod(l)) + im * rand(prod(l))

p = NFMT(basis_vect, N, M)
p.x = X
p.fhat = fhat

f3 = p * nfmt_get_coefficient_array(fhat, p)

L = nfmt_get_LinearMap(basis_vect, collect(N), X)

f4 = L * fhat

#NFFT3.nfmt_trafo(p)
nfmt_trafo(p)
f2 = p.f

F = getMat(N, M, X,basis_vect)

f1 = F * fhat

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

error_vector = f1 - f3
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

error_vector = f1 - f4
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

# --- real-valued backend: when every dimension is cos/alg (no exp
# --- dimension) and fhat/f is a real Float64 array, NFMT realises the mixed
# --- transform with a plain NFCT instead of an NFFT (half the data, real
# --- arithmetic throughout) - see _nfmt_all_cos_alg()/_nfmt_ensure_backend!
# --- in src/NFMT.jl. A complex fhat on the same basis_vect still uses the
# --- general NFFT path, and a plan that already picked the NFCT backend
# --- must reject complex data afterwards.
N_real = (8, 4)
basis_real = ("cos", "alg")
M_real = 60

X_real = rand(2, M_real)
a_real = sum(getindex.([NFFT3.BASES], basis_real) .> 0)
nf_real = prod(N_real) ÷ (2^a_real)
fhat_real = rand(nf_real)

plan_real = NFMT(basis_real, N_real, M_real)
plan_real.x = X_real
plan_real.fhat = fhat_real
@test plan_real.backend isa NFCT{2}

F_real = getMat(N_real, M_real, X_real, basis_real)
f1_real = real.(F_real * complex.(fhat_real))

nfmt_trafo(plan_real)
@test eltype(plan_real.f) == Float64
@test norm(plan_real.f - f1_real) / norm(f1_real) < 10^(-6)

NFFT3.nfmt_trafo_direct(plan_real)
@test norm(plan_real.f - f1_real) / norm(f1_real) < 10^(-9)

f1_real_adj = real.(F_real' * complex.(plan_real.f))

nfmt_adjoint(plan_real)
@test eltype(plan_real.fhat) == Float64
@test norm(plan_real.fhat - f1_real_adj) / norm(f1_real_adj) < 10^(-6)

NFFT3.nfmt_adjoint_direct(plan_real)
@test norm(plan_real.fhat - f1_real_adj) / norm(f1_real_adj) < 10^(-9)

# complex fhat on an all-cos/alg basis_vect still uses the NFFT backend
plan_mixed_dtype = NFMT(basis_real, N_real, M_real)
plan_mixed_dtype.x = X_real
fhat_complex = rand(nf_real) + im * rand(nf_real)
plan_mixed_dtype.fhat = fhat_complex
@test plan_mixed_dtype.backend isa NFFT{2}
nfmt_trafo(plan_mixed_dtype)
f1_complex = F_real * fhat_complex
@test norm(plan_mixed_dtype.f - f1_complex) / norm(f1_complex) < 10^(-6)

# once a plan has picked the real NFCT backend, complex data must be rejected
@test_throws ErrorException (plan_real.fhat = fhat_complex)

f3 = nfmt_get_coefficient_vector(p' * p.f)

f4 = L' * p.f

#NFFT3.nfmt_adjoint(p)
nfmt_adjoint(p)
f2 = p.fhat

f1 = F' * p.f

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

error_vector = f1 - f3
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

error_vector = f1 - f4
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)
