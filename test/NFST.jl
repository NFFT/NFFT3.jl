N = (16, 8, 4)
M = 10000
X = 0.5 .* rand(3, M)

fhat = rand(prod(collect(N) .- 1))

p = NFST(N, M)
p.x = X
p.fhat = fhat

f3 = p * nfst_get_coefficient_array(fhat, p)

L = nfst_get_LinearMap(collect(N), X)

f4 = L * fhat

NFFT3.nfst_trafo(p)
f2 = p.f

I = [[j; i; k] for k = 1:N[3]-1, i = 1:N[2]-1, j = 1:N[1]-1]
I = vec(I)

F = [
    sin(2 * pi * X[:, j][1] * I[l][1]) *
    sin(2 * pi * X[:, j][2] * I[l][2]) *
    sin(2 * pi * X[:, j][3] * I[l][3]) for j = 1:M, l = 1:prod(collect(N) .- 1)
]

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

@test E_2 < 10^(-8)
@test E_infty < 10^(-8)

f3 = nfst_get_coefficient_vector(p' * p.f)

f4 = L' * p.f

NFFT3.nfst_adjoint(p)
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

@test E_2 < 10^(-8)
@test E_infty < 10^(-8)

# Error tests
@test_throws DomainError NFST((-1, 2), M)

@test_throws DomainError NFST(N, -1)

@test_throws DomainError NFST(N, M, (-1, 2, 3), Int32(8), f1_default, f2_default)

# Question: What if just one integer n_i in the tuple is smaler than N_i
@test_throws DomainError NFST(N, M, (16, 8, 4), Int32(8), f1_default, f2_default)

@test_throws DomainError NFST(N, M, (18, 9, 6), Int32(8), f1_default, f2_default)

@test_throws DomainError NFST(N, M, (18, 10, 6), Int32(-8), f1_default, f2_default)

@test_throws DomainError NFST((-16, 8, 4), M, (18, 10, 6), Int32(8), f1_default, f2_default)