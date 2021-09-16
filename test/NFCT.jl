N = (16, 8, 4)
M = 10000

X = 0.5 .* rand(3, M)
fhat = rand(prod(N))

p = NFCT(N, M)
p.x = X
p.fhat = fhat

NFFT3.nfct_trafo(p)
f2 = p.f

I = [[j; i; k] for k = 0:N[3]-1, i = 0:N[2]-1, j = 0:N[1]-1]
I = vec(I)

F = [
    cos(2 * pi * X[:, j][1] * I[l][1]) *
    cos(2 * pi * X[:, j][2] * I[l][2]) *
    cos(2 * pi * X[:, j][3] * I[l][3]) for j = 1:M, l = 1:prod(N)
]

f1 = F * fhat

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

NFFT3.nfct_adjoint(p)
f2 = p.fhat

f1 = F' * p.f

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

# DomainError tests
@test_throws DomainError NFCT((-1, 2), M)

@test_throws DomainError NFCT(N, -1)

@test_throws DomainError NFCT(N, M, (-1, 2, 3), Int32(8), f1_default, f2_default)

# Question: What if just one integer n_i in the tuple is smaler than N_i
@test_throws DomainError NFCT(N, M, (16, 8, 4), Int32(8), f1_default, f2_default)

@test_throws DomainError NFCT(N, M, (18, 9, 6), Int32(8), f1_default, f2_default)

@test_throws DomainError NFCT(N, M, (18, 10, 6), Int32(-8), f1_default, f2_default)

@test_throws DomainError NFCT((-16, 8, 4), M, (18, 10, 6), Int32(8), f1_default, f2_default)

# error tests

N = (16, 8, 4)
M = 10000

p = NFCT(N, M)

p.init_done = false

# @test_throws error nfct_finalize_plan(p)