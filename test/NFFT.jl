N = (16, 8, 4)
M = 10000

X = rand(3, M) .- 0.5
fhat = rand(prod(N)) + im * rand(prod(N))

p = NFFT(N, M)
p.x = X
p.fhat = fhat


NFFT3.trafo(p)
f2 = p.f

I = [[j; i; k] for k = -N[3]/2:N[3]/2-1, i = -N[2]/2:N[2]/2-1, j = -N[1]/2:N[1]/2-1]
I = vec(I)

F = [exp(-2 * pi * im * sum(X[:, j]' * I[l])) for j = 1:M, l = 1:prod(N)]

f1 = F * fhat

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

NFFT3.adjoint(p)
f2 = p.fhat

f1 = F' * p.f

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

# Error tests
@test_throws DomainError NFFT((-1, 2), M)

@test_throws DomainError NFFT((3, 2), M)

@test_throws DomainError NFFT(N, -1)

@test_throws DomainError NFFT(N, M, (-1, 2, 3), Int32(8), f1_default, f2_default)

# Question: What if just one integer n_i in the tuple is smaler than N_i
@test_throws DomainError NFFT(N, M, (16, 8, 4), Int32(8), f1_default, f2_default)

@test_throws DomainError NFFT(N, M, (18, 9, 6), Int32(8), f1_default, f2_default)

@test_throws DomainError NFFT(N, M, (18, 10, 6), Int32(-8), f1_default, f2_default)

@test_throws DomainError NFFT((-16, 8, 4), M, (18, 10, 6), Int32(8), f1_default, f2_default)

NFFT3.finalize_plan(p)

N = (16,)
M = 100

p2 = NFFT(N, M)

@test_throws "NFFT not initialized." NFFT3.finalize_plan(p2)

X = rand(3, M) .- 0.5

@test_throws "type NFFT has no field X" p2.X = X

N = (16,)
M = 100

p = NFFT(N, M)

X = rand(M) .- 0.5
fhat = rand(prod(N)) + im * rand(prod(N))

p = NFFT(N, M)
NFFT3.init(p)
p.x = X
p.fhat = fhat

NFFT3.nfft_trafo(p)

@test_logs (:warn,"You can't modify the C pointer to the NFFT plan.") p.plan = p2.plan