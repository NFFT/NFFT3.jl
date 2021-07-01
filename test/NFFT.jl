N = (16, 8, 4)
M = 10000

X = rand(3, M) .- 0.5
fhat = rand(prod(N)) + im * rand(prod(N))

p = NFFT(N, M)
p.x = X
p.fhat = fhat

NFFT3.nfft_trafo(p)
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

NFFT3.nfft_adjoint(p)
f2 = p.fhat

f1 = F' * p.f

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)