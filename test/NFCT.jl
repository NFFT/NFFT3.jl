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