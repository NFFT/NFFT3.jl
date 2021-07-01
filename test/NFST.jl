fhat = rand(prod(collect(N) .- 1))

p = NFST(N, M)
p.x = X
p.fhat = fhat

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

NFFT3.nfst_adjoint(p)
f2 = p.fhat

f1 = F' * p.f

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)