using NFFT3
using Test
using LinearAlgebra

N = (16, 8, 4)
M = 10000

### NFFT Test ###

X = rand(3, M) .- 0.5
fhat = rand(prod(N)) + im * rand(prod(N))

#create plan
p = NFFT(N, M)
p.x = X
p.fhat = fhat

NFFT3.nfft_trafo(p)
f2 = p.f

#generate correctly ordered index set
I = [[j; i; k] for k = -N[3]/2:N[3]/2-1, i = -N[2]/2:N[2]/2-1, j = -N[1]/2:N[1]/2-1]
I = vec(I)

#generate Fourier matrix
F = [exp(-2 * pi * im * sum(X[:, j]' * I[l])) for j = 1:M, l = 1:prod(N)]

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F * fhat

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

NFFT3.nfft_adjoint(p)
f2 = p.fhat

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F' * p.f

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

### NFCT Test ###

X = 0.5 .* rand(3, M)
fhat = rand(prod(N))

#create plan
p = NFCT(N, M)
p.x = X
p.fhat = fhat

NFFT3.nfct_trafo(p)
f2 = p.f

#generate correctly ordered index set
I = [[j; i; k] for k = 0:N[3]-1, i = 0:N[2]-1, j = 0:N[1]-1]
I = vec(I)

#generate Fourier matrix
F = [
    cos(2 * pi * X[:, j][1] * I[l][1]) *
    cos(2 * pi * X[:, j][2] * I[l][2]) *
    cos(2 * pi * X[:, j][3] * I[l][3]) for j = 1:M, l = 1:prod(N)
]

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F * fhat

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

NFFT3.nfct_adjoint(p)
f2 = p.fhat

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F' * p.f

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

### NFST Test ###

fhat = rand(prod(collect(N) .- 1))

#create plan
p = NFST(N, M)
p.x = X
p.fhat = fhat

NFFT3.nfst_trafo(p)
f2 = p.f

#generate correctly ordered index set
I = [[j; i; k] for k = 1:N[3]-1, i = 1:N[2]-1, j = 1:N[1]-1]
I = vec(I)

#generate Fourier matrix
F = [
    sin(2 * pi * X[:, j][1] * I[l][1]) *
    sin(2 * pi * X[:, j][2] * I[l][2]) *
    sin(2 * pi * X[:, j][3] * I[l][3]) for j = 1:M, l = 1:prod(collect(N) .- 1)
]

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F * fhat

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

NFFT3.nfst_adjoint(p)
f2 = p.fhat

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F' * p.f

error_vector = f1 - f2
E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(fhat, 1)

@test E_2 < 10^(-10)
@test E_infty < 10^(-10)

### fastsum Test ###

d = 2
N = 20000
M = 20000
kernel = "multiquadric"
c = 1 / sqrt(N)
p = 8
flags = 0
m = p
n = 256
eps_I = p / n
eps_B = max(1 / 16, p / n)
nn = 2 * n

# create a Plan-Object in Julia
plan = FASTSUM(d, N, M, n, p, kernel, c, eps_I, eps_B, nn, m)

# generate source nodes in circle of radius 0.25-eps_B/2
r = sqrt.(rand(N)) .* (0.25 - eps_B / 2)
phi = rand(N) .* (2 * pi)
X = [(r .* cos.(phi)) (r .* sin.(phi))]
plan.x = X

# generate coefficients alpha_k
alpha = rand(N) + im * rand(N)
plan.alpha = alpha

# generate target nodes in circle of radius 0.25-eps_B/2
r = sqrt.(rand(M)) .* (0.25 - eps_B / 2)
phi = rand(M) .* (2 * pi)
Y = [(r .* cos.(phi)) (r .* sin.(phi))]
plan.y = Y

NFFT3.fastsum_trafo(plan)
f1 = copy(plan.f)

NFFT3.fastsum_trafo_exact(plan)
f2 = copy(plan.f)

error_vector = f1 - f2

E_2 = norm(error_vector) / norm(f1)
E_infty = norm(error_vector, Inf) / norm(plan.alpha, 1)

@test E_2 < 10^(-5)
@test E_infty < 10^(-5)
