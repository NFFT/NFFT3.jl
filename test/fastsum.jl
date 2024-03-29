d = 2
N = 20000
M = 20000
kernel = "multiquadric"
c = 1 / sqrt(N)
eps_B = 1/16

# create a Plan-Object in Julia
plan = FASTSUM(d, N, M, kernel, c)

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

# Error tests
@test_throws DomainError FASTSUM(d, -1, M, kernel, c)

@test_throws DomainError FASTSUM(d, N, -1, kernel, c)

@test_throws DomainError FASTSUM(d, N, M, kernel, c, -1)

@test_throws DomainError FASTSUM(d, N, M, kernel, c, 256, 8, 256/8, 1/16, 512, -1)

cv = Vector{Float64}(undef, 1)
cv[1] = Float64(c)

@test_throws DomainError FASTSUM(d, N, M, 1, 1, kernel, cv, 256/8, 1/16, 512, 512, -1, 1, UInt32(0))

@test_throws DomainError FASTSUM(d, N, M, 1, 1, kernel, cv, 256/8, 1/16, 512, 512, 1, -1, UInt32(0))