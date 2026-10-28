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

# Regression test for d=1: the x/y setters used to validate the distance-to-0
# bound via maximum(val) instead of maximum(abs.(val)), so out-of-bounds
# *negative* nodes were never rejected (the d>1 branch already used norm.(),
# which is non-negative by construction, so this only affected d=1).
d1 = 1
N1 = 200
M1 = 200
kernel1 = "gaussian"
c1 = 1 / sqrt(N1)

plan1 = FASTSUM(d1, N1, M1, kernel1, c1)

x1 = (rand(N1) .- 0.5) .* (0.25 - eps_B / 2)
plan1.x = x1

alpha1 = rand(N1) + im * rand(N1)
plan1.alpha = alpha1

y1 = (rand(M1) .- 0.5) .* (0.25 - eps_B / 2)
plan1.y = y1

NFFT3.fastsum_trafo(plan1)
f1_1d = copy(plan1.f)

NFFT3.fastsum_trafo_exact(plan1)
f2_1d = copy(plan1.f)

error_vector_1d = f1_1d - f2_1d
E_2_1d = norm(error_vector_1d) / norm(f1_1d)
@test E_2_1d < 10^(-5)

# A sufficiently negative node must still be rejected, same as a positive one.
plan2 = FASTSUM(d1, N1, M1, kernel1, c1)
x_bad = fill(-(0.25 - eps_B / 2) - 0.01, N1)
@test_throws ErrorException (plan2.x = x_bad)