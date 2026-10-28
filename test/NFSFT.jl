N = 8
M = 100

# pseudo-random nodes
X = zeros(2, M)
X[1, :] = rand(M) .- 0.5
X[2, :] = rand(M) .* 0.5

# test init and setting x
plan = NFSFT(N, M)
plan.x = X

# Regression test: jnfsft_set_x (julia/nfsft/libnfsftjulia.c) expects raw
# angles and divides by 2*pi itself to arrive at the plan's native
# coordinates - x must come back unchanged (in the documented normalized
# [-0.5,0.5)x[0,1) convention), not scaled down by another 2*pi.
@test plan.x ≈ X

# generate pseudo-random spherical Fourier coefficients
fhat = zeros(ComplexF64, (2 * plan.N + 2)^2)
for k = 0:N
    for n = (-k):k
        fhat[nfsft_index(plan, k, n)] = (rand() - 0.5) + im * (rand() - 0.5)
    end
end
plan.fhat = fhat

# test trafo direct
NFFT3.trafo_direct(plan)
f1 = copy(plan.f)

# test trafo
# NOTE: nfsft_default does not include NFSFT_PRESERVE_F_HAT, so trafo()
# overwrites plan.fhat in place with internal FPT/FFT scratch data as a side
# effect (see kernel/nfsft/nfsft.c: nfsft_trafo() aliases f_hat_intern =
# f_hat instead of copying unless that flag is set). This is expected
# behaviour of the underlying library, not a bug - it also means
# trafo_direct() must run BEFORE trafo() on a shared plan, as done here,
# since trafo_direct() only ever reads fhat.
NFFT3.trafo(plan)
f2 = copy(plan.f)

# test adjoint direct
NFFT3.adjoint_direct(plan)
f3 = zeros(ComplexF64, length(plan.fhat))
for k = 0:N
    for n = (-k):k
        idx = nfsft_index(plan, k, n)
        f3[idx] = plan.fhat[idx]
    end
end

# test fast approximate adjoint
NFFT3.adjoint(plan)
f4 = zeros(ComplexF64, length(plan.fhat))
for k = 0:N
    for n = (-k):k
        idx = nfsft_index(plan, k, n)
        f4[idx] = plan.fhat[idx]
    end
end

# calculate the error vectors
error_vector_traf = f1 - f2
error_vector_adj = f3 - f4

E_2_traf = norm(error_vector_traf) / norm(f1)
E_infty_traf = norm(error_vector_traf, Inf) / norm(fhat, 1)
@test E_2_traf < 10^(-9)
@test E_infty_traf < 10^(-9)

E_2_adj = norm(error_vector_adj) / norm(f3)
E_infty_adj = norm(error_vector_adj, Inf) / norm(f2, 1)
@test E_2_adj < 10^(-9)
@test E_infty_adj < 10^(-9)

# Regression test: for N < NFSFT_BREAK_EVEN (5) with default flags (no
# NFSFT_PRESERVE_F_HAT), nfsft_trafo()'s internal direct-algorithm fallback
# only *reads* f_hat_intern, so fhat is never actually touched - a subsequent
# trafo_direct() call must not be refused as "corrupted".
N_small = 3
plan_small = NFSFT(N_small, M)
plan_small.x = X
fhat_small = zeros(ComplexF64, (2 * plan_small.N + 2)^2)
for k = 0:N_small
    for n = (-k):k
        fhat_small[nfsft_index(plan_small, k, n)] = (rand() - 0.5) + im * (rand() - 0.5)
    end
end
plan_small.fhat = fhat_small
NFFT3.trafo(plan_small)
@test !plan_small.fhat_corrupted_by_fast_trafo
NFFT3.trafo_direct(plan_small)  # must not throw

# Regression test: a custom f1 that omits the NFSFT_MALLOC_* bits must not
# crash nfsft_init()/the setters below - they are always force-added (mirrors
# NFFT{D}'s own constructor forcing MALLOC_X/MALLOC_F_HAT/MALLOC_F).
plan_custom_f1 = NFSFT(N_small, M, NFSFT_PRESERVE_F_HAT)
plan_custom_f1.x = X
plan_custom_f1.fhat = fhat_small
NFFT3.trafo(plan_custom_f1)
NFFT3.trafo_direct(plan_custom_f1)

# Error tests
@test_throws DomainError NFSFT(-1, M)
@test_throws DomainError NFSFT(N, -1)
@test_throws ErrorException NFSFT(N, M, NFSFT_EQUISPACED)
