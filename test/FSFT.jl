# NOTE on FSFT and *_direct(): nfsft_trafo_direct()/nfsft_adjoint_direct()
# (shared with NFSFT) unconditionally read the node array for node positions.
# FSFT has no x setter of its own (the whole point of NFSFT_EQUISPACED is
# that the grid is implicit in N) - fsft_init() sets it to the equispaced
# grid itself (mirroring the "Set equispaced nodes" code in
# nfsft_init_guru(), kernel/nfsft/nfsft.c), since the bundled library does
# not reliably do this on its own (verified while porting this from the
# equivalent pyNFFT3 fix: two separately-constructed FSFT plans with the
# same N disagreed on trafo_direct() output in roughly half of repeated
# runs without this fix, 0/30 with it). That makes trafo_direct()/
# adjoint_direct() trustworthy cross-plan references for any N.

# IMPORTANT ordering note: NFSFT/FSFT share one process-global "wisdom"
# precompute (see nfsft_precompute()/wisdom.initialized in
# kernel/nfsft/nfsft.c). It is computed once per process and sized for
# whichever bandwidth N is used *first* - nfsft_precompute() is a no-op on
# every later call, even for a larger N. Any later plan whose N exceeds that
# first N silently gets NaN from trafo()/adjoint() (wisdom.N_MAX check). So
# the larger bandwidth must run first in this file.

function fsft_random_fhat(plan, N)
    fhat = zeros(ComplexF64, (2 * plan.N + 2)^2)
    for k = 0:N
        for n = -k:k
            fhat[fsft_index(plan, k, n)] = (rand() - 0.5) + im * (rand() - 0.5)
        end
    end
    return fhat
end

function fsft_valid_entries(plan, N)
    # Only indices (k, n) with |n| <= k <= N correspond to actual spherical
    # Fourier coefficients; without NFSFT_ZERO_F_HAT (not part of the
    # default flags), adjoint()/adjoint_direct() leave every other entry in
    # the padded N_total-sized fhat array as unspecified scratch content, so
    # comparisons must only look at these indices.
    out = zeros(ComplexF64, length(plan.fhat))
    for k = 0:N
        for n = -k:k
            idx = fsft_index(plan, k, n)
            out[idx] = plan.fhat[idx]
        end
    end
    return out
end

# --- at/above NFSFT_BREAK_EVEN (5): trafo()/adjoint() use a distinct,
# --- FFT-based fast algorithm here, genuinely independent of trafo_direct().
N = 6

plan = FSFT(N)
fhat = fsft_random_fhat(plan, N)
plan.fhat = fhat

NFFT3.trafo_direct(plan)
f1 = copy(plan.f)

# nfsft_default does not set NFSFT_PRESERVE_F_HAT, so fsft_trafo() overwrites
# plan.fhat in place with internal FPT/FFT scratch data (expected behaviour
# of the underlying library) - so trafo_direct() (which only reads fhat)
# must run before trafo() on a shared plan, as done here.
NFFT3.trafo(plan)
f2 = copy(plan.f)
@test norm(f1 - f2) / norm(f1) < 10^(-9)

NFFT3.adjoint_direct(plan)
f3 = fsft_valid_entries(plan, N)

NFFT3.adjoint(plan)
f4 = fsft_valid_entries(plan, N)
@test norm(f3 - f4) / norm(f3) < 10^(-9)

# --- below NFSFT_BREAK_EVEN: trafo()/adjoint() literally call
# --- trafo_direct()/adjoint_direct() internally (see kernel/nfsft/nfsft.c),
# --- so this is an exact-agreement check, not just an accuracy check. Must
# --- run after the N=6 block above (see ordering note).
N_small = 3

plan_small = FSFT(N_small)
fhat_small = fsft_random_fhat(plan_small, N_small)
plan_small.fhat = fhat_small

NFFT3.trafo_direct(plan_small)
f1_small = copy(plan_small.f)

NFFT3.trafo(plan_small)
f2_small = copy(plan_small.f)
@test f1_small == f2_small

NFFT3.adjoint_direct(plan_small)
f3_small = fsft_valid_entries(plan_small, N_small)

NFFT3.adjoint(plan_small)
f4_small = fsft_valid_entries(plan_small, N_small)
@test f3_small == f4_small

# Error tests
@test_throws DomainError FSFT(-1)
@test_throws ErrorException FSFT(N, 0, nfsft_default)
