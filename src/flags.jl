@doc raw"""
    PRE_PHI_HUT

precompute and store values ``\hat{\phi}(k)`` of the Fourier transform of the window function ``\hat{\phi}``.
"""
PRE_PHI_HUT = UInt32(1) << 0
@doc raw"""
    FG_PSI

fast Gaussian gridding.
"""
FG_PSI = UInt32(1) << 1
@doc raw"""
    PRE_LIN_PSI

linear interpolation of the window function from a lookup table.
"""
PRE_LIN_PSI = UInt32(1) << 2
@doc raw"""
    PRE_FG_PSI

fast Gaussian gridding.
"""
PRE_FG_PSI = UInt32(1) << 3
@doc raw"""
    PRE_PSI

precomputation based on tensor product structure of the window function.
"""
PRE_PSI = UInt32(1) << 4
@doc raw"""
    PRE_FULL_PSI

calculate and store all values ``\tilde{\psi}(x_j - \frac{1}{n} \odot l)``.
"""
PRE_FULL_PSI = UInt32(1) << 5
@doc raw"""
    MALLOC_X

allocate memory for node ``x_j``.
"""
MALLOC_X = UInt32(1) << 6
@doc raw"""
    MALLOC_F_HAT

allocate memory for coefficient``\hat{f}_k``.
"""
MALLOC_F_HAT = UInt32(1) << 7
@doc raw"""
    MALLOC_F

allocate memory for approximate function value ``f_j``.
"""
MALLOC_F = UInt32(1) << 8
@doc raw"""
    FFT_OUT_OF_PLACE

FFTW uses disjoint input/output vector.
"""
FFT_OUT_OF_PLACE = UInt32(1) << 9
@doc raw"""
    FFTW_INIT

initialize FFTW plan.
"""
FFTW_INIT = UInt32(1) << 10
@doc raw"""
    NFFT_SORT_NODES

internal sorting of the nodes ``x_j`` that may increase performance.
"""
NFFT_SORT_NODES = UInt32(1) << 11
@doc raw"""
    NFFT_OMP_BLOCKWISE_ADJOINT

blockwise calculation for adjoint NFFT in the case of OpenMP.
"""
NFFT_OMP_BLOCKWISE_ADJOINT = UInt32(1) << 12
@doc raw"""
    NFCT_SORT_NODES

internal sorting of the nodes ``x_j`` that may increase performance.
"""
NFCT_SORT_NODES = UInt32(1) << 11
@doc raw"""
    NFCT_OMP_BLOCKWISE_ADJOINT

blockwise calculation for adjoint NFFT in the case of OpenMP.
"""
NFCT_OMP_BLOCKWISE_ADJOINT = UInt32(1) << 12
@doc raw"""
    NFST_SORT_NODES

internal sorting of the nodes ``x_j`` that may increase performance.
"""
NFST_SORT_NODES = UInt32(1) << 11
@doc raw"""
    NFST_OMP_BLOCKWISE_ADJOINT

blockwise calculation for adjoint NFFT in the case of OpenMP.
"""
NFST_OMP_BLOCKWISE_ADJOINT = UInt32(1) << 12
@doc raw"""
    PRE_ONE_PSI
"""
PRE_ONE_PSI = (PRE_LIN_PSI | PRE_FG_PSI | PRE_PSI | PRE_FULL_PSI)

# FFTW flags
@doc raw"""
    FFTW_MEASURE

find optimal plan by executing several FFTs and compare times.
"""
FFTW_MEASURE = UInt32(0)
@doc raw"""
    FFTW_DESTROY_INPUT

an out-of-place transform is allowed to overwrite the input array with arbitrary data.
"""
FFTW_DESTROY_INPUT = UInt32(1) << 0
@doc raw"""
    FFTW_UNALIGNED

the algorithm may not impose any unusual alignment requirements on the input/output arrays (not necessary in most context).
"""
FFTW_UNALIGNED = UInt32(1) << 1
@doc raw"""
    FFTW_CONSERVE_MEMORY

conserving memory.
"""
FFTW_CONSERVE_MEMORY = UInt32(1) << 2
@doc raw"""
    FFTW_EXHAUSTIVE

behaves like FFTW_PATIENT with an even wider range of tests.
"""
FFTW_EXHAUSTIVE = UInt32(1) << 3
@doc raw"""
    FFTW_PRESERVE_INPUT

input vector is preserved and unchanged.
"""
FFTW_PRESERVE_INPUT = UInt32(1) << 4
@doc raw"""
    FFTW_PATIENT

behaves like FFTW_MEASURE with a wider range of tests.
"""
FFTW_PATIENT = UInt32(1) << 5
@doc raw"""
    FFTW_ESTIMATE

use simple heuristic instead of measurements to pick a plan.
"""
FFTW_ESTIMATE = UInt32(1) << 6
@doc raw"""
    FFTW_WISDOM_ONLY

a plan is only created if wisdom from tests is available.
"""
FFTW_WISDOM_ONLY = UInt32(1) << 21

# default flag values
@doc raw"""
    f1_default_1d
"""
f1_default_1d = UInt32(
    PRE_PHI_HUT |
    PRE_PSI |
    MALLOC_X |
    MALLOC_F_HAT |
    MALLOC_F |
    FFTW_INIT |
    FFT_OUT_OF_PLACE,
)
@doc raw"""
    f1_default
"""
f1_default = UInt32(
    PRE_PHI_HUT |
    PRE_PSI |
    MALLOC_X |
    MALLOC_F_HAT |
    MALLOC_F |
    FFTW_INIT |
    FFT_OUT_OF_PLACE |
    NFCT_SORT_NODES |
    NFCT_OMP_BLOCKWISE_ADJOINT,
)
@doc raw"""
    f1_default
"""
f2_default = UInt32(FFTW_ESTIMATE | FFTW_DESTROY_INPUT)

# default window cut off
default_window_cut_off =
    ccall(("nfft_get_default_window_cut_off", lib_path_nfft), Int64, ())

BASES = Dict("exp"=>0,"cos"=>1,"alg"=>2)

# NFSFT flags
@doc raw"""
    NFSFT_NORMALIZED

multiply Fourier coefficients with normalization weight.
"""
NFSFT_NORMALIZED = UInt32(1) << 0
@doc raw"""
    NFSFT_USE_NDFT

the fast NFSFT algorithms will use internally the exact but usually slower direct NDFT algorithm in favor of the fast but approximative NFFT algorithm.
"""
NFSFT_USE_NDFT = UInt32(1) << 1
@doc raw"""
    NFSFT_USE_DPT

the fast NFSFT algorithms will use internally the usually slower direct DPT algorithm in favor of the fast FPT algorithm.
"""
NFSFT_USE_DPT = UInt32(1) << 2
@doc raw"""
    NFSFT_MALLOC_X

the init methods will allocate memory and nfsft_finalize will free the array x.
"""
NFSFT_MALLOC_X = UInt32(1) << 3
@doc raw"""
    NFSFT_MALLOC_F_HAT

the init methods will allocate memory and nfsft_finalize will free the array fhat.
"""
NFSFT_MALLOC_F_HAT = UInt32(1) << 5
@doc raw"""
    NFSFT_MALLOC_F

the init methods will allocate memory and nfsft_finalize will free the array f.
"""
NFSFT_MALLOC_F = UInt32(1) << 6
@doc raw"""
    NFSFT_PRESERVE_F_HAT

guarantee that during an execution of nfsft_trafo_direct or nfsft_trafo, the content of fhat remains unchanged.
"""
NFSFT_PRESERVE_F_HAT = UInt32(1) << 7
@doc raw"""
    NFSFT_PRESERVE_X

guarantee that during an execution of nfsft_trafo_direct, nfsft_trafo, nfsft_adjoint_direct or nfsft_adjoint, the content of x remains unchanged.
"""
NFSFT_PRESERVE_X = UInt32(1) << 8
@doc raw"""
    NFSFT_PRESERVE_F

guarantee that during an execution of nfsft_adjoint_direct or nfsft_adjoint, the content of f remains unchanged.
"""
NFSFT_PRESERVE_F = UInt32(1) << 9
@doc raw"""
    NFSFT_DESTROY_F_HAT

explicitly allow that during an execution of nfsft_trafo_direct or nfsft_trafo, the content of fhat may be changed.
"""
NFSFT_DESTROY_F_HAT = UInt32(1) << 10
@doc raw"""
    NFSFT_DESTROY_X

explicitly allow that during an execution of nfsft_trafo_direct, nfsft_trafo, nfsft_adjoint_direct or nfsft_adjoint, the content of x may be changed.
"""
NFSFT_DESTROY_X = UInt32(1) << 11
@doc raw"""
    NFSFT_DESTROY_F

explicitly allow that during an execution of nfsft_adjoint_direct or nfsft_adjoint, the content of f may be changed.
"""
NFSFT_DESTROY_F = UInt32(1) << 12
@doc raw"""
    NFSFT_NO_DIRECT_ALGORITHM

the transforms nfsft_trafo_direct and nfsft_adjoint_direct do not work, saving some memory for precomputed data.
"""
NFSFT_NO_DIRECT_ALGORITHM = UInt32(1) << 13
@doc raw"""
    NFSFT_NO_FAST_ALGORITHM

the transforms nfsft_trafo and nfsft_adjoint do not work, saving memory for precomputed data.
"""
NFSFT_NO_FAST_ALGORITHM = UInt32(1) << 14
@doc raw"""
    NFSFT_ZERO_F_HAT

the transforms nfsft_adjoint and nfsft_adjoint_direct set all unused entries in fhat not corresponding to spherical Fourier coefficients to zero.
"""
NFSFT_ZERO_F_HAT = UInt32(1) << 16
@doc raw"""
    NFSFT_EQUISPACED

use the equispaced FFT (FSFT) instead of the NFFT.
"""
NFSFT_EQUISPACED = UInt32(1) << 17

@doc raw"""
    nfsft_default

Default NFSFT flags.

Note: without NFSFT_PRESERVE_F_HAT (not part of this default), the
underlying C library's fast nfsft_trafo() overwrites the plan's fhat buffer
in place with internal FPT/FFT scratch data as a side effect (it treats
fhat as reusable workspace unless told to copy it first) - so fhat no
longer holds the original coefficients after a trafo() call. This is
expected behaviour of the underlying library, not a bug; pass
NFSFT_PRESERVE_F_HAT explicitly if you need fhat to survive a trafo() call.
See kernel/nfsft/nfsft.c in the NFFT/nfft repository (jnfsft_trafo ->
nfsft_trafo) for the underlying code.
"""
nfsft_default = UInt32(NFSFT_MALLOC_X | NFSFT_MALLOC_F | NFSFT_MALLOC_F_HAT)
@doc raw"""
    nfsft_nfft_default

Default NFFT flags for NFSFT calculations.
"""
nfsft_nfft_default = UInt32(PRE_PHI_HUT | PRE_PSI | FFTW_INIT | NFFT_OMP_BLOCKWISE_ADJOINT)
@doc raw"""
    nfsft_default_nfft_cut_off

The default NFFT cutoff parameter for NFSFT calculations.
"""
nfsft_default_nfft_cut_off = 6
@doc raw"""
    nfsft_default_threshold

The default threshold for the FPT.
"""
nfsft_default_threshold = 1000