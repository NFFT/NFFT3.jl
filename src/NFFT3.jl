module NFFT3

export NFFT, NFCT, NFST, FASTSUM

# file ending for OS
ending = ".so"

if Sys.iswindows()
    ending = ".dll"
elseif Sys.isapple()
    ending = ".dylib"
end

const lib_path_nfft = string(@__DIR__, "/libnfftjulia", ending)
const lib_path_nfct = string(@__DIR__, "/libnfctjulia", ending)
const lib_path_nfst = string(@__DIR__, "/libnfstjulia", ending)
const lib_path_fastsum = string(@__DIR__, "/libfastsumjulia", ending)

# flags

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
NFCT_SORT_NODES = UInt32(1) << 11
NFCT_OMP_BLOCKWISE_ADJOINT = UInt32(1) << 12
NFST_SORT_NODES = UInt32(1) << 11
NFST_OMP_BLOCKWISE_ADJOINT = UInt32(1) << 12
PRE_ONE_PSI = (PRE_LIN_PSI | PRE_FG_PSI | PRE_PSI | PRE_FULL_PSI)

# FFTW flags
FFTW_MEASURE = UInt32(0)
FFTW_DESTROY_INPUT = UInt32(1) << 0
FFTW_UNALIGNED = UInt32(1) << 1
FFTW_CONSERVE_MEMORY = UInt32(1) << 2
FFTW_EXHAUSTIVE = UInt32(1) << 3
FFTW_PRESERVE_INPUT = UInt32(1) << 4
FFTW_PATIENT = UInt32(1) << 5
FFTW_ESTIMATE = UInt32(1) << 6
FFTW_WISDOM_ONLY = UInt32(1) << 21

# default flag values
f1_default_1d = UInt32(
    PRE_PHI_HUT |
    PRE_PSI |
    MALLOC_X |
    MALLOC_F_HAT |
    MALLOC_F |
    FFTW_INIT |
    FFT_OUT_OF_PLACE,
)
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
f2_default = UInt32(FFTW_ESTIMATE | FFTW_DESTROY_INPUT)

# default window cut off
default_window_cut_off =
    ccall(("nfft_get_default_window_cut_off", lib_path_nfft), Int64, ())

include("NFFT.jl")
include("NFCT.jl")
include("NFST.jl")
include("fastsum.jl")

end # module
