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
PRE_PHI_HUT = UInt32(1) << 0
FG_PSI = UInt32(1) << 1
PRE_LIN_PSI = UInt32(1) << 2
PRE_FG_PSI = UInt32(1) << 3
PRE_PSI = UInt32(1) << 4
PRE_FULL_PSI = UInt32(1) << 5
MALLOC_X = UInt32(1) << 6
MALLOC_F_HAT = UInt32(1) << 7
MALLOC_F = UInt32(1) << 8
FFT_OUT_OF_PLACE = UInt32(1) << 9
FFTW_INIT = UInt32(1) << 10
NFFT_SORT_NODES = UInt32(1) << 11
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
