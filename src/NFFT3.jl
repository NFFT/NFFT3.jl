"""
`NFFT3.jl` – 
"""
module NFFT3


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


include("NFFT.jl")
include("NFCT.jl")
include("NFST.jl")
include("fastsum.jl")
include("flags.jl")

# plan structures
export NFFT, NFCT, NFST, FASTSUM

# flags 
export PRE_PHI_HUT,
    FG_PSI,
    PRE_LIN_PSI,
    PRE_FG_PSI,
    PRE_PSI,
    PRE_FULL_PSI,
    MALLOC_X,
    MALLOC_F_HAT,
    MALLOC_F,
    FFT_OUT_OF_PLACE,
    FFTW_INIT,
    NFFT_SORT_NODES,
    NFFT_OMP_BLOCKWISE_ADJOINT,
    NFCT_SORT_NODES,
    NFCT_OMP_BLOCKWISE_ADJOINT,
    NFST_SORT_NODES,
    NFST_OMP_BLOCKWISE_ADJOINT,
    PRE_ONE_PSI,
    FFTW_MEASURE,
    FFTW_DESTROY_INPUT,
    FFTW_UNALIGNED,
    FFTW_CONSERVE_MEMORY,
    FFTW_EXHAUSTIVE,
    FFTW_PRESERVE_INPUT,
    FFTW_PATIENT,
    FFTW_ESTIMATE,
    FFTW_WISDOM_ONLY,
    f1_default_1d,
    f1_default,
    f2_default,
    default_window_cut_off

end # module
