"""
`NFFT3.jl` – Nonequispaced Fast Fourier Transform
"""
module NFFT3

using Aqua
using CpuId

if cpufeature(:AVX2)
    ending = ".so"

    if Sys.iswindows()
        ending = ".dll"
    elseif Sys.isapple()
        ending = ".dylib"
    end

    lib_path_nfft_str = string(@__DIR__, "/libnfftjuliaavx2", ending)
    lib_path_nfct_str = string(@__DIR__, "/libnfctjuliaavx2", ending)
    lib_path_nfst_str = string(@__DIR__, "/libnfstjuliaavx2", ending)
    lib_path_fastsum_str = string(@__DIR__, "/libfastsumjuliaavx2", ending)
else
    ending = ".so"

    if Sys.iswindows()
        error("NFFT needs a CPU that supports the AVX2 instruction set expansion.")
    elseif Sys.isapple()
        error("NFFT needs a CPU that supports the AVX2 instruction set expansion.")
    end
    
    lib_path_nfft_str = string(@__DIR__, "/libnfftjulia", ending)
    lib_path_nfct_str = string(@__DIR__, "/libnfctjulia", ending)
    lib_path_nfst_str = string(@__DIR__, "/libnfstjulia", ending)
    lib_path_fastsum_str = string(@__DIR__, "/libfastsumjulia", ending)
end


const lib_path_nfft = lib_path_nfft_str
const lib_path_nfct = lib_path_nfct_str
const lib_path_nfst = lib_path_nfst_str
const lib_path_fastsum = lib_path_fastsum_str

avx2 = cpufeature(:AVX2)


include("NFFT.jl")
include("NFCT.jl")
include("NFST.jl")
include("NFMT.jl")
include("fastsum.jl")
include("flags.jl")



# plan structures
export NFFT, NFCT, NFST, FASTSUM, NFMT

# functions
export nfft_finalize_plan,
    nfft_init,
    nfft_trafo,
    nfft_adjoint,
    nfft_trafo_direct,
    nfft_adjoint_direct,
    nfct_finalize_plan,
    nfct_init,
    nfct_trafo,
    nfct_adjoint,
    nfct_transposed,
    nfct_trafo_direct,
    nfct_adjoint_direct,
    nfct_transposed_direct,
    nfst_finalize_plan,
    nfst_init,
    nfst_trafo,
    nfst_adjoint,
    nfst_trafo_direct,
    nfst_adjoint_direct,
    fastsum_finalize_plan,
    fastsum_init,
    fastsum_trafo,
    fastsum_trafo_exact,
    finalize_plan,
    nfmt_finalize_plan,
    nfmt_init,
    nfmt_trafo,
    nfmt_adjoint,
    init,
    trafo,
    adjoint,
    trafo_direct,
    adjoint_direct,
    trafo_exact

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
