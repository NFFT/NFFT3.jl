"""
`NFFT3.jl` – Nonequispaced Fast Fourier Transform
"""
module NFFT3

using Aqua
using CpuId
using LinearAlgebra
using LinearMaps

ending = ".so"
path = "/lib/"
glibcversion = ""

if Sys.iswindows()
    ending = ".dll"
    path = "\\lib\\"
elseif Sys.isapple()
    ending = ".dylib"
else
    glibcversion = "glibc2.40/"
    if VersionNumber(unsafe_string(@ccall string(@__DIR__, path, "glibc-version.so").glibc_version()::Cstring)) < v"2.35"
        glibcversion = "glibc2.22/"
    end
end

if cpufeature(:AVX2)
    flag = "AVX2/"
elseif cpufeature(:AVX)
    flag = "AVX/"
else
    flag = "SSE2/"
end

lib_path_nfft = string(@__DIR__, path, flag, glibcversion, "libnfftjulia", ending)
lib_path_nfct = string(@__DIR__, path, flag, glibcversion, "libnfctjulia", ending)
lib_path_nfst = string(@__DIR__, path, flag, glibcversion, "libnfstjulia", ending)
lib_path_fastsum = string(@__DIR__, path, flag, glibcversion, "libfastsumjulia", ending)

include("NFFT.jl")
include("NFCT.jl")
include("NFST.jl")
include("NFMT.jl")
include("fastsum.jl")
include("flags.jl")

kernels = Dict(
    "gaussian" => 1,
    "multiquadric" => 1,
    "inverse_multiquadric" => 1,
    "logarithm" => 0,
    "thinplate_spline" => 0,
    "one_over_square" => 0,
    "one_over_modulus" => 0,
    "one_over_x" => 0,
    "inverse_multiquadric3" => 1,
    "sinc_kernel" => 1,
    "cosc" => 1,
    "kcot" => 1,
    "one_over_cube" => 0,
    "log_sin" => 1,
    "laplacian_rbf" => 1,
    "der_laplacian_rbf" => 1,
    "xx_gaussian" => 1,
    "absx" => 0,
)

# plan structures
export NFFT, NFCT, NFST, FASTSUM, NFMT

# functions
export nfft_finalize_plan,
    nfft_init,
    nfft_trafo,
    nfft_adjoint,
    nfft_trafo_direct,
    nfft_adjoint_direct,
    nfft_get_LinearMap,
    nfft_get_coefficient_vector,
    nfft_get_coefficient_array,
    nfct_finalize_plan,
    nfct_init,
    nfct_trafo,
    nfct_adjoint,
    nfct_transposed,
    nfct_trafo_direct,
    nfct_adjoint_direct,
    nfct_transposed_direct,
    nfct_get_LinearMap,
    nfct_get_coefficient_vector,
    nfct_get_coefficient_array,
    nfst_finalize_plan,
    nfst_init,
    nfst_trafo,
    nfst_adjoint,
    nfst_trafo_direct,
    nfst_adjoint_direct,
    nfst_get_LinearMap,
    nfst_get_coefficient_vector,
    nfst_get_coefficient_array,
    nfmt_finalize_plan,
    nfmt_trafo,
    nfmt_adjoint,
    nfmt_get_LinearMap,
    nfmt_get_coefficient_vector,
    nfmt_get_coefficient_array,
    fastsum_finalize_plan,
    fastsum_init,
    fastsum_trafo,
    fastsum_trafo_exact,
    finalize_plan,
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
