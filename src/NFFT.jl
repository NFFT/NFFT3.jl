# dummy struct for C
mutable struct nfft_plan end

# NFFT plan struct
@doc raw"""
    NFFT{D}

A NFFT (nonequispaced fast Fourier transform) plan, where D is the dimension. 

The classic FFT (Fast Fourier transform) algorithm computes the discrete Fourier transform

```math
f_j \colon = \sum^{\frac{N}{2} - 1}_{k = - \frac{N}{2}} \hat{f}_k e^{-2 \pi \mathrm{i} \frac{k_j}{N}}
```

for ``j = - \frac{N}{2}, \dots, \frac{N}{2} - 1`` and given complex coefficients ``\hat{f}_{k} \in \mathbb{C}``. Using a divide and conquer approach, the number of floating point operations is reduced from ``\mathcal {O}(N^2)`` for a straightforward computation to only ``\mathcal {O}(N \log N)``. 
However, two shortcomings of traditional schemes are the need for equispaced sampling and the restriction to the system of complex exponential functions. The NFFT overcomes the need for an equispaced sampling grid. Considering a D-dimensional trigonometric polynomial

```math
f \colon \mathbb{T}^D \to \mathbb{C}, \; f(x) \colon = \sum_{k \in I_{N}} \hat{f}_k e^{-2 \pi \mathrm{i} \mathbf{k} \cdot \mathbf{x}}
```

with an index set ``I_N \colon = \{ k \in \mathbb{Z}^{D} \colon - \frac{N_i}{2} \leq k_i \leq \frac{N_i}{2} - 1, \, i = 0, \cdots, D-1 \}`` where ``N \in 2 \mathbb{N}^{D}`` is the multibandlimit. 
The NDFT (non uniform discrete fourier transform) is its evaluation at ``M \in 2 \mathbb{N}`` nonequispaced points ``x_j \in \mathbb{T}^D`` for ``j = 0, 1, \cdots, M``,

```math
f(x_j) \colon = \sum_{k \in I_{N}} \hat{f}_k e^{-2 \pi \mathrm{i} \mathbf{k} \cdot \mathbf{x_j}}
```

with given coefficients ``\hat{f}_k \in \mathbb{C}`` where we identify the smooth manifold of the torus ``\mathbb{T}`` with ``[−1/2, 1/2)``. The NFFT is an algorithm for the fast evaluation of the NDFT and the adjoint problem, the fast evaluation of the adjoint NDFT

```math
\hat{h}_k \colon = \sum^{M-1}_{j = 0} f_j e^{-2 \pi \mathrm{i} \mathbf{k} \cdot \mathbf{x_j}}, \; k \in I_N
```

for given coefficients ``f_j \in \mathbb{C}``. In general, the adjoint NDFT is not the inverse transform of the NDFT.
The available NFFT3 library [^KeinerKunisPotts] provides C routines for the NFFT (for computing the NDFT in one or more dimensions, of arbitrary input size, and of complex data), applications such as the fast evaluation of sums

```math
g(y_j) \colon = \sum^{N}_{k = 1} \alpha_k K(\lVert y_j - x_k \rVert_2), \; j = 1, \dots, M
```

for given coefficients ``\alpha_k \in \mathbb{C}``, nodes ``x_k, y_j \in \mathbb{R}^D`` and a radial kernel function ``K \colon [0, \infty) \to [0, \infty)``, and generalizations such as the NNFFT for nonequispaced nodes in time and frequency domain.

# Fields
* `N` - the multibandlimit.
* `M` - the number of nodes.
* `n` - the oversampling per dimension.
* `m` - the window size. Larger m means more accuracy but also more computational costs. 
* `f1` - the NFFT flags.
* `f2` - the FFTW flags.
* `init_done` - indicates if the plan is initialized.
* `finalized` - indicates if the plan is finalized.
* `x` - the nodes.
* `f` - the function values.
* `fhat` - the Fourier coefficients.
* `plan` - plan (C pointer).

# Constructor
    NFFT{D}(N::NTuple{D,Int32},M::Int32,n::NTuple{D,Int32},m::Int32,f1::UInt32,f2::UInt32) where D

# Additional Constructor
    NFFT(N::NTuple{D,Int32},M::Int32,n::NTuple{D,Int32},m::Int32,f1::UInt32,f2::UInt32) where {D}
    NFFT(N::NTuple{D,Int32},M::Int32) where {D}

[^Schmischke2018]:
    > Schmischke, Michael 
    > Nonequispaced Fast Fourier Transform (NFFT) Interface for Julia.
    > 2018
    > url: https://arxiv.org/abs/1810.09891

[^PlonkaPottsSteidelTasche2018]:
    > Plonka, Gerlind and Potts, Daniel and Steidl, Gabriele and Tasche, Manfred 
    > Numerical Fourier Analysis
    > 2018

[^KeinerKunisPotts]:
    > J. Keiner, S. Kunis, and D. Potts
    > NFFT 3.0, C subroutine library
    > url: http://www.tu-chemnitz.de/~potts/nfft.
"""
mutable struct NFFT{D}
    N::NTuple{D,Int32}      # bandwidth tuple
    M::Int32                # number of nodes
    n::NTuple{D,Int32}      # oversampling per dimension
    m::Int32                # windows size
    f1::UInt32              # NFFT flags
    f2::UInt32              # FFTW flags
    init_done::Bool         # bool for plan init
    finalized::Bool         # bool for finalizer
    x::Ref{Float64}         # nodes
    f::Ref{ComplexF64}      # function values
    fhat::Ref{ComplexF64}   # Fourier coefficients
    plan::Ref{nfft_plan}    # plan (C pointer)
    function NFFT{D}(
        N::NTuple{D,Int32},
        M::Int32,
        n::NTuple{D,Int32},
        m::Int32,
        f1::UInt32,
        f2::UInt32,
    ) where {D}
        # create plan object
        new(N, M, n, m, f1, f2, false, false)
    end
end

# additional constructor for easy use [NFFT((N,N),M) instead of NFFT{2}((N,N),M)]
function NFFT(N::NTuple{D,Integer}, M::Integer) where {D}
    if any(x -> x <= 0, N)
        error("Every entry of N has to be an even, positive integer.")
    end

    if sum(N .% 2) != 0
        error("Every entry of N has to be an even, positive integer.")
    end

    if M <= 0
        error("M has to be a positive integer.")
    end

    # convert N to vector for passing it over to C
    Nv = collect(N)

    # default oversampling
    n = Array{Int32}(2 .^ (ceil.(log.(Nv) / log(2)) .+ 1))
    n = NTuple{D,Int32}(n)

    # default NFFT flags
    f1 = UInt32(0)

    if D > 1
        f1 = f1_default
    else
        f1 = f1_default_1d
    end

    NFFT{D}(NTuple{D,Int32}(N), Int32(M), n, Int32(8), f1, f2_default)
end

function NFFT(
    N::NTuple{D,Integer},
    M::Integer,
    n::NTuple{D,Integer},
    m::Integer = Int32(default_window_cut_off),
    f1::UInt32 = (D > 1 ? f1_default : f1_default_1d),
    f2::UInt32 = f2_default,
) where {D}
    # safety checks
    if any(x -> x <= 0, N)
        error("Every entry of N has to be an even, positive integer.")
    end

    if sum(N .% 2) != 0
        error("Every entry of N has to be an even, positive integer.")
    end

    if M <= 0
        error("M has to be a positive integer.")
    end

    if any(x -> x <= 0, n)
        error("Every entry of n has to be an even integer.")
    end

    if n <= N
        error("Every entry of n has to be larger than the corresponding entry in N.")
    end

    if sum(n .% 2) != 0
        error("Every entry of n has to be an even integer.")
    end

    if m <= 0
        error("m has to be a positive integer.")
    end

    NFFT{D}(
        NTuple{D,Int32}(N),
        Int32(M),
        NTuple{D,Int32}(n),
        Int32(m),
        (f1 | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT),
        f2,
    )
end

# finalizer
@doc raw"""
    nfft_finalize_plan(P)

destroys a NFFT plan structure.

# Input
* `P` - a NFFT plan structure.

# See also
[`NFFT{D}`](@ref), [`nfft_init`](@ref)
"""
function nfft_finalize_plan(P::NFFT{D}) where {D}
    if !P.init_done
        error("NFFT not initialized.")
    end

    if !P.finalized
        Core.setfield!(P, :finalized, true)
        ccall(("jnfft_finalize", lib_path_nfft), Nothing, (Ref{nfft_plan},), P.plan)
    end
end

function finalize_plan(P::NFFT{D}) where {D}
    return nfft_finalize_plan(P)
end

# allocate plan memory and init with D,N,M,n,m,f1,f2
@doc raw"""
    nfft_init(P)

intialises a transform plan.

# Input
* `P` - a NFFT plan structure.

# See also
[`NFFT{D}`](@ref), [`nfft_finalize_plan`](@ref)
"""
function nfft_init(P::NFFT{D}) where {D}
    # convert N and n to vectors for passing them over to C
    Nv = collect(P.N)
    n = collect(P.n)

    # call init for memory allocation
    ptr = ccall(("jnfft_alloc", lib_path_nfft), Ptr{nfft_plan}, ())

    # set pointer
    Core.setfield!(P, :plan, ptr)

    # initialize values
    ccall(
        ("jnfft_init", lib_path_nfft),
        Nothing,
        (Ref{nfft_plan}, Int32, Ref{Int32}, Int32, Ref{Int32}, Int32, UInt32, UInt32),
        ptr,
        D,
        Nv,
        P.M,
        n,
        P.m,
        P.f1,
        P.f2,
    )
    Core.setfield!(P, :init_done, true)
    finalizer(nfft_finalize_plan, P)
end

function init(P::NFFT{D}) where {D}
    return nfft_init(P)
end

# overwrite dot notation for plan struct in order to use C memory
function Base.setproperty!(P::NFFT{D}, v::Symbol, val) where {D}
    # init plan if not done [usually with setting nodes]
    if !P.init_done
        nfft_init(P)
    end

    # prevent bad stuff from happening
    if P.finalized
        error("NFFT already finalized")
    end

    # setting nodes, verification of correct size dxM
    if v == :x
        if D == 1
            if typeof(val) != Vector{Float64}
                error("x has to be a Float64 vector.")
            end
            if size(val)[1] != P.M
                error("x has to be a Float64 vector of length M.")
            end
        else
            if typeof(val) != Array{Float64,2}
                error("x has to be a Float64 matrix.")
            end
            if size(val)[1] != D || size(val)[2] != P.M
                error("x has to be a Float64 matrix of size dxM.")
            end
        end
        ptr = ccall(
            ("jnfft_set_x", lib_path_nfft),
            Ptr{Float64},
            (Ref{nfft_plan}, Ref{Cdouble}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)
        # setting values
    elseif v == :f
        if typeof(val) != Array{ComplexF64,1}
            error("f has to be a ComplexFloat64 vector.")
        end
        if size(val)[1] != P.M
            error("f has to be a ComplexFloat64 vector of size M.")
        end
        ptr = ccall(
            ("jnfft_set_f", lib_path_nfft),
            Ptr{ComplexF64},
            (Ref{nfft_plan}, Ref{ComplexF64}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)
        # setting Fourier coefficients
    elseif v == :fhat
        if typeof(val) != Array{ComplexF64,1}
            error("fhat has to be a ComplexFloat64 vector.")
        end
        l = prod(P.N)
        if size(val)[1] != l
            error("fhat has to be a ComplexFloat64 vector of size prod(N).")
        end
        ptr = ccall(
            ("jnfft_set_fhat", lib_path_nfft),
            Ptr{ComplexF64},
            (Ref{nfft_plan}, Ref{ComplexF64}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)
        # prevent modification of NFFT plan pointer
    elseif v == :plan
        @warn "You can't modify the C pointer to the NFFT plan."
    elseif v == :num_threads
        @warn "You can't currently modify the number of threads."
    elseif v == :init_done
        @warn "You can't modify this flag."
    elseif v == :N
        @warn "You can't modify the bandwidth, please create an additional plan."
    elseif v == :M
        @warn "You can't modify the number of nodes, please create an additional plan."
    elseif v == :n
        @warn "You can't modify the oversampling parameter, please create an additional plan."
    elseif v == :m
        @warn "You can't modify the window size, please create an additional plan."
    elseif v == :f1
        @warn "You can't modify the NFFT flags, please create an additional plan."
    elseif v == :f2
        @warn "You can't modify the FFTW flags, please create an additional plan."
        # handle other set operations the default way
    else
        Core.setfield!(P, v, val)
    end
end

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(P::NFFT{D}, v::Symbol) where {D}
    if v == :x
        if !isdefined(P, :x)
            error("x is not set.")
        end
        ptr = Core.getfield(P, :x)
        if D == 1
            return unsafe_wrap(Vector{Float64}, ptr, P.M)             # get nodes from C memory and convert to Julia type
        else
            return unsafe_wrap(Matrix{Float64}, ptr, (D, Int64(P.M)))  # get nodes from C memory and convert to Julia type
        end
    elseif v == :num_threads
        return ccall(("nfft_get_num_threads", lib_path_nfft), Int64, ())
    elseif v == :f
        if !isdefined(P, :f)
            error("f is not set.")
        end
        ptr = Core.getfield(P, :f)
        return unsafe_wrap(Vector{ComplexF64}, ptr, P.M)  # get function values from C memory and convert to Julia type
    elseif v == :fhat
        if !isdefined(P, :fhat)
            error("fhat is not set.")
        end
        ptr = Core.getfield(P, :fhat)
        return unsafe_wrap(Vector{ComplexF64}, ptr, prod(P.N)) # get Fourier coefficients from C memory and convert to Julia type
    else
        return Core.getfield(P, v)
    end
end

# nfft trafo direct [call with NFFT.trafo_direct outside module]
@doc raw"""
    nfft_trafo_direct(P)

computes a NFFT.

# Input
* `P` - a NFFT plan structure.

# See also
[`NFFT{D}`](@ref), [`trafo`](@ref)
"""
function nfft_trafo_direct(P::NFFT{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFFT already finalized")
    end

    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end

    if !isdefined(P, :x)
        error("x has not been set.")
    end

    ptr = ccall(
        ("jnfft_trafo_direct", lib_path_nfft),
        Ptr{ComplexF64},
        (Ref{nfft_plan},),
        P.plan,
    )
    Core.setfield!(P, :f, ptr)
end

function trafo_direct(P::NFFT{D}) where {D}
    return nfft_trafo_direct(P)
end

# adjoint trafo direct [call with NFFT.adjoint_direct outside module]
@doc raw"""
    nfft_adjoint_direct(P)

computes an adjoint NFFT.

# Input
* `P` - a NFFT plan structure.

# See also
[`NFFT{D}`](@ref), [`adjoint`](@ref)
"""
function nfft_adjoint_direct(P::NFFT{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFFT already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(
        ("jnfft_adjoint_direct", lib_path_nfft),
        Ptr{ComplexF64},
        (Ref{nfft_plan},),
        P.plan,
    )
    Core.setfield!(P, :fhat, ptr)
end

function adjoint_direct(P::NFFT{D}) where {D}
    return nfft_adjoint_direct(P)
end

# nfft trafo [call with NFFT.trafo outside module]
@doc raw"""
    nfft_trafo(P)

computes a NFFT.

# Input
* `P` - a NFFT plan structure.

# See also
[`NFFT{D}`](@ref), [`nfft_trafo_direct`](@ref)
"""
function nfft_trafo(P::NFFT{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFFT already finalized")
    end
    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(("jnfft_trafo", lib_path_nfft), Ptr{ComplexF64}, (Ref{nfft_plan},), P.plan)
    Core.setfield!(P, :f, ptr)
end

function trafo(P::NFFT{D}) where {D}
    return nfft_trafo(P)
end

# adjoint trafo [call with NFFT.adjoint outside module]
@doc raw"""
    nfft_adjoint(P)

computes an adjoint NFFT.

# Input
* `P` - a NFFT plan structure.

# See also
[`NFFT{D}`](@ref), [`nfft_adjoint_direct`](@ref)
"""
function nfft_adjoint(P::NFFT{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFFT already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr =
        ccall(("jnfft_adjoint", lib_path_nfft), Ptr{ComplexF64}, (Ref{nfft_plan},), P.plan)
    Core.setfield!(P, :fhat, ptr)
end

function adjoint(P::NFFT{D}) where {D}
    return nfft_adjoint(P)
end