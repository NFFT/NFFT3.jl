# dummy struct for C
mutable struct nfst_plan end

# NFST plan struct
@doc raw"""
    NFST{D}

A NFST (Nonequispaced fast sine transform) plan, where D is the dimension. 

The NFST realizes a direct and fast computation of the discrete nonequispaced sine transform. The aim is to compute

```math
f^c (x) \colon = \sum_{k \in I_N} \hat{f}^{c}_{k} \sin (k \dot x), \quad x \in \mathbb{R}^D
```

at given (nonequidistant) knots ``x_k \in [0, \pi ]^D, \; k = 0, \cdots, M-1``, coefficients ``\hat{f}^{c}_{k} \in \mathbb{R}``, ``k \in I_N \colon = \{ k \in \mathbb{Z}^{D} \colon 0 \leq k_i \leq N_i, \, \forall i = 1, \cdots, D\}`` for some multibandlimit vector ``N \in \mathbb{N}^{D}``. The transposed (adjoined) problem reads as

```math
\hat{h}_k \colon = \sum^{M-1}_{j = 0} f_j e^{-2 \pi \mathrm{i} \mathbf{k} \cdot \mathbf{x_j}}, \; k \in I_N
```

for given knots ``x_k \in [0, \pi ]^D, \; k = 0, \cdots, M-1``, and coefficients fj ∈ C, j ∈ I
l
M.

# Fields
* `N` - the multibandlimit of the trigonometric polynomial f.
* `M` - the number of nodes.
* `n` - the oversampling per dimension.
* `m` - the window size. Larger m means more accuracy but also more computational costs. 
* `f1` - the NFST flags.
* `f2` - the FFTW flags.
* `init_done` - indicates if the plan is initialized.
* `finalized` - indicates if the plan is finalized.
* `x` - the nodes.
* `f` - the function values.
* `fhat` - the Fourier coefficients.
* `plan`

# Constructor
    NFST{D}(N::NTuple{D,Int32},M::Int32,n::NTuple{D,Int32},m::Int32,f1::UInt32,f2::UInt32) where {D}

# Additional Constructor
    NFST(N::NTuple{D,Int32},M::Int32,n::NTuple{D,Int32},m::Int32,f1::UInt32,f2::UInt32) where {D}
    NFST(N::NTuple{D,Int32},M::Int32) where {D}

# See also
[`NFFT`](@ref)
"""
mutable struct NFST{D}
    N::NTuple{D,Int32}      # bandwidth tuple
    M::Int32                # number of nodes
    n::NTuple{D,Int32}      # oversampling per dimension
    m::Int32                # windows size
    f1::UInt32              # NFST flags
    f2::UInt32              # FFTW flags
    init_done::Bool         # bool for plan init
    finalized::Bool    # bool for finalizer
    x::Ref{Float64}         # nodes
    f::Ref{Float64}      # function values
    fhat::Ref{Float64}   # Fourier coefficients
    plan::Ref{nfst_plan}    # plan (C pointer)
    function NFST{D}(
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

# additional constructor for easy use [NFST((N,N),M) instead of NFST{2}((N,N),M)]
function NFST(N::NTuple{D,Integer}, M::Integer) where {D}
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

    # default NFST flags
    f1 = UInt32(0)

    if D > 1
        f1 = f1_default
    else
        f1 = f1_default_1d
    end

    NFST{D}(NTuple{D,Int32}(N), Int32(M), n, Int32(8), f1, f2_default)
end

function NFST(
    N::NTuple{D,Integer},
    M::Integer,
    n::NTuple{D,Integer},
    m::Integer = Int32(8),
    f1::UInt32 = (D > 1 ? f1_default : f1_default_1d),
    f2::UInt32 = f2_default,
) where {D}

    # safety checks
    if any(x -> x <= 0, N)
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

    NFST{D}(
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
    nfst_finalize_plan(P)

destroys a NFST plan structure.

# Input
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_init`](@ref)
"""
function nfst_finalize_plan(P::NFST{D}) where {D}
    if !P.init_done
        error("NFST not initialized.")
    end

    if !P.finalized
        Core.setfield!(P, :finalized, true)
        ccall(("jnfst_finalize", lib_path_nfst), Nothing, (Ref{nfst_plan},), P.plan)
    end
end

function finalize_plan(P::NFST{D}) where {D}
    return nfst_finalize_plan(P)
end

# allocate plan memory and init with D,N,M,n,m,f1,f2
@doc raw"""
    nfst_init(P)

intialises a transform plan.

# Input
* `p` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_finalize_plan`](@ref)
"""
function nfst_init(P::NFST{D}) where {D}
    # convert N and n to vectors for passing them over to C
    Nv = collect(P.N)
    n = collect(P.n)

    # call init for memory allocation
    ptr = ccall(("jnfst_alloc", lib_path_nfst), Ptr{nfst_plan}, ())

    # set pointer
    Core.setfield!(P, :plan, ptr)

    # initialize values
    ccall(
        ("jnfst_init", lib_path_nfst),
        Nothing,
        (Ref{nfst_plan}, Int32, Ref{Int32}, Int32, Ref{Int32}, Int32, UInt32, UInt32),
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
    finalizer(nfst_finalize_plan, P)
end

function init(P::NFST{D}) where {D}
    return nfst_init(P)
end

# overwrite dot notation for plan struct in order to use C memory
function Base.setproperty!(P::NFST{D}, v::Symbol, val) where {D}
    # init plan if not done [usually with setting nodes]
    if !P.init_done
        nfst_init(P)
    end

    # prevent bad stuff from happening
    if P.finalized
        error("NFST already finalized")
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
            ("jnfst_set_x", lib_path_nfst),
            Ptr{Float64},
            (Ref{nfst_plan}, Ref{Cdouble}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)

        # setting values
    elseif v == :f
        if typeof(val) != Array{Float64,1}
            error("f has to be a Float64 vector.")
        end
        if size(val)[1] != P.M
            error("f has to be a Float64 vector of size M.")
        end
        ptr = ccall(
            ("jnfst_set_f", lib_path_nfst),
            Ptr{Float64},
            (Ref{nfst_plan}, Ref{Float64}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)

        # setting Fourier coefficients
    elseif v == :fhat
        if typeof(val) != Array{Float64,1}
            error("fhat has to be a Float64 vector.")
        end
        l = prod(collect(P.N) .- 1)
        if size(val)[1] != l
            error("fhat has to be a Float64 vector of size prod(N-1).")
        end
        ptr = ccall(
            ("jnfst_set_fhat", lib_path_nfst),
            Ptr{Float64},
            (Ref{nfst_plan}, Ref{Float64}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)

        # prevent modification of NFST plan pointer
    elseif v == :plan
        @warn "You can't modify the C pointer to the NFST plan."
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
        @warn "You can't modify the NFST flags, please create an additional plan."
    elseif v == :f2
        @warn "You can't modify the FFTW flags, please create an additional plan."
        # handle other set operations the default way
    else
        Core.setfield!(P, v, val)
    end
end

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(P::NFST{D}, v::Symbol) where {D}
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
        return ccall(("nfft_get_num_threads", lib_path_nfst), Int64, ())
    elseif v == :f
        if !isdefined(P, :f)
            error("f is not set.")
        end
        ptr = Core.getfield(P, :f)
        return unsafe_wrap(Vector{Float64}, ptr, P.M)  # get function values from C memory and convert to Julia type
    elseif v == :fhat
        if !isdefined(P, :fhat)
            error("fhat is not set.")
        end
        ptr = Core.getfield(P, :fhat)
        return unsafe_wrap(Vector{Float64}, ptr, prod(collect(P.N) .- 1)) # get Fourier coefficients from C memory and convert to Julia type
    else
        return Core.getfield(P, v)
    end
end

# nfst trafo direct [call with NFST.trafo_direct outside module]
@doc raw"""
    nfst_trafo_direct(P)

computes a NFST.

# Input
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_trafo`](@ref)
"""
function nfst_trafo_direct(P::NFST{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFST already finalized")
    end

    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end

    if !isdefined(P, :x)
        error("x has not been set.")
    end

    ptr = ccall(
        ("jnfst_trafo_direct", lib_path_nfst),
        Ptr{Float64},
        (Ref{nfst_plan},),
        P.plan,
    )
    Core.setfield!(P, :f, ptr)
end

function trafo_direct(P::NFST{D}) where {D}
    return nfst_trafo_direct(P)
end

# adjoint trafo direct [call with NFST.adjoint_direct outside module]
@doc raw"""
    nfst_adjoint_direct(P)

computes an adjoint NFST.

# Input
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_adjoint`](@ref)
"""
function nfst_adjoint_direct(P::NFST{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFST already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(
        ("jnfst_adjoint_direct", lib_path_nfst),
        Ptr{Float64},
        (Ref{nfst_plan},),
        P.plan,
    )
    Core.setfield!(P, :fhat, ptr)
end

function adjoint_direct(P::NFST{D}) where {D}
    return nfst_adjoint_direct(P)
end

# nfst trafo [call with NFST.trafo outside module]
@doc raw"""
    nfst_trafo(P)

computes a NFST.

# Input
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_trafo_direct`](@ref)
"""
function nfst_trafo(P::NFST{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFST already finalized")
    end
    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(("jnfst_trafo", lib_path_nfst), Ptr{Float64}, (Ref{nfst_plan},), P.plan)
    Core.setfield!(P, :f, ptr)
end

function trafo(P::NFST{D}) where {D}
    return nfst_trafo(P)
end

# adjoint trafo [call with NFST.adjoint outside module]
@doc raw"""
    nfst_adjoint(P)

computes an adjoint NFST.

# Input
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_adjoint_direct`](@ref)
"""
function nfst_adjoint(P::NFST{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFST already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(("jnfst_adjoint", lib_path_nfst), Ptr{Float64}, (Ref{nfst_plan},), P.plan)
    Core.setfield!(P, :fhat, ptr)
end

function adjoint(P::NFST{D}) where {D}
    return nfst_adjoint(P)
end

