# dummy struct for C
mutable struct nfct_plan end

# NFCT plan struct
@doc raw"""
    NFCT{D}

A NFCT (nonequispaced fast cosine transform) plan, where D is the dimension. 

The NFCT realizes a direct and fast computation of the discrete nonequispaced cosine transform. The aim is to compute

```math
f^c(\pmb{x}_j) = \sum_{\pmb{k} \in I_{\pmb{N},\mathrm{c}}^D} \hat{f}_{\pmb{k}}^c \, \cos(2\pi \, \pmb{k} \odot \pmb{x}_j)
```

at given arbitrary knots ``\pmb{x}_j \in [0,0.5]^D, j = 1, \cdots, M``, for coefficients ``\hat{f}^{c}_{\pmb{k}} \in \mathbb{R}``, ``\pmb{k} \in I_{\pmb{N},\mathrm{c}}^D \coloneqq \left\{ \pmb{k} \in \mathbb{Z}^D: 1 \leq k_i \leq N_i - 1, \, i = 1,2,\ldots,D \right\}``, and a multibandlimit vector ``\pmb{N} \in \mathbb{N}^{D}``. Note that we define ``\cos(\pmb{k} \circ \pmb{x}) \coloneqq \prod_{i=1}^D \cos(k_i \cdot x_i)``. The transposed problem reads as

```math
\hat{h}^c_{\pmb{k}} = \sum_{j=1}^M f^c_j \, \cos(2\pi \, \pmb{k} \odot \pmb{x}_j)
```

for the frequencies ``\pmb{k} \in I_{\pmb{N},\mathrm{c}}^D`` with given coefficients ``f^c_j \in \mathbb{R}, j = 1,2,\ldots,M``.

# Fields
* `N` - the multibandlimit ``(N_1, N_2, \ldots, N_D)`` of the trigonometric polynomial ``f^s``.
* `M` - the number of nodes.
* `n` - the oversampling ``(n_1, n_2, \ldots, n_D)`` per dimension.
* `m` - the window size. A larger m results in more accuracy but also a higher computational cost. 
* `f1` - the NFCT flags.
* `f2` - the FFTW flags.
* `init_done` - indicates if the plan is initialized.
* `finalized` - indicates if the plan is finalized.
* `x` - the nodes ``x_j \in [0,0.5]^D, \, j = 1, \ldots, M``.
* `f` - the values ``f^c(\pmb{x}_j)`` for the NFCT or the coefficients ``f_j^c \in \mathbb{R}, j = 1, \ldots, M,`` for the transposed NFCT.
* `fhat` - the Fourier coefficients ``\hat{f}_{\pmb{k}}^c \in \mathbb{R}`` for the NFCT or the values ``\hat{h}_{\pmb{k}}^c, \pmb{k} \in I_{\pmb{N},\mathrm{c}}^D,`` for the adjoint NFCT.
* `plan` - plan (C pointer).

# Constructor
    NFCT{D}( N::NTuple{D,Int32}, M::Int32, n::NTuple{D,Int32}, m::Int32, f1::UInt32, f2::UInt32 ) where {D}

# Additional Constructor
    NFCT( N::NTuple{D,Int32}, M::Int32, n::NTuple{D,Int32}, m::Int32, f1::UInt32, f2::UInt32) where {D}
    NFCT( N::NTuple{D,Int32}, M::Int32) where {D}

# See also
[`NFCT`](@ref)
"""
mutable struct NFCT{D}
    N::NTuple{D,Int32}      # bandwidth tuple
    M::Int32                # number of nodes
    n::NTuple{D,Int32}      # oversampling per dimension
    m::Int32                # windows size
    f1::UInt32              # NFCT flags
    f2::UInt32              # FFTW flags
    init_done::Bool         # bool for plan init
    finalized::Bool    # bool for finalizer
    x::Ref{Float64}         # nodes
    f::Ref{Float64}      # function values
    fhat::Ref{Float64}   # Fourier coefficients
    plan::Ref{nfct_plan}    # plan (C pointer)
    function NFCT{D}(
        N::NTuple{D,Int32},
        M::Int32,
        n::NTuple{D,Int32},
        m::Int32,
        f1::UInt32,
        f2::UInt32,
    ) where {D}
        if any(x -> x <= 0, N)
            throw(DomainError(N, "argument must be a positive integer")) 
        end

        if sum(N .% 2) != 0
            throw(DomainError(N, "argument must be an even integer")) 
        end

        if M <= 0
            throw(DomainError(M, "argument must be a positive integer"))
        end

        if any(x -> x <= 0, n)
            throw(DomainError(n, "argument must be a positive integer")) 
        end

        if n <= N
            throw(DomainError(n, "argument must fulfil n_i > N_i")) 
        end

        if sum(n .% 2) != 0
            throw(DomainError(n, "argument must be an even integer")) 
        end

        if m <= 0
            throw(DomainError(m, "argument must be a positive integer")) 
        end

        new(N, M, n, m, f1, f2, false, false)
    end
end

# additional constructor for easy use [NFCT((N,N),M) instead of NFCT{2}((N,N),M)]
function NFCT(N::NTuple{D,Integer}, M::Integer) where {D}
    if any(x -> x <= 0, N)
            throw(DomainError(N, "argument must be a positive integer")) 
    end

    # convert N to vector for passing it over to C
    Nv = collect(N)

    # default oversampling
    n = Array{Int32}(2 .^ (ceil.(log.(Nv) / log(2)) .+ 1))
    n = NTuple{D,Int32}(n)

    # default NFCT flags
    f1 = UInt32(0)

    if D > 1
        f1 = f1_default
    else
        f1 = f1_default_1d
    end

    NFCT{D}(NTuple{D,Int32}(N), Int32(M), n, Int32(8), f1, f2_default)
end

function NFCT(
    N::NTuple{D,Integer},
    M::Integer,
    n::NTuple{D,Integer},
    m::Integer = Int32(8),
    f1::UInt32 = (D > 1 ? f1_default : f1_default_1d),
    f2::UInt32 = f2_default,
) where {D}
    NFCT{D}(
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
    nfct_finalize_plan(P::NFCT{D})

destroys a NFCT plan structure.

# Input
* `P` - a NFCT plan structure.

# See also
[`NFCT{D}`](@ref), [`nfct_init`](@ref)
"""
function nfct_finalize_plan(P::NFCT{D}) where {D}
    if !P.init_done
        error("NFCT not initialized.")
    end

    if !P.finalized
        Core.setfield!(P, :finalized, true)
        ccall(("jnfct_finalize", lib_path_nfct), Nothing, (Ref{nfct_plan},), P.plan)
    end
end

function finalize_plan(P::NFCT{D}) where {D}
    return nfct_finalize_plan(P)
end

# allocate plan memory and init with D,N,M,n,m,f1,f2
@doc raw"""
    nfct_init(P)

intialises the NFCT plan in C. This function does not have to be called by the user.

# Input
* `P` - a NFCT plan structure.

# See also
[`NFCT{D}`](@ref), [`nfct_finalize_plan`](@ref)
"""
function nfct_init(P::NFCT{D}) where {D}
    # convert N and n to vectors for passing them over to C
    Nv = collect(P.N)
    n = collect(P.n)

    # call init for memory allocation
    ptr = ccall(("jnfct_alloc", lib_path_nfct), Ptr{nfct_plan}, ())

    # set pointer
    Core.setfield!(P, :plan, ptr)

    # initialize values
    ccall(
        ("jnfct_init", lib_path_nfct),
        Nothing,
        (Ref{nfct_plan}, Int32, Ref{Int32}, Int32, Ref{Int32}, Int32, UInt32, UInt32),
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
    finalizer(nfct_finalize_plan, P)
end

function init(P::NFCT{D}) where {D}
    return nfct_init(P)
end

# overwrite dot notation for plan struct in order to use C memory
function Base.setproperty!(P::NFCT{D}, v::Symbol, val) where {D}
    # init plan if not done [usually with setting nodes]
    if !P.init_done
        nfct_init(P)
    end

    # prevent bad stuff from happening
    if P.finalized
        error("NFCT already finalized")
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
            ("jnfct_set_x", lib_path_nfct),
            Ptr{Float64},
            (Ref{nfct_plan}, Ref{Cdouble}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)

        # setting values
    elseif v == :f
        if !isa( val, Array{<:Real,1} )
            error("f has to be a vector of real numbers.")
        end
        if size(val)[1] != P.M
            error("f has to be a Float64 vector of size M.")
        end
        f_real = convert(Vector{Float64},val)
        ptr = ccall(
            ("jnfct_set_f", lib_path_nfct),
            Ptr{Float64},
            (Ref{nfct_plan}, Ref{Float64}),
            P.plan,
            f_real,
        )
        Core.setfield!(P, v, ptr)

        # setting Fourier coefficients
    elseif v == :fhat
        if !isa( val, Array{<:Real,1} )
            error("fhat has to be a vector of real numbers.")
        end
        l = prod(P.N)
        if size(val)[1] != l
            error("fhat has to be a Float64 vector of size prod(N).")
        end
        fhat_real = convert(Vector{Float64},val)
        ptr = ccall(
            ("jnfct_set_fhat", lib_path_nfct),
            Ptr{Float64},
            (Ref{nfct_plan}, Ref{Float64}),
            P.plan,
            fhat_real,
        )
        Core.setfield!(P, v, ptr)
        # prevent modification of NFCT plan pointer
    elseif v == :plan
        @warn "You can't modify the C pointer to the NFCT plan."
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
        @warn "You can't modify the NFCT flags, please create an additional plan."
    elseif v == :f2
        @warn "You can't modify the FFTW flags, please create an additional plan."
        # handle other set operations the default way
    else
        Core.setfield!(P, v, val)
    end
end

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(P::NFCT{D}, v::Symbol) where {D}
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
        return ccall(("nfft_get_num_threads", lib_path_nfct), Int64, ())
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
        return unsafe_wrap(Vector{Float64}, ptr, prod(P.N)) # get Fourier coefficients from C memory and convert to Julia type
    else
        return Core.getfield(P, v)
    end
end

# nfct trafo direct [call with NFCT.trafo_direct outside module]
@doc raw"""
    nfct_trafo_direct(P)

computes the NDCT via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}}^c \in \mathbb{R}, \pmb{k} \in I_{\pmb{N},\mathrm{c}}^D,`` in `P.fhat`.

# Input
* `P` - a NFCT plan structure.

# See also
[`NFCT{D}`](@ref), [`nfct_trafo`](@ref)
"""
function nfct_trafo_direct(P::NFCT{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFCT already finalized")
    end

    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end

    if !isdefined(P, :x)
        error("x has not been set.")
    end

    ptr = ccall(
        ("jnfct_trafo_direct", lib_path_nfct),
        Ptr{Float64},
        (Ref{nfct_plan},),
        P.plan,
    )
    Core.setfield!(P, :f, ptr)
end

function trafo_direct(P::NFCT{D}) where {D}
    return nfct_trafo_direct(P)
end

# adjoint trafo direct [call with NFCT.adjoint_direct outside module]
@doc raw"""
    nfct_transposed_direct(P)

computes the transposed NDCT via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j^c \in \mathbb{R}, j =1,2,\dots,M,`` in `P.f`.

# Input
* `P` - a NFCT plan structure.

# See also
[`NFCT{D}`](@ref), [`nfct_transposed`](@ref)
"""
function nfct_transposed_direct(P::NFCT{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFCT already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(
        ("jnfct_adjoint_direct", lib_path_nfct),
        Ptr{Float64},
        (Ref{nfct_plan},),
        P.plan,
    )
    Core.setfield!(P, :fhat, ptr)
end

function nfct_adjoint_direct(P::NFCT{D}) where {D}
    return nfct_transposed_direct(P)
end

function adjoint_direct(P::NFCT{D}) where {D}
    return nfct_adjoint_direct(P)
end

# nfct trafo [call with NFCT.trafo outside module]
@doc raw"""
    nfct_trafo(P)

computes the NDCT via the fast NFCT algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}}^c \in \mathbb{R}, \pmb{k} \in I_{\pmb{N},\mathrm{c}}^D,`` in `P.fhat`.

# Input
* `P` - a NFCT plan structure.

# See also
[`NFCT{D}`](@ref), [`nfct_trafo_direct`](@ref)
"""
function nfct_trafo(P::NFCT{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFCT already finalized")
    end
    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(("jnfct_trafo", lib_path_nfct), Ptr{Float64}, (Ref{nfct_plan},), P.plan)
    Core.setfield!(P, :f, ptr)
end

function trafo(P::NFCT{D}) where {D}
    return nfct_trafo(P)
end

# adjoint trafo [call with NFCT.adjoint outside module]
@doc raw"""
    nfct_transposed(P)

computes the transposed NDCT via the fast transposed NFCT algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j^c \in \mathbb{R}, j =1,2,\dots,M,`` in `P.f`.

# Input
* `P` - a NFCT plan structure.

# See also
[`NFCT{D}`](@ref), [`nfct_transposed_direct`](@ref)
"""
function nfct_transposed(P::NFCT{D}) where {D}
    # prevent bad stuff from happening
    if P.finalized
        error("NFCT already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(("jnfct_adjoint", lib_path_nfct), Ptr{Float64}, (Ref{nfct_plan},), P.plan)
    Core.setfield!(P, :fhat, ptr)
end

function nfct_adjoint(P::NFCT{D}) where {D}
    return nfct_transposed(P)
end

function adjoint(P::NFCT{D}) where {D}
    return nfct_adjoint(P)
end

@doc raw"""
    nfct_get_LinearMap(N::Vector{Int}, X::Array{Float64}; n, m::Integer, f1::UInt32, f2::UInt32)::LinearMap

gives an linear map which computes the NFCT.

# Input
* `N` - the multibandlimit ``(N_1, N_2, \ldots, N_D)`` of the trigonometric polynomial ``f``.
* `X` - the nodes X.
* `n` - the oversampling ``(n_1, n_2, \ldots, n_D)`` per dimension.
* `m` - the window size. A larger m results in more accuracy but also a higher computational cost. 
* `f1` - the NFCT flags.
* `f2` - the FFTW flags.

# See also
[`NFCT{D}`](@ref)
"""
function nfct_get_LinearMap(
    N::Vector{Int},
    X::Array{Float64};
    n = undef,
    m::Integer = 5,
    f1::UInt32 = (size(X, 1) > 1 ? f1_default : f1_default_1d),
    f2::UInt32 = f2_default,
)::LinearMap
    if size(X, 1) == 1
        X = vec(X)
        d = 1
        M = length(X)
    else
        (d, M) = size(X)
    end

    if N == []
        return LinearMap{Float64}(fhat -> fill(fhat[1], M), f -> [sum(f)], M, 1)
    end

    N = Tuple(N)
    if n == undef
        n = Tuple(2 * collect(N))
    end

    plan = NFCT(N, M, n, m, f1, f2)
    plan.x = X

    function trafo(fhat::Vector{Float64})::Vector{Float64}
        plan.fhat = fhat
        nfct_trafo(plan)
        return plan.f
    end

    function adjoint(f::Vector{Float64})::Vector{Float64}
        plan.f = f
        nfct_adjoint(plan)
        return plan.fhat
    end

    N = prod(N)
    return LinearMap{Float64}(trafo, adjoint, M, N)
end

@doc raw"""
    nfct_get_coefficient_vector(fhat::Array{Float64})::Vector{Float64}

reshapes an coefficient array to an vector for multiplication with the linear map of the NFCT.

# Input
* `fhat` - the Fourier coefficients.

# See also
[`NFCT{D}`](@ref), [`nfct_get_LinearMap`](@ref)
"""
function nfct_get_coefficient_vector(fhat::Array{Float64})::Vector{Float64}
    N = size(fhat)
    return vec(permutedims(fhat,length(N):-1:1))
end

@doc raw"""
    nfct_get_coefficient_array(fhat::Vector{Float64},P::NFCT{D})::Array{Float64} where {D}

reshapes an coefficient vector returned from a linear map of the NFCT to an array.

# Input
* `fhat` - the Fourier coefficients.
* `P` - a NFCT plan structure.

# See also
[`NFCT{D}`](@ref), [`nfct_get_LinearMap`](@ref)
"""
function nfct_get_coefficient_array(fhat::Vector{Float64},P::NFCT{D})::Array{Float64} where {D}
    return permutedims(reshape(fhat,reverse(P.N)),length(P.N):-1:1)
end

@doc raw"""
    nfct_get_coefficient_array(fhat::Vector{Float64},N::Vector{Int64})::Array{Float64}

reshapes an coefficient vector returned from a linear map of the NFCT to an array.

# Input
* `fhat` - the Fourier coefficients.
* `N` - the multibandlimit ``(N_1, N_2, \ldots, N_D)`` of the trigonometric polynomial ``f``.

# See also
[`NFCT{D}`](@ref), [`nfct_get_LinearMap`](@ref)
"""
function nfct_get_coefficient_array(fhat::Vector{Float64},N::Vector{Int64})::Array{Float64}
    N = Tuple(N)
    return permutedims(reshape(fhat,reverse(N)),length(N):-1:1)
end

@doc raw"""
    *(plan::NFCT{D}, fhat::Array{Float64})::Vector{Float64}

This function defines the multiplication of an NFCT plan with an coefficient array.
"""
function Base.:*(plan::NFCT{D}, fhat::Array{Float64})::Vector{Float64} where {D}
    if !isdefined(plan,:x)
        error("x is not set.")
    end
    plan.fhat = nfct_get_coefficient_vector(fhat)
    nfct_trafo(plan)
    return plan.f
end

struct Adjoint_NFCT{D}
    plan::NFCT{D}
end

Adjoint_NFCT{D}(plan::NFCT{D}) where {D} = plan

@doc raw"""
    adjoint(plan::NFCT{D})::Adjoint_NFCT{D}

This function defines the adjoint operator for an NFCT.
"""
function Base.adjoint(plan::NFCT{D})::Adjoint_NFCT{D} where {D}
    return Adjoint_NFCT(plan)
end

@doc raw"""
    *(plan::Adjoint_NFCT{D}, f::Vector{Float64})::Array{Float64}

This function defines the multiplication of an adjoint NFCT plan with an vector of function values.
"""
function Base.:*(plan::Adjoint_NFCT{D}, f::Vector{Float64})::Array{Float64} where {D}
    if !isdefined(plan.plan,:x)
        error("x is not set.")
    end
    plan.plan.f = f
    nfct_adjoint(plan.plan)
    return nfct_get_coefficient_array(plan.plan.fhat, plan.plan)
end
