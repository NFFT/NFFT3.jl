# dummy struct for C
mutable struct nfst_plan end

# NFST plan struct
@doc raw"""
    NFST{D}

A NFST (nonequispaced fast sine transform) plan, where D is the dimension. 

The NFST realizes a direct and fast computation of the discrete nonequispaced sine transform. The aim is to compute

```math
f^s(\pmb{x}_j) = \sum_{\pmb{k} \in I_{\pmb{N},\mathrm{s}}^D} \hat{f}_{\pmb{k}}^s \, \sin(2\pi \, \pmb{k} \odot \pmb{x}_j)
```

at given arbitrary knots ``\pmb{x}_j \in [0,0.5]^D, j = 1, \cdots, M``, for coefficients ``\hat{f}^{s}_{\pmb{k}} \in \mathbb{R}``, ``\pmb{k} \in I_{\pmb{N},\mathrm{s}}^D \coloneqq \left\{ \pmb{k} \in \mathbb{Z}^D: 1 \leq k_i \leq N_i - 1, \, i = 1,2,\ldots,D \right\}``, and a multibandlimit vector ``\pmb{N} \in \mathbb{N}^{D}``. Note that we define ``\sin(\pmb{k} \circ \pmb{x}) \coloneqq \prod_{i=1}^D \sin(k_i \cdot x_i)``. The transposed problem reads as

```math
\hat{h}^s_{\pmb{k}} = \sum_{j=1}^M f^s_j \, \sin(2\pi \, \pmb{k} \odot \pmb{x}_j)
```

for the frequencies ``\pmb{k} \in I_{\pmb{N},\mathrm{s}}^D`` with given coefficients ``f^s_j \in \mathbb{R}, j = 1,2,\ldots,M``.

# Fields
* `N` - the multibandlimit ``(N_1, N_2, \ldots, N_D)`` of the trigonometric polynomial ``f^s``.
* `M` - the number of nodes.
* `n` - the oversampling ``(n_1, n_2, \ldots, n_D)`` per dimension.
* `m` - the window size. A larger m results in more accuracy but also a higher computational cost. 
* `f1` - the NFST flags.
* `f2` - the FFTW flags.
* `init_done` - indicates if the plan is initialized.
* `finalized` - indicates if the plan is finalized.
* `x` - the nodes ``x_j \in [0,0.5]^D, \, j = 1, \ldots, M``.
* `f` - the values ``f^s(\pmb{x}_j)`` for the NFST or the coefficients ``f_j^s \in \mathbb{R}, j = 1, \ldots, M,`` for the transposed NFST.
* `fhat` - the Fourier coefficients ``\hat{f}_{\pmb{k}}^s \in \mathbb{R}`` for the NFST or the values ``\hat{h}_{\pmb{k}}^s, \pmb{k} \in I_{\pmb{N},\mathrm{s}}^D,`` for the adjoint NFST.
* `plan` - plan (C pointer).

# Constructor
    NFST{D}( N::NTuple{D,Int32}, M::Int32, n::NTuple{D,Int32}, m::Int32, f1::UInt32, f2::UInt32 ) where {D}

# Additional Constructor
    NFST( N::NTuple{D,Int32}, M::Int32, n::NTuple{D,Int32}, m::Int32, f1::UInt32, f2::UInt32) where {D}
    NFST( N::NTuple{D,Int32}, M::Int32) where {D}

# See also
[`NFST`](@ref)
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
        if any(x -> x <= 0, N)
            throw(DomainError(N, "argument must be a positive integer")) 
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

# additional constructor for easy use [NFST((N,N),M) instead of NFST{2}((N,N),M)]
function NFST(N::NTuple{D,Integer}, M::Integer) where {D}
    if any(x -> x <= 0, N)
        throw(DomainError(N, "argument must be a positive integer")) 
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

intialises the NFST plan in C. This function does not have to be called by the user.

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
        if !isa( val, Array{<:Real,1} )
            error("f has to be a vector of real numbers.")
        end
        if size(val)[1] != P.M
            error("f has to be a Float64 vector of size M.")
        end
        f_real = convert(Vector{Float64},val)
        ptr = ccall(
            ("jnfst_set_f", lib_path_nfst),
            Ptr{Float64},
            (Ref{nfst_plan}, Ref{Float64}),
            P.plan,
            f_real,
        )
        Core.setfield!(P, v, ptr)
        # setting Fourier coefficients
    elseif v == :fhat
        if !isa( val, Array{<:Real,1} )
            error("fhat has to be a vector of real numbers.")
        end
        l = prod(collect(P.N) .- 1)
        if size(val)[1] != l
            error("fhat has to be a Float64 vector of size prod(N-1).")
        end
        fhat_real = convert(Vector{Float64},val)
        ptr = ccall(
            ("jnfst_set_fhat", lib_path_nfst),
            Ptr{Float64},
            (Ref{nfst_plan}, Ref{Float64}),
            P.plan,
            fhat_real,
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

@doc raw"""
    nfst_trafo_direct(P)

computes the NDST via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}}^s \in \mathbb{R}, \pmb{k} \in I_{\pmb{N},s}^D,`` in `P.fhat`.

# Input
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_trafo`](@ref)
"""
function nfst_trafo_direct(P::NFST{D}) where {D}
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

@doc raw"""
    nfst_transposed_direct(P)

computes the transposed NDST via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j^s \in \mathbb{R}, j =1,2,\dots,M,`` in `P.f`.

# Input
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_transposed`](@ref)
"""
function nfst_transposed_direct(P::NFST{D}) where {D}
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

function nfst_adjoint_direct(P::NFST{D}) where {D}
    return nfst_transposed_direct(P)
end

function adjoint_direct(P::NFST{D}) where {D}
    return nfst_adjoint_direct(P)
end

@doc raw"""
    nfst_trafo(P)

computes the NDST via the fast NFST algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}}^s \in \mathbb{R}, \pmb{k} \in I_{\pmb{N},s}^D,`` in `P.fhat`.

# Input
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_trafo_direct`](@ref)
"""
function nfst_trafo(P::NFST{D}) where {D}
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

@doc raw"""
    nfst_transposed(P)

computes the transposed NDST via the fast transposed NFST algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j^s \in \mathbb{R}, j =1,2,\dots,M,`` in `P.f`.

# Input
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_transposed_direct`](@ref)
"""
function nfst_transposed(P::NFST{D}) where {D}
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

function nfst_adjoint(P::NFST{D}) where {D}
    return nfst_transposed(P)
end

function adjoint(P::NFST{D}) where {D}
    return nfst_adjoint(P)
end

@doc raw"""
    nfst_get_LinearMap(N::Vector{Int}, X::Array{Float64}; n, m::Integer, f1::UInt32, f2::UInt32)::LinearMap

gives an linear map which computes the NFST.

# Input
* `N` - the multibandlimit ``(N_1, N_2, \ldots, N_D)`` of the trigonometric polynomial ``f``.
* `X` - the nodes X.
* `n` - the oversampling ``(n_1, n_2, \ldots, n_D)`` per dimension.
* `m` - the window size. A larger m results in more accuracy but also a higher computational cost. 
* `f1` - the NFST flags.
* `f2` - the FFTW flags.

# See also
[`NFST{D}`](@ref)
"""
function nfst_get_LinearMap(
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

    plan = NFST(N, M, n, m, f1, f2)
    plan.x = X

    function trafo(fhat::Vector{Float64})::Vector{Float64}
        plan.fhat = fhat
        nfst_trafo(plan)
        return plan.f
    end

    function adjoint(f::Vector{Float64})::Vector{Float64}
        plan.f = f
        nfst_adjoint(plan)
        return plan.fhat
    end

    N = prod(N)
    return LinearMap{Float64}(trafo, adjoint, M, N)
end

@doc raw"""
    nfst_get_coefficient_vector(fhat::Array{Float64})::Vector{Float64}

reshapes an coefficient array to an vector for multiplication with the linear map of the NFST.

# Input
* `fhat` - the Fourier coefficients.

# See also
[`NFST{D}`](@ref), [`nfst_get_LinearMap`](@ref)
"""
function nfst_get_coefficient_vector(fhat::Array{Float64})::Vector{Float64}
    N = size(fhat)
    return vec(permutedims(fhat,length(N):-1:1))
end

@doc raw"""
    nfst_get_coefficient_array(fhat::Vector{Float64},P::NFST{D})::Array{Float64} where {D}

reshapes an coefficient vector returned from a linear map of the NFST to an array.

# Input
* `fhat` - the Fourier coefficients.
* `P` - a NFST plan structure.

# See also
[`NFST{D}`](@ref), [`nfst_get_LinearMap`](@ref)
"""
function nfst_get_coefficient_array(fhat::Vector{Float64},P::NFST{D})::Array{Float64} where {D}
    N = P.N .- 1
    return permutedims(reshape(fhat,reverse(P.N)),length(P.N):-1:1)
end

@doc raw"""
    nfst_get_coefficient_array(fhat::Vector{Float64},N::Vector{Int64})::Array{Float64}

reshapes an coefficient vector returned from a linear map of the NFST to an array.

# Input
* `fhat` - the Fourier coefficients.
* `N` - the multibandlimit ``(N_1, N_2, \ldots, N_D)`` of the trigonometric polynomial ``f``.

# See also
[`NFST{D}`](@ref), [`nfst_get_LinearMap`](@ref)
"""
function nfst_get_coefficient_array(fhat::Vector{Float64},N::Vector{Int64})::Array{Float64}
    N = Tuple(N.-1)
    return permutedims(reshape(fhat,reverse(N)),length(N):-1:1)
end

@doc raw"""
    *(plan::NFST{D}, fhat::Array{Float64})::Vector{Float64}

This function defines the multiplication of an NFST plan with an coefficient array.
"""
function Base.:*(plan::NFST{D}, fhat::Array{Float64})::Vector{Float64} where {D}
    if !isdefined(plan,:x)
        error("x is not set.")
    end
    plan.fhat = nfst_get_coefficient_vector(fhat)
    nfst_trafo(plan)
    return plan.f
end

struct Adjoint_NFST{D}
    plan::NFST{D}
end

Adjoint_NFST{D}(plan::NFST{D}) where {D} = plan

@doc raw"""
    adjoint(plan::NFST{D})::Adjoint_NFST{D}

This function defines the adjoint operator for an NFST.
"""
function Base.adjoint(plan::NFST{D})::Adjoint_NFST{D} where {D}
    return Adjoint_NFST(plan)
end

@doc raw"""
    *(plan::Adjoint_NFST{D}, f::Vector{Float64})::Array{Float64}

This function defines the multiplication of an adjoint NFST plan with an vector of function values.
"""
function Base.:*(plan::Adjoint_NFST{D}, f::Vector{Float64})::Array{Float64} where {D}
    if !isdefined(plan.plan,:x)
        error("x is not set.")
    end
    plan.plan.f = f
    nfst_adjoint(plan.plan)
    return nfst_get_coefficient_array(plan.plan.fhat, plan.plan)
end
