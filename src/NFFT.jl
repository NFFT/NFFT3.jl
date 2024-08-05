# dummy struct for C
mutable struct nfft_plan end

# NFFT plan struct
@doc raw"""
    NFFT{D}

A NFFT (nonequispaced fast Fourier transform) plan, where D is the dimension. 

Considering a D-dimensional trigonometric polynomial

```math
f \colon \mathbb{T}^D \to \mathbb{C}, \; f(\pmb{x}) \colon = \sum_{\pmb{k} \in I_{\pmb{N}}^D} \hat{f}_{\pmb{k}} \, \mathrm{e}^{-2 \pi \mathrm{i} \, \pmb{k} \cdot \pmb{x}}
```

with an index set ``I_{\pmb{N}}^D \coloneqq \left\{ \pmb{k} \in \mathbb{Z}^D: - \frac{N_i}{2} \leq k_i \leq \frac{N_i}{2} - 1, \, i = 1,2,\ldots,D \right\}`` where ``\pmb{N} \in (2\mathbb{N})^{D}`` is the multibandlimit. 
The NDFT (nonequispaced discrete Fourier transform) is its evaluation at ``M \in \mathbb{N}`` arbitrary points ``\pmb{x}_j \in [-0.5,0.5)^D`` for ``j = 1, \ldots, M``,

```math
f(\pmb{x}_j) \colon = \sum_{\pmb{k} \in I^D_{\pmb{N}}} \hat{f}_{\pmb{k}} \, \mathrm{e}^{-2 \pi \mathrm{i} \, \pmb{k} \cdot \pmb{x}_j}
```

with given coefficients ``\hat{f}_{\pmb{k}} \in \mathbb{C}``. The NFFT is an algorithm for the fast evaluation of the NDFT and the adjoint problem, the fast evaluation of the adjoint NDFT

```math
\hat{h}_{\pmb{k}} \coloneqq \sum^{M}_{j = 1} f_j \, \mathrm{e}^{-2 \pi \mathrm{i} \, \pmb{k} \cdot \pmb{x}_j}, \, \pmb{k} \in I_{\pmb{N}}^D,
```

for given coefficients ``f_j \in \mathbb{C}, j =1,2,\ldots,M``. Note that in general, the adjoint NDFT is not the inverse transform of the NDFT.

# Fields
* `N` - the multibandlimit ``(N_1, N_2, \ldots, N_D)`` of the trigonometric polynomial ``f``.
* `M` - the number of nodes.
* `n` - the oversampling ``(n_1, n_2, \ldots, n_D)`` per dimension.
* `m` - the window size. A larger m results in more accuracy but also a higher computational cost. 
* `f1` - the NFFT flags.
* `f2` - the FFTW flags.
* `init_done` - indicates if the plan is initialized.
* `finalized` - indicates if the plan is finalized.
* `x` - the nodes ``\pmb{x}_j \in [-0.5,0.5)^D, j = 1, \ldots, M``.
* `f` - the values ``f(\pmb{x}_j)`` for the NFFT or the coefficients ``f_j \in \mathbb{C}, j = 1, \ldots, M,`` for the adjoint NFFT.
* `fhat` - the Fourier coefficients ``\hat{f}_{\pmb{k}} \in \mathbb{C}`` for the NFFT or the values ``\hat{h}_{\pmb{k}}, \pmb{k} \in I_{\pmb{N}}^D,`` for the adjoint NFFT.
* `plan` - plan (C pointer).

# Constructor
    NFFT{D}( N::NTuple{D,Int32}, M::Int32, n::NTuple{D,Int32}, m::Int32, f1::UInt32, f2::UInt32 ) where D

# Additional Constructor
    NFFT( N::NTuple{D,Int32}, M::Int32, n::NTuple{D,Int32}, m::Int32, f1::UInt32, f2::UInt32 ) where {D}
    NFFT( N::NTuple{D,Int32}, M::Int32 ) where {D}
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

# additional constructor for easy use [NFFT((N,N),M) instead of NFFT{2}((N,N),M)]
function NFFT(
    N::NTuple{D,Integer}, 
    M::Integer;
    m::Integer = Int32(default_window_cut_off),
    f1::UInt32 = (D > 1 ? f1_default : f1_default_1d),
    f2::UInt32 = f2_default,
) where {D}
    if any(x -> x <= 0, N)
        throw(DomainError(N, "argument must be a positive integer")) 
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

intialises the NFFT plan in C. This function does not have to be called by the user.

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
        if !isa( val, Array{<:Number,1} )
            error("f has to be a vector of numbers.")
        end
        if size(val)[1] != P.M
            error("f has to be a ComplexFloat64 vector of size M.")
        end
        f_complex = convert(Vector{ComplexF64},val)
        ptr = ccall(
            ("jnfft_set_f", lib_path_nfft),
            Ptr{ComplexF64},
            (Ref{nfft_plan}, Ref{ComplexF64}),
            P.plan,
            f_complex,
        )
        Core.setfield!(P, v, ptr)
        # setting Fourier coefficients
    elseif v == :fhat
        if !isa( val, Array{<:Number,1} )
            error("fhat has to be a vector of numbers.")
        end
        l = prod(P.N)
        if size(val)[1] != l
            error("fhat has to be a ComplexFloat64 vector of size prod(N).")
        end
        fhat_complex = convert(Vector{ComplexF64},val)
        ptr = ccall(
            ("jnfft_set_fhat", lib_path_nfft),
            Ptr{ComplexF64},
            (Ref{nfft_plan}, Ref{ComplexF64}),
            P.plan,
            fhat_complex,
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

computes the NDFT via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}} \in \mathbb{C}, \pmb{k} \in I_{\pmb{N}}^D,`` in `P.fhat`.

# Input
* `P` - a NFFT plan structure.

# See also
[`NFFT{D}`](@ref), [`nfft_trafo`](@ref)
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

computes the adjoint NDFT via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j \in \mathbb{C}, j =1,2,\dots,M,`` in `P.f`.

# Input
* `P` - a NFFT plan structure.

# See also
[`NFFT{D}`](@ref), [`nfft_adjoint`](@ref)
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

computes the NDFT via the fast NFFT algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}} \in \mathbb{C}, \pmb{k} \in I_{\pmb{N}}^D,`` in `P.fhat`.

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

computes the adjoint NDFT via the fast adjoint NFFT algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j \in \mathbb{C}, j =1,2,\dots,M,`` in `P.f`.

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

function nfft_get_LinearMap(
    bandwidths::Vector{Int},
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

    if bandwidths == []
        return LinearMap{ComplexF64}(fhat -> fill(fhat[1], M), f -> [sum(f)], M, 1)
    end

    N = Tuple(bandwidths)
    if n == undef
        n = Tuple(2 * collect(N))
    end

    plan = NFFT(N, M, n, m, f1, f2)
    plan.x = X

    function trafo(fhat::Vector{ComplexF64})::Vector{ComplexF64}
        plan.fhat = fhat
        nfft_trafo(plan)
        return plan.f
    end

    function adjoint(f::Vector{ComplexF64})::Vector{ComplexF64}
        plan.f = f
        nfft_adjoint(plan)
        return plan.fhat
    end

    N = prod(bandwidths)
    return LinearMap{ComplexF64}(trafo, adjoint, M, N)
end

function nfft_get_coefficient_vector(fhat::Array{ComplexF64})::Vector{ComplexF64}
    N = size(fhat)
    return vec(permutedims(fhat,length(N):-1:1))
end

function nfft_get_coefficient_array(fhat::Vector{ComplexF64},P::NFFT{D})::Array{ComplexF64} where {D}
    return permutedims(reshape(fhat,reverse(N)),length(N):-1:1)
end

function nfft_get_coefficient_array(fhat::Vector{ComplexF64},N::Vector{Int64})::Array{ComplexF64}
    N = Tuple(N)
    return permutedims(reshape(fhat,reverse(N)),length(N):-1:1)
end

function Base.:*(plan::NFFT{D}, fhat::Array{ComplexF64})::Vector{ComplexF64} where {D}
    if !isdefined(plan,:x)
        error("x is not set.")
    end
    plan.fhat = nfft_get_coefficient_vector(fhat)
    nfft_trafo(plan)
    return plan.f
end

struct Adjoint_NFFT{D}
    plan::NFFT{D}
end

Adjoint_NFFT{D}(plan::NFFT{D}) where {D} = plan

function Base.adjoint(plan::NFFT{D})::Adjoint_NFFT{D} where {D}
    return Adjoint_NFFT(plan)
end

function Base.:*(plan::Adjoint_NFFT{D}, f::Vector{ComplexF64})::Array{ComplexF64} where {D}
    if !isdefined(plan.plan,:x)
        error("x is not set.")
    end
    plan.plan.f = f
    nfft_adjoint(plan.plan)
    return nfft_get_coefficient_array(plan.plan.fhat, plan.plan)
end
