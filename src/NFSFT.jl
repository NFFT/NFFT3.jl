# dummy struct for C
mutable struct nfsft_plan end

@doc raw"""
    NFSFT

A NFSFT (nonequispaced fast spherical Fourier transform) plan.

The NFSFT realizes a direct and fast computation of the discrete nonequispaced spherical Fourier transform. The aim is to compute

```math
f(\pmb{x}_j) = \sum_{k=0}^N \sum_{n=-k}^k \hat{f}_k^n \, Y_k^n(\pmb{x}_j)
```

at given arbitrary nodes ``\pmb{x}_j`` on the unit sphere, ``j = 1, \ldots, M``, for spherical Fourier coefficients ``\hat{f}_k^n \in \mathbb{C}``, where ``Y_k^n`` are the spherical harmonics of degree ``k`` and order ``n``. The transposed problem reads as

```math
\hat{h}_k^n = \sum_{j=1}^M f_j \, \overline{Y_k^n(\pmb{x}_j)}
```

for given coefficients ``f_j \in \mathbb{C}, j = 1, \ldots, M``.

A node ``\pmb{x} = (x_0, x_1)`` encodes spherical coordinates via ``x_0 = \varphi/(2\pi) \in [-0.5, 0.5)`` (azimuth) and ``x_1 = \vartheta/\pi \in [0, 1)`` (polar angle).

# Fields
* `N` - the bandwidth (maximum degree ``k``).
* `M` - the number of nodes.
* `f1` - the NFSFT flags.
* `f2` - the NFFT flags used internally.
* `nfft_cutoff` - the NFFT window cutoff used internally.
* `init_done` - indicates if the plan is initialized.
* `finalized` - indicates if the plan is finalized.
* `x` - the nodes ``\pmb{x}_j \in [-0.5, 0.5) \times [0, 1), \, j = 1, \ldots, M``.
* `f` - the values ``f(\pmb{x}_j)`` for the NFSFT or the coefficients ``f_j`` for the adjoint NFSFT.
* `fhat` - the spherical Fourier coefficients ``\hat{f}_k^n`` for the NFSFT, or the values ``\hat{h}_k^n`` for the adjoint NFSFT.
* `plan` - plan (C pointer).

# Constructor
    NFSFT(N::Integer, M::Integer, f1::UInt32=nfsft_default, f2::UInt32=nfsft_nfft_default, nfft_cutoff::Integer=nfsft_default_nfft_cut_off)

# See also
[`FSFT`](@ref) for the equispaced variant.
"""
mutable struct NFSFT
    N::Int32                # bandwidth
    M::Int32                # number of nodes
    f1::UInt32               # NFSFT flags
    f2::UInt32               # NFFT flags used internally
    nfft_cutoff::Int32
    init_done::Bool
    finalized::Bool
    fhat_corrupted_by_fast_trafo::Bool
    x::Ref{Float64}
    f::Ref{ComplexF64}
    fhat::Ref{ComplexF64}
    plan::Ref{nfsft_plan}
    function NFSFT(
        N::Integer,
        M::Integer,
        f1::UInt32 = nfsft_default,
        f2::UInt32 = nfsft_nfft_default,
        nfft_cutoff::Integer = nfsft_default_nfft_cut_off,
    )
        if N <= 0
            throw(DomainError(N, "argument must be a positive integer"))
        end
        if M <= 0
            throw(DomainError(M, "argument must be a positive integer"))
        end
        if (f1 & NFSFT_EQUISPACED) != 0
            error("Cannot use NFSFT_EQUISPACED flag. Use FSFT instead")
        end
        new(Int32(N), Int32(M), f1, f2, Int32(nfft_cutoff), false, false, false)
    end
end

@doc raw"""
    nfsft_finalize_plan(P::NFSFT)

destroys a NFSFT plan structure.

# See also
[`NFSFT`](@ref), [`nfsft_init`](@ref)
"""
function nfsft_finalize_plan(P::NFSFT)
    if !P.init_done
        error("NFSFT not initialized.")
    end
    if !P.finalized
        Core.setfield!(P, :finalized, true)
        ccall(("jnfsft_finalize", lib_path_nfsft), Nothing, (Ref{nfsft_plan},), P.plan)
    end
end

function finalize_plan(P::NFSFT)
    return nfsft_finalize_plan(P)
end

@doc raw"""
    nfsft_init(P::NFSFT)

initializes the NFSFT plan in C. This function does not have to be called by the user.

# See also
[`NFSFT`](@ref), [`nfsft_finalize_plan`](@ref)
"""
function nfsft_init(P::NFSFT)
    ptr = ccall(("jnfsft_alloc", lib_path_nfsft), Ptr{nfsft_plan}, ())
    Core.setfield!(P, :plan, ptr)

    ccall(
        ("jnfsft_init", lib_path_nfsft),
        Nothing,
        (Ref{nfsft_plan}, Int32, Int32, UInt32, UInt32, Int32),
        ptr,
        P.N,
        P.M,
        P.f1,
        P.f2,
        P.nfft_cutoff,
    )
    Core.setfield!(P, :init_done, true)
    finalizer(nfsft_finalize_plan, P)
end

function init(P::NFSFT)
    return nfsft_init(P)
end

# overwrite dot notation for plan struct in order to use C memory
function Base.setproperty!(P::NFSFT, v::Symbol, val)
    if !P.init_done
        nfsft_init(P)
    end
    if P.finalized
        error("NFSFT already finalized")
    end

    if v == :x
        if typeof(val) != Matrix{Float64}
            error("x has to be a Float64 matrix.")
        end
        if size(val)[1] != 2 || size(val)[2] != P.M
            error("x has to be a Float64 matrix of size 2xM.")
        end
        ptr = ccall(
            ("jnfsft_set_x", lib_path_nfsft),
            Ptr{Float64},
            (Ref{nfsft_plan}, Ref{Cdouble}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)
    elseif v == :f
        if typeof(val) != Vector{ComplexF64}
            error("f has to be a ComplexF64 vector.")
        end
        if size(val)[1] != P.M
            error("f has to be a ComplexF64 vector of length M.")
        end
        ptr = ccall(
            ("jnfsft_set_f", lib_path_nfsft),
            Ptr{ComplexF64},
            (Ref{nfsft_plan}, Ref{ComplexF64}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)
    elseif v == :fhat
        N_total = (2 * P.N + 2)^2
        if typeof(val) != Vector{ComplexF64}
            error("fhat has to be a ComplexF64 vector.")
        end
        if size(val)[1] != N_total
            error("fhat has to be a ComplexF64 vector of size N_total = (2N+2)^2.")
        end
        ptr = ccall(
            ("jnfsft_set_fhat", lib_path_nfsft),
            Ptr{ComplexF64},
            (Ref{nfsft_plan}, Ref{ComplexF64}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)
        Core.setfield!(P, :fhat_corrupted_by_fast_trafo, false)
    elseif v == :plan
        @warn "You can't modify the C pointer to the NFSFT plan."
    elseif v == :num_threads
        @warn "You can't currently modify the number of threads."
    elseif v == :init_done
        @warn "You can't modify this flag."
    elseif v == :N || v == :M || v == :f1 || v == :f2 || v == :nfft_cutoff
        @warn "You can't modify this parameter, please create an additional plan."
    else
        Core.setfield!(P, v, val)
    end
end

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(P::NFSFT, v::Symbol)
    if v == :x
        if !isdefined(P, :x)
            error("x is not set.")
        end
        ptr = Core.getfield(P, :x)
        return unsafe_wrap(Matrix{Float64}, ptr, (2, Int64(P.M)))
    elseif v == :num_threads
        return ccall(("nfft_get_num_threads", lib_path_nfsft), Int64, ())
    elseif v == :f
        if !isdefined(P, :f)
            error("f is not set.")
        end
        ptr = Core.getfield(P, :f)
        return unsafe_wrap(Vector{ComplexF64}, ptr, P.M)
    elseif v == :fhat
        if !isdefined(P, :fhat)
            error("fhat is not set.")
        end
        ptr = Core.getfield(P, :fhat)
        return unsafe_wrap(Vector{ComplexF64}, ptr, (2 * P.N + 2)^2)
    else
        return Core.getfield(P, v)
    end
end

@doc raw"""
    nfsft_index(P::NFSFT, k::Integer, n::Integer)

calculates the (1-based, Julia-indexable) index into `fhat` for the spherical Fourier coefficient of degree `k` and order `n`.

# See also
[`NFSFT`](@ref)
"""
function nfsft_index(P::NFSFT, k::Integer, n::Integer)
    return (2 * P.N + 2) * (P.N - n + 1) + (P.N + k + 1) + 1
end

@doc raw"""
    nfsft_trafo(P::NFSFT)

computes the NFSFT using the fast approximate transform for the provided nodes in `P.x` and coefficients in `P.fhat`.

# See also
[`NFSFT`](@ref), [`nfsft_trafo_direct`](@ref)
"""
function nfsft_trafo(P::NFSFT)
    if P.finalized
        error("NFSFT already finalized")
    end
    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr =
        ccall(("jnfsft_trafo", lib_path_nfsft), Ptr{ComplexF64}, (Ref{nfsft_plan},), P.plan)
    Core.setfield!(P, :f, ptr)
    # Without NFSFT_PRESERVE_F_HAT the C library overwrites fhat in place
    # with internal FPT/FFT scratch data as a side effect of trafo() (see
    # nfsft_trafo() in kernel/nfsft/nfsft.c) - remember that so
    # nfsft_trafo_direct() can refuse to trust it afterwards.
    if (P.f1 & NFSFT_PRESERVE_F_HAT) == 0
        Core.setfield!(P, :fhat_corrupted_by_fast_trafo, true)
    end
end

function trafo(P::NFSFT)
    return nfsft_trafo(P)
end

@doc raw"""
    nfsft_trafo_direct(P::NFSFT)

computes the NDSFT via direct (naive) transformation for the provided nodes in `P.x` and coefficients in `P.fhat`.

# See also
[`NFSFT`](@ref), [`nfsft_trafo`](@ref)
"""
function nfsft_trafo_direct(P::NFSFT)
    if P.finalized
        error("NFSFT already finalized")
    end
    if P.fhat_corrupted_by_fast_trafo
        error(
            "nfsft_trafo_direct() cannot be trusted on this plan: nfsft_trafo() " *
            "already overwrote fhat in place with internal scratch data (expected " *
            "behaviour of the underlying library without the NFSFT_PRESERVE_F_HAT " *
            "flag). Set fhat again, pass NFSFT_PRESERVE_F_HAT explicitly, or create " *
            "a new NFSFT plan.",
        )
    end
    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(
        ("jnfsft_trafo_direct", lib_path_nfsft),
        Ptr{ComplexF64},
        (Ref{nfsft_plan},),
        P.plan,
    )
    Core.setfield!(P, :f, ptr)
end

function trafo_direct(P::NFSFT)
    return nfsft_trafo_direct(P)
end

@doc raw"""
    nfsft_adjoint(P::NFSFT)

computes the adjoint NFSFT using the fast approximate transform for the provided nodes in `P.x` and coefficients in `P.f`.

# See also
[`NFSFT`](@ref), [`nfsft_adjoint_direct`](@ref)
"""
function nfsft_adjoint(P::NFSFT)
    if P.finalized
        error("NFSFT already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(
        ("jnfsft_adjoint", lib_path_nfsft),
        Ptr{ComplexF64},
        (Ref{nfsft_plan},),
        P.plan,
    )
    Core.setfield!(P, :fhat, ptr)
    # nfsft_adjoint() writes a fresh fhat (it does not read/reuse the
    # previous one), so any earlier trafo()-induced corruption no longer
    # matters.
    Core.setfield!(P, :fhat_corrupted_by_fast_trafo, false)
end

function adjoint(P::NFSFT)
    return nfsft_adjoint(P)
end

@doc raw"""
    nfsft_adjoint_direct(P::NFSFT)

computes the adjoint NDSFT via direct (naive) transformation for the provided nodes in `P.x` and coefficients in `P.f`.

# See also
[`NFSFT`](@ref), [`nfsft_adjoint`](@ref)
"""
function nfsft_adjoint_direct(P::NFSFT)
    if P.finalized
        error("NFSFT already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    if !isdefined(P, :x)
        error("x has not been set.")
    end
    ptr = ccall(
        ("jnfsft_adjoint_direct", lib_path_nfsft),
        Ptr{ComplexF64},
        (Ref{nfsft_plan},),
        P.plan,
    )
    Core.setfield!(P, :fhat, ptr)
    Core.setfield!(P, :fhat_corrupted_by_fast_trafo, false)
end

function adjoint_direct(P::NFSFT)
    return nfsft_adjoint_direct(P)
end
