# dummy struct for C
mutable struct fsft_plan end

@doc raw"""
    FSFT

A FSFT (fast spherical Fourier transform on an equispaced grid) plan.

FSFT is the `NFSFT_EQUISPACED` variant of [`NFSFT`](@ref): the nodes are not
user-supplied, they are the implicit equispaced grid of size `(2N+2)*(N+2)`
that `nfsft_init_guru()` always uses for `NFSFT_EQUISPACED` plans (the `M`
argument is ignored for exactly this reason - `FSFT` computes it itself).

Note: `nfsft_trafo_direct()`/`nfsft_adjoint_direct()` (shared with [`NFSFT`](@ref)),
including when called internally by `trafo()`/`adjoint()` below the
`NFSFT_BREAK_EVEN` bandwidth, unconditionally read the plan's node array for
node positions - but since FSFT's grid is implicit, there is no public `x`
setter. `fsft_init()` therefore sets the equispaced grid itself via
`jnfsft_set_x()`, mirroring the formula in `nfsft_init_guru()`
(`kernel/nfsft/nfsft.c`): without this, `trafo_direct()`/`adjoint_direct()`
read whatever was left in memory for the node array, which is not reliably
the equispaced grid in the bundled library build.

# Fields
* `N` - the bandwidth (maximum degree ``k``).
* `M` - the number of equispaced grid nodes, `(2N+2)*(N+2)` (computed internally; any `M` passed to the constructor is ignored).
* `f1` - the NFSFT flags (`NFSFT_EQUISPACED` is always included).
* `f2` - the NFFT flags used internally.
* `nfft_cutoff` - the NFFT window cutoff used internally.
* `f` - the values ``f(\pmb{x}_j)`` for the FSFT or the coefficients ``f_j`` for the adjoint FSFT.
* `fhat` - the spherical Fourier coefficients for the FSFT, or the values for the adjoint FSFT.
* `plan` - plan (C pointer).

# Constructor
    FSFT(N::Integer, M::Integer=0, f1::UInt32=nfsft_default|NFSFT_EQUISPACED, f2::UInt32=nfsft_nfft_default, nfft_cutoff::Integer=nfsft_default_nfft_cut_off)

# See also
[`NFSFT`](@ref) for the nonequispaced variant.
"""
mutable struct FSFT
    N::Int32
    M::Int32
    f1::UInt32
    f2::UInt32
    nfft_cutoff::Int32
    init_done::Bool
    finalized::Bool
    fhat_corrupted_by_fast_trafo::Bool
    f::Ref{ComplexF64}
    fhat::Ref{ComplexF64}
    plan::Ref{fsft_plan}
    function FSFT(
        N::Integer,
        M::Integer = 0,
        f1::UInt32 = nfsft_default | NFSFT_EQUISPACED,
        f2::UInt32 = nfsft_nfft_default,
        nfft_cutoff::Integer = nfsft_default_nfft_cut_off,
    )
        if N <= 0
            throw(DomainError(N, "argument must be a positive integer"))
        end
        if (f1 & NFSFT_EQUISPACED) == 0
            error("Must use NFSFT_EQUISPACED flag. Use NFSFT instead")
        end
        # fsft_init() below unconditionally copies into the native f/fhat
        # buffers (and sets the equispaced grid via jnfsft_set_x), so the C
        # library must always allocate them regardless of which flags the
        # caller passed - see the identical reasoning in NFSFT.jl's
        # constructor. Without this, e.g.
        # FSFT(N, 0, NFSFT_EQUISPACED | NFSFT_PRESERVE_F_HAT) silently drops
        # the allocation flags and crashes on the first jnfsft_set_x call.
        f1 |= NFSFT_MALLOC_X | NFSFT_MALLOC_F | NFSFT_MALLOC_F_HAT
        # The underlying C library ignores the M passed to its init function
        # for equispaced plans and always allocates/uses the fixed grid size
        # (2N+2)(N+2) instead - mirror that here so M actually matches the
        # real buffer size.
        M_grid = (2 * N + 2) * (N + 2)
        new(Int32(N), Int32(M_grid), f1, f2, Int32(nfft_cutoff), false, false, false)
    end
end

function fsft_finalize_plan(P::FSFT)
    if !P.init_done
        error("FSFT not initialized.")
    end
    if !P.finalized
        Core.setfield!(P, :finalized, true)
        ccall(("jnfsft_finalize", lib_path_nfsft), Nothing, (Ref{fsft_plan},), P.plan)
    end
end

function finalize_plan(P::FSFT)
    return fsft_finalize_plan(P)
end

function fsft_init(P::FSFT)
    ptr = ccall(("jnfsft_alloc", lib_path_nfsft), Ptr{fsft_plan}, ())
    Core.setfield!(P, :plan, ptr)

    ccall(
        ("jnfsft_init", lib_path_nfsft),
        Nothing,
        (Ref{fsft_plan}, Int32, Int32, UInt32, UInt32, Int32),
        ptr,
        P.N,
        P.M,
        P.f1,
        P.f2,
        P.nfft_cutoff,
    )
    Core.setfield!(P, :init_done, true)
    finalizer(fsft_finalize_plan, P)

    # Set the equispaced grid ourselves - see the struct docstring above and
    # the "Set equispaced nodes" comment in nfsft_init_guru() this mirrors:
    #   x[2*(i*(N+2)+j)]   = (i-N-1)/(2N+2)   for i in 0..2N+1
    #   x[2*(i*(N+2)+j)+1] =      j/(2N+2)    for j in 0..N+1
    # jnfsft_set_x expects its input in "raw angle" units (it divides by
    # 2*pi itself, with a wraparound for the first component) rather than
    # the plan's own normalized units directly, so the grid values below are
    # first scaled back up by 2*pi to compensate.
    N = P.N
    grid = Matrix{Float64}(undef, 2, P.M)
    for i = 0:(2*N+1), j = 0:(N+1)
        m = i * (N + 2) + j + 1  # 1-based column index into grid
        grid[1, m] = 2.0 * pi * (i - N - 1) / (2.0 * N + 2)
        grid[2, m] = 2.0 * pi * j / (2.0 * N + 2)
    end
    ccall(
        ("jnfsft_set_x", lib_path_nfsft),
        Ptr{Float64},
        (Ref{fsft_plan}, Ref{Cdouble}),
        P.plan,
        grid,
    )
end

function init(P::FSFT)
    return fsft_init(P)
end

function Base.setproperty!(P::FSFT, v::Symbol, val)
    if !P.init_done
        fsft_init(P)
    end
    if P.finalized
        error("FSFT already finalized")
    end

    if v == :f
        if typeof(val) != Vector{ComplexF64}
            error("f has to be a ComplexF64 vector.")
        end
        if size(val)[1] != P.M
            error("f has to be a ComplexF64 vector of length M.")
        end
        ptr = ccall(
            ("jnfsft_set_f", lib_path_nfsft),
            Ptr{ComplexF64},
            (Ref{fsft_plan}, Ref{ComplexF64}),
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
            (Ref{fsft_plan}, Ref{ComplexF64}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)
        Core.setfield!(P, :fhat_corrupted_by_fast_trafo, false)
    elseif v == :plan
        @warn "You can't modify the C pointer to the FSFT plan."
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

function Base.getproperty(P::FSFT, v::Symbol)
    if v == :num_threads
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
    fsft_index(P::FSFT, k::Integer, n::Integer)

calculates the (1-based, Julia-indexable) index into `fhat` for the spherical Fourier coefficient of degree `k` and order `n`.

# See also
[`FSFT`](@ref)
"""
function fsft_index(P::FSFT, k::Integer, n::Integer)
    return (2 * P.N + 2) * (P.N - n + 1) + (P.N + k + 1) + 1
end

@doc raw"""
    fsft_trafo(P::FSFT)

computes the FSFT using the fast approximate transform for the coefficients in `P.fhat`.

# See also
[`FSFT`](@ref), [`fsft_trafo_direct`](@ref)
"""
function fsft_trafo(P::FSFT)
    if P.finalized
        error("FSFT already finalized")
    end
    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end
    ptr =
        ccall(("jnfsft_trafo", lib_path_nfsft), Ptr{ComplexF64}, (Ref{fsft_plan},), P.plan)
    Core.setfield!(P, :f, ptr)
    # Same underlying nfsft_trafo() as NFSFT: without NFSFT_PRESERVE_F_HAT it
    # overwrites fhat in place with internal FPT/FFT scratch data - but only
    # for P.N >= NFSFT_BREAK_EVEN or with NFSFT_NORMALIZED set, see the
    # identical reasoning in NFSFT.jl's nfsft_trafo().
    if (P.f1 & NFSFT_PRESERVE_F_HAT) == 0 &&
       (P.N >= NFSFT_BREAK_EVEN || (P.f1 & NFSFT_NORMALIZED) != 0)
        Core.setfield!(P, :fhat_corrupted_by_fast_trafo, true)
    end
end

function trafo(P::FSFT)
    return fsft_trafo(P)
end

@doc raw"""
    fsft_trafo_direct(P::FSFT)

computes the NDSFT via direct (naive) transformation for the coefficients in `P.fhat`.

# See also
[`FSFT`](@ref), [`fsft_trafo`](@ref)
"""
function fsft_trafo_direct(P::FSFT)
    if P.finalized
        error("FSFT already finalized")
    end
    if P.fhat_corrupted_by_fast_trafo
        error(
            "fsft_trafo_direct() cannot be trusted on this plan: fsft_trafo() " *
            "already overwrote fhat in place with internal scratch data (expected " *
            "behaviour of the underlying library without the NFSFT_PRESERVE_F_HAT " *
            "flag). Set fhat again, pass NFSFT_PRESERVE_F_HAT explicitly, or create " *
            "a new FSFT plan.",
        )
    end
    if !isdefined(P, :fhat)
        error("fhat has not been set.")
    end
    ptr = ccall(
        ("jnfsft_trafo_direct", lib_path_nfsft),
        Ptr{ComplexF64},
        (Ref{fsft_plan},),
        P.plan,
    )
    Core.setfield!(P, :f, ptr)
end

function trafo_direct(P::FSFT)
    return fsft_trafo_direct(P)
end

@doc raw"""
    fsft_adjoint(P::FSFT)

computes the adjoint FSFT using the fast approximate transform for the coefficients in `P.f`.

# See also
[`FSFT`](@ref), [`fsft_adjoint_direct`](@ref)
"""
function fsft_adjoint(P::FSFT)
    if P.finalized
        error("FSFT already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    ptr = ccall(
        ("jnfsft_adjoint", lib_path_nfsft),
        Ptr{ComplexF64},
        (Ref{fsft_plan},),
        P.plan,
    )
    Core.setfield!(P, :fhat, ptr)
    Core.setfield!(P, :fhat_corrupted_by_fast_trafo, false)
end

function adjoint(P::FSFT)
    return fsft_adjoint(P)
end

@doc raw"""
    fsft_adjoint_direct(P::FSFT)

computes the adjoint NDSFT via direct (naive) transformation for the coefficients in `P.f`.

# See also
[`FSFT`](@ref), [`fsft_adjoint`](@ref)
"""
function fsft_adjoint_direct(P::FSFT)
    if P.finalized
        error("FSFT already finalized")
    end
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    ptr = ccall(
        ("jnfsft_adjoint_direct", lib_path_nfsft),
        Ptr{ComplexF64},
        (Ref{fsft_plan},),
        P.plan,
    )
    Core.setfield!(P, :fhat, ptr)
    Core.setfield!(P, :fhat_corrupted_by_fast_trafo, false)
end

function adjoint_direct(P::FSFT)
    return fsft_adjoint_direct(P)
end
