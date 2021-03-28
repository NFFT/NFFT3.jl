# dummy struct for C
mutable struct nfct_plan end

@doc raw"""
    NFCT{D}

A NFCT (Nonequispaced fast cosine transform) plan, where D is the dimension. 

The NFCT realizes a direct and fast computation of the discrete nonequispaced cosine transform. The aim is to compute

```math
f^c (x) = \sum_{k \in I_N} \hat{f}^{c}_{k} \cos (k \dot x), \quad x \in \mathbb{R}^D
```

at given (nonequidistant) knots ``x_k \in [0, \pi ]^D, \; k = 0, \cdots, M-1``, coefficients ``\hat{f}^{c}_{k} \in \mathbb{R}``, `k \in I_N \colon = \{ k \in \mathbb{Z}^{D} \colon 0 \leq k_i \leq N_i, \, \forall i = 1, \cdots, D\}`` for some multibandlimit vector ``N \in \mathbb{N}^{D}``. The transposed (adjoined) problem reads as

```math
\hat{h}_k \colon = \sum_{j \in I_M} f_j \cos (k \dot x_j), \quad k \in I^{D}_N
```

for given knots ``x_k \in [0, \pi ]^D, \; k = 0, \cdots, M-1``, and coefficients fj ∈ C, j ∈ I
l
M.

## Fields
* `N` - the multibandlimit of the trigonometric polynomial f.
* `M` - the number of nodes.
* `n` - the oversampling per dimension.
* `m` - the window size. Larger m means more accuracy but also more computational costs. 
* `f1` - the NFCT flags.
* `f2` - the FFTW flags.
* `init_done` - indicates if the plan is initialized.
* `finalized` - indicates if the plan is finalized.
* `x` - the nodes.
* `f` - the function values.
* `fhat` - the Fourier coefficients.
* `plan`

# Constructor
    NFCT{D}(N::NTuple{D,Int32},M::Int32,n::NTuple{D,Int32},m::Int32,f1::UInt32,f2::UInt32) where {D}

# See also
[`NFFT`](@ref)

[^PlonkaPottsSteidelTasche2018]:
    > Plonka, Gerlind and Potts, Daniel and Steidl, Gabriele and Tasche, Manfred 
    > Numerical Fourier Analysis
    > 2018

[^KeinerKunisPotts]:
    > J. Keiner, S. Kunis, and D. Potts
    > NFFT 3.0, C subroutine library
    > url: http://www.tu-chemnitz.de/~potts/nfft.
"""
# NFCT plan struct
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
        # create plan object
        new(N, M, n, m, f1, f2, false, false)
    end
end

# additional constructor for easy use [NFCT((N,N),M) instead of NFCT{2}((N,N),M)]
function NFCT(N::NTuple{D,Integer}, M::Integer) where {D}
    if any(x -> x <= 0, N)
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
function finalize_plan(P::NFCT{D}) where {D}
    if !P.init_done
        error("NFCT not initialized.")
    end

    if !P.finalized
        Core.setfield!(P, :finalized, true)
        ccall(("jnfct_finalize", lib_path_nfct), Nothing, (Ref{nfct_plan},), P.plan)
    end
end

# allocate plan memory and init with D,N,M,n,m,f1,f2
function nfct_init(p::NFCT{D}) where {D}
    # convert N and n to vectors for passing them over to C
    Nv = collect(p.N)
    n = collect(p.n)

    # call init for memory allocation
    ptr = ccall(("jnfct_alloc", lib_path_nfct), Ptr{nfct_plan}, ())

    # set pointer
    Core.setfield!(p, :plan, ptr)

    # initialize values
    ccall(
        ("jnfct_init", lib_path_nfct),
        Nothing,
        (Ref{nfct_plan}, Int32, Ref{Int32}, Int32, Ref{Int32}, Int32, UInt32, UInt32),
        ptr,
        D,
        Nv,
        p.M,
        n,
        p.m,
        p.f1,
        p.f2,
    )
    Core.setfield!(p, :init_done, true)
    finalizer(finalize_plan, p)
end

# overwrite dot notation for plan struct in order to use C memory
function Base.setproperty!(p::NFCT{D}, v::Symbol, val) where {D}
    # init plan if not done [usually with setting nodes]
    if !p.init_done
        nfct_init(p)
    end

    # prevent bad stuff from happening
    if p.finalized
        error("NFCT already finalized")
    end

    # setting nodes, verification of correct size dxM
    if v == :x
        if D == 1
            if typeof(val) != Vector{Float64}
                error("x has to be a Float64 vector.")
            end
            if size(val)[1] != p.M
                error("x has to be a Float64 vector of length M.")
            end
        else
            if typeof(val) != Array{Float64,2}
                error("x has to be a Float64 matrix.")
            end
            if size(val)[1] != D || size(val)[2] != p.M
                error("x has to be a Float64 matrix of size dxM.")
            end
        end
        ptr = ccall(
            ("jnfct_set_x", lib_path_nfct),
            Ptr{Float64},
            (Ref{nfct_plan}, Ref{Cdouble}),
            p.plan,
            val,
        )
        Core.setfield!(p, v, ptr)

        # setting values
    elseif v == :f
        if typeof(val) != Array{Float64,1}
            error("f has to be a Float64 vector.")
        end
        if size(val)[1] != p.M
            error("f has to be a Float64 vector of size M.")
        end
        ptr = ccall(
            ("jnfct_set_f", lib_path_nfct),
            Ptr{Float64},
            (Ref{nfct_plan}, Ref{Float64}),
            p.plan,
            val,
        )
        Core.setfield!(p, v, ptr)

        # setting Fourier coefficients
    elseif v == :fhat
        if typeof(val) != Array{Float64,1}
            error("fhat has to be a Float64 vector.")
        end
        l = prod(p.N)
        if size(val)[1] != l
            error("fhat has to be a Float64 vector of size prod(N).")
        end
        ptr = ccall(
            ("jnfct_set_fhat", lib_path_nfct),
            Ptr{Float64},
            (Ref{nfct_plan}, Ref{Float64}),
            p.plan,
            val,
        )
        Core.setfield!(p, v, ptr)

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
        Core.setfield!(p, v, val)
    end
end

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(p::NFCT{D}, v::Symbol) where {D}
    if v == :x
        if !isdefined(p, :x)
            error("x is not set.")
        end
        ptr = Core.getfield(p, :x)
        if D == 1
            return unsafe_wrap(Vector{Float64}, ptr, p.M)             # get nodes from C memory and convert to Julia type
        else
            return unsafe_wrap(Matrix{Float64}, ptr, (D, Int64(p.M)))  # get nodes from C memory and convert to Julia type
        end
    elseif v == :num_threads
        return ccall(("nfft_get_num_threads", lib_path_nfct), Int64, ())
    elseif v == :f
        if !isdefined(p, :f)
            error("f is not set.")
        end
        ptr = Core.getfield(p, :f)
        return unsafe_wrap(Vector{Float64}, ptr, p.M)  # get function values from C memory and convert to Julia type
    elseif v == :fhat
        if !isdefined(p, :fhat)
            error("fhat is not set.")
        end
        ptr = Core.getfield(p, :fhat)
        return unsafe_wrap(Vector{Float64}, ptr, prod(p.N)) # get Fourier coefficients from C memory and convert to Julia type
    else
        return Core.getfield(p, v)
    end
end

# nfct trafo direct [call with NFCT.trafo_direct outside module]
function trafo_direct(P::NFCT{D}) where {D}
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

# adjoint trafo direct [call with NFCT.adjoint_direct outside module]
function adjoint_direct(P::NFCT{D}) where {D}
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

# nfct trafo [call with NFCT.trafo outside module]
function trafo(P::NFCT{D}) where {D}
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

# adjoint trafo [call with NFCT.adjoint outside module]
function adjoint(P::NFCT{D}) where {D}
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
