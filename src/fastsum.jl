mutable struct fastsum_plan end

@doc raw"""
    FASTSUM

The fast summation algorithm evaluates the function

```math
f(\pmb{y}) \coloneqq \sum^{N}_{k = 1} \alpha_k \, \mathscr{K}(\pmb{y} - \pmb{x}_k) = \sum^{N}_{k = 1} \alpha_k \, K(\lVert \pmb{y} - \pmb{x}_k \rVert_2)
```

for given arbitrary source knots ``\pmb{x}_k \in \mathbb{R}^d, k = 1,2, \cdots, N`` and a given kernel function ``\mathscr{K}(\cdot) = K (\lVert \cdot \rVert_2), \; \pmb{x} \in \mathbb{R}^d``, 
which is an even, real univariate function which is infinitely differentiable at least in ``\mathbb{R} \setminus \{ 0 \}``. 
If ``K`` is infinitely differentiable at zero as well, then ``\mathscr{K}`` is defined on ``\mathbb{R}^d`` and is called 
nonsingular kernel function. The evaluation is done at ``M`` different points ``\pmb{y}_j \in \mathbb{R}^d, j = 1, \cdots, M``. 

# Fields
* `d` - dimension.
* `N` - number of source nodes.
* `M` - number of target nodes.
* `n` - expansion degree.
* `p` - degree of smoothness.
* `kernel` - name of kernel function ``K``.
* `c` - kernel parameters.
* `eps_I` - inner boundary.
* `eps_B` - outer boundary.
* `nn_x` - oversampled nn in x.
* `nn_y` - oversampled nn in y.
* `m_x` - NFFT-cutoff in x.
* `m_y` - NFFT-cutoff in y.
* `init_done` - bool for plan init.
* `finalized` - bool for finalizer.
* `flags` - flags.
* `x` - source nodes.
* `y` - target nodes.
* `alpha` - source coefficients.
* `f` - target evaluations.
* `plan` - plan (C pointer).

# Constructor
    FASTSUM( d::Integer, N::Integer, M::Integer, n::Integer, p::Integer, kernel::String, c::Vector{<:Real}, eps_I::Real, eps_B::Real, nn_x::Integer, nn_y::Integer, m_x::Integer, m_y::Integer, flags::UInt32 )

# Additional Constructor
    FASTSUM( d::Integer, N::Integer, M::Integer, n::Integer, p::Integer, kernel::String, c::Real, eps_I::Real, eps_B::Real, nn::Integer, m::Integer )

# See also
[`NFFT`](@ref)
"""
mutable struct FASTSUM
    d::Integer # dimension
    N::Integer # number of source nodes
    M::Integer # number of target nodes
    n::Integer # expansion degree
    p::Integer # degree of smoothness
    kernel::String # name of kernel
    c::Vector{Float64} # kernel parameters
    eps_I::Real # inner boundary
    eps_B::Real # outer boundary
    nn_x::Integer # oversampled nn in x
    nn_y::Integer # oversampled nn in y
    m_x::Integer # NFFT-cutoff in x
    m_y::Integer # NFFT-cutoff in y
    init_done::Bool # bool for plan init
    finalized::Bool # bool for finalizer
    flags::UInt32 # flags
    x::Ref{Float64} # source nodes
    y::Ref{Float64} # target nodes
    alpha::Ref{ComplexF64} # source coefficients
    f::Ref{ComplexF64} # target evaluations
    plan::Ref{fastsum_plan} # plan (C pointer)

    function FASTSUM(
        d::Integer,
        N::Integer,
        M::Integer,
        n::Integer,
        p::Integer,
        kernel::String,
        c::Vector{<:Real},
        eps_I::Real,
        eps_B::Real,
        nn_x::Integer,
        nn_y::Integer,
        m_x::Integer,
        m_y::Integer,
        flags::UInt32,
    )
        if N <= 0
            throw(DomainError(M, "argument must be a positive integer"))
        end

        if M <= 0
            throw(DomainError(M, "argument must be a positive integer"))
        end

        if n <= 0
            throw(DomainError(n, "argument must be a positive integer"))
        end

        if m_x <= 0 
            throw(DomainError(m_x, "argument must be a positive integer"))
        end

        if m_y <= 0 
            throw(DomainError(m_y, "argument must be a positive integer"))
        end

        new(
            d,
            N,
            M,
            n,
            p,
            kernel,
            Vector{Float64}(c),
            eps_I,
            eps_B,
            nn_x,
            nn_y,
            m_x,
            m_y,
            false,
            false,
            flags,
        )
    end
end #struct fastsumplan

function FASTSUM(
    d::Integer,
    N::Integer,
    M::Integer,
    kernel::String,
    c::Real,
    n::Integer = 256,
    p::Integer = 8,
    eps_I::Real = 256/8,
    eps_B::Real = 1/16,
    nn::Integer = 512,
    m::Integer = 8,
)
    cv = Vector{Float64}(undef, 1)
    cv[1] = Float64(c)

    FASTSUM(d, N, M, n, p, kernel, cv, eps_I, eps_B, nn, nn, m, m, UInt32(0))
end #constructor

@doc raw"""
    fastsum_init(P)

intialises a transform plan.

# Input
* `P` - a FASTSUM plan structure.

# See also
[`FASTSUM{D}`](@ref), [`fastsum_finalize_plan`](@ref)
"""
function fastsum_init(P::FASTSUM)

    ptr = ccall(("jfastsum_alloc", lib_path_fastsum), Ptr{fastsum_plan}, ())
    Core.setfield!(P, :plan, ptr)

    code = ccall(
        ("jfastsum_init", lib_path_fastsum),
        Int64,
        (
            Ref{fastsum_plan},
            Int32,
            Cstring,
            Ref{Float64},
            UInt32,
            Int32,
            Int32,
            Float64,
            Float64,
            Int32,
            Int32,
            Int32,
            Int32,
            Int32,
            Int32,
        ),
        ptr,
        Int32(P.d),
        P.kernel,
        P.c,
        P.flags,
        Int32(P.n),
        Int32(P.p),
        Float64(P.eps_I),
        Float64(P.eps_B),
        Int32(P.N),
        Int32(P.M),
        Int32(P.nn_x),
        Int32(P.nn_y),
        Int32(P.m_x),
        Int32(P.m_y),
    )

    if code == 1
        error("Unkown kernel.")
    end

    Core.setfield!(P, :init_done, true)
    finalizer(fastsum_finalize_plan, P)

end #fastsum_init

function init(P::FASTSUM)
    return fastsum_init(P)
end #fastsum_init

@doc raw"""
    fastsum_finalize_plan(P)

destroys a FASTSUM plan structure.

# Input
* `P` - a FASTSUM plan structure.

# See also
[`FASTSUM{D}`](@ref), [`fastsum_init`](@ref)
"""
function fastsum_finalize_plan(P::FASTSUM)

    if !P.init_done
        error("FASTSUM not initialized.")
    end

    if !P.finalized
        ccall(
            ("jfastsum_finalize", lib_path_fastsum),
            Nothing,
            (Ref{fastsum_plan},),
            P.plan,
        )
        Core.setfield!(P, :finalized, true)
    end

end #fastsum_finalize_plan

function finalize_plan(P::FASTSUM)
    return fastsum_finalize_plan(P)
end #fastsum_finalize_plan

function Base.setproperty!(P::FASTSUM, v::Symbol, val)

    if !P.init_done
        fastsum_init(P)
    end

    if P.finalized
        error("FASTSUM already finalized")
    end

    # edit source nodes
    if v == :x

        if P.d == 1
            if typeof(val) != Vector{Float64}
                error("x has to be a Float64 vector.")
            end
            if size(val)[1] != P.N
                error("x has to be a Float64 vector of length N.")
            end
        else # => D >1
            if typeof(val) != Array{Float64,2}
                error("x has to be a Float64 matrix.")
            end
            if size(val)[1] != P.N || size(val)[2] != P.d
                error("x has to be a Float64 matrix of size N.")
            end
        end

        ptr = ccall(
            ("jfastsum_set_x", lib_path_fastsum),
            Ptr{Float64},
            (Ref{fastsum_plan}, Ref{Cdouble}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)

        # edit target nodes
    elseif v == :y

        if P.d == 1
            if typeof(val) != Vector{Float64}
                error("y has to be a Float64 vector.")
            end
            if size(val)[1] != P.M
                error("y has to be a Float64 vector of length M.")
            end
        else # => D > 1
            if typeof(val) != Array{Float64,2}
                error("y has to be a Float64 matrix.")
            end
            if size(val)[1] != P.M || size(val)[2] != P.d
                error("y has to be a Float64 matrix of size M.")
            end
        end

        ptr = ccall(
            ("jfastsum_set_y", lib_path_fastsum),
            Ptr{Float64},
            (Ref{fastsum_plan}, Ref{Cdouble}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)

        # edit source coefficients
    elseif v == :alpha

        if typeof(val) != Vector{ComplexF64}
            error("alpha has to be a ComplexF64 vector.")
        end
        if size(val)[1] != P.N
            error("alpha has to be a ComplexF64 vector of length N.")
        end

        ptr = ccall(
            ("jfastsum_set_alpha", lib_path_fastsum),
            Ptr{ComplexF64},
            (Ref{fastsum_plan}, Ref{ComplexF64}),
            P.plan,
            val,
        )
        Core.setfield!(P, v, ptr)

    elseif v == :M
        @warn("You can't modify the number of target nodes.")
    elseif v == :N
        @warn("You can't modify the number of source nodes.")
    elseif v == :n
        @warn("You can't modify the expansion degree.")
    elseif v == :m
        @warn("You can't modify the cut-off parameter.")
    elseif v == :p
        @warn("You can't modify the degree of smoothness.")
    elseif v == :kernel
        @warn("You can't modify the kernel.")
    elseif v == :c
        @warn("You can't modify the kernel parameters.")
    elseif v == :eps_I
        @warn("You can't modify the inner boundary.")
    elseif v == :eps_B
        @warn("You can't modify the outer boundary.")
    elseif v == :plan
        @warn("You can't modify the pointer to the fastsum plan.")
    else
        Core.setfield!(P, v, val)
    end

end # Base.setproperty!

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(P::FASTSUM, v::Symbol)

    if v == :x

        if !isdefined(P, :x)
            error("x is not set.")
        end

        ptr = Core.getfield(P, :x)

        if P.d == 1
            return unsafe_wrap(Vector{Float64}, ptr, P.N) # get source nodes from C memory and convert to Julia type
        else
            return unsafe_wrap(Matrix{Float64}, ptr, (P.d, P.N))  # get source nodes from C memory and convert to Julia type
        end

    elseif v == :y

        if !isdefined(P, :y)
            error("y is not set.")
        end

        ptr = Core.getfield(P, :y)

        if P.d == 1
            return unsafe_wrap(Vector{Float64}, ptr, P.M)
        else
            return unsafe_wrap(Matrix{Float64}, ptr, (P.d, P.M))
        end

    elseif v == :alpha

        if !isdefined(P, :alpha)
            error("alpha is not set.")
        end

        ptr = Core.getfield(P, :alpha)
        return unsafe_wrap(Vector{ComplexF64}, ptr, P.N) # get coefficients from C memory and convert to Julia type

    elseif v == :f

        if !isdefined(P, :f)
            error("f is not set.")
        end

        ptr = Core.getfield(P, :f)
        return unsafe_wrap(Vector{ComplexF64}, ptr, P.M)  # get function values from C memory and convert to Julia type

    else
        return Core.getfield(P, v)
    end
end # Base.getproperty

@doc raw"""
    fastsum_trafo(P)

fast NFFT-based summation.

# Input
* `P` - a FASTSUM plan structure.

# See also
[`FASTSUM{D}`](@ref), [`fastsum_trafo_exact`](@ref)
"""
function fastsum_trafo(P::FASTSUM)

    if P.finalized
        error("FASTSUM already finalized.")
    end

    if !isdefined(P, :x)
        error("x has not been set.")
    end

    if !isdefined(P, :y)
        error("y has not been set.")
    end

    if !isdefined(P, :alpha)
        error("alpha has not been set.")
    end

    ptr = ccall(
        ("jfastsum_trafo", lib_path_fastsum),
        Ptr{ComplexF64},
        (Ref{fastsum_plan},),
        P.plan,
    )
    Core.setfield!(P, :f, ptr)

end #trafo

function trafo(P::FASTSUM)
    return fastsum_trafo(P)
end #trafo

@doc raw"""
    fastsum_trafo_exact(P)

direct computation of sums.

# Input
* `P` - a FASTSUM plan structure.

# See also
[`FASTSUM{D}`](@ref), [`fastsum_trafo`](@ref)
"""
function fastsum_trafo_exact(P::FASTSUM)

    if P.finalized
        error("FASTSUM already finalized.")
    end

    if !isdefined(P, :x)
        error("x has not been set.")
    end

    if !isdefined(P, :y)
        error("y has not been set.")
    end

    if !isdefined(P, :alpha)
        error("alpha has not been set.")
    end

    ptr = ccall(
        ("jfastsum_exact", lib_path_fastsum),
        Ptr{ComplexF64},
        (Ref{fastsum_plan},),
        P.plan,
    )
    Core.setfield!(P, :f, ptr)

end #trafo

function trafo_exact(P::FASTSUM)
    return fastsum_trafo_exact(P)
end #trafo