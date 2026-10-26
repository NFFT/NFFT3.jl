using NFFT3

# NFMT plan struct
@doc raw"""
    NFMT{D}

A NFMT (nonequispaced fast mixed transform) plan, where D is the dimension.

The NFMT realizes a direct and fast computation of the discrete nonequispaced mixed transform. The aim is to compute

$$f^{\pmb{d}}(\pmb{x}_j) \coloneqq \sum_{\pmb{k} \in I_{\pmb{N},\pmb{d}}^d} \hat{f}_{\pmb{k}}^{\pmb{d}} \, \phi_{\pmb{k}}^{\pmb{d}}(\pmb{x}_j)$$

with

$$\phi_{\pmb{k}}^{\pmb{d}}(\pmb{x})=\prod_{j=1}^d\begin{cases}1,&k_j=0\\\exp(2\pi\mathrm{i}k_jx_j),&d_j=\exp,\;k_j\neq0\\
\sqrt{2}\cos(\pi k_jx_j),&d_j=\cos,\;k_j\neq0\\
\sqrt{2}\cos(k_j\arccos(2x_j-1)),&d_j=\mathrm{alg},\;k_j\neq0\end{cases}$$

for $\pmb{d}\in\{\exp,\cos,\mathrm{alg}\}^d$ at given arbitrary knots ``\pmb{x}_j \in [0,1]^D, j = 1, \cdots, M``, for coefficients ``\hat{f}^{c}_{\pmb{k}}``,

$$\pmb{k}\inI_{\pmb{N},\pmb{d}}^d \coloneqq \overset{d}{\underset{j=1}{\vphantom{\mathop{\raisebox{-.5ex}{\hbox{\huge{$\times$}}}}}⨉}}\begin{cases}\Big\{-\frac{N_j}{2},-\frac{N_j}{2}+1,\ldots,\frac{N_j}{2}\Big\},&d_j=\exp\\\Big\{0,1,\ldots,\frac{N_j}{2}\Big\},&d_j\neq\exp\end{cases},$$

and a multibandlimit vector ``\pmb{N} \in (2\mathbb{N})^{D}``. The transposed problem reads as

$$\hat{h}^{\pmb{d}}_{\pmb{k}} = \sum_{j=1}^M f^{\pmb{d}}_j \, \phi_{\pmb{k}}^{\pmb{d}}(\pmb{x}_j)$$

for the frequencies ``\pmb{k} \in I_{\pmb{N},\pmb{d}}^D`` with given coefficients ``f^{\pmb{d}}_j, j = 1,2,\ldots,M``.

Internally this reduces the mixed transform to a plain [`NFFT`](@ref) or
[`NFCT`](@ref): node coordinates are remapped per dimension (`cos`/`alg`
nodes via a scale/arccos substitution), and the one-sided mixed-basis
coefficients are expanded into/folded back from the backend's own
coefficient array. `fhat`/`f` may be assigned as either a `ComplexF64` or a
real `Float64` array; the backend is chosen lazily from whichever is
assigned first: a real array picks the cheaper NFCT backend, but only when
every dimension uses a `cos`/`alg` basis (the only case where the result is
guaranteed real) - a complex array, or any `exp` dimension, always picks the
NFFT backend (upcasting real input to `ComplexF64` if needed). Once a plan
has picked the real NFCT backend, assigning it complex data raises an
error - create a new plan instead. The NFCT backend additionally requires
`N/2` to be even in every dimension (NFCT's own even-bandwidth requirement);
otherwise use an `N` that is a multiple of 4, or pass complex data to fall
back to the NFFT backend.

# Fields
* `basis_vect` - tuple with the bases (`exp`, `cos`, `alg`) for each dimension.
* `N` - the multibandlimit ``(N_1, N_2, \ldots, N_D)``.
* `M` - the number of nodes.
* `n` - the oversampling per dimension, or `nothing` to let the backend pick its own default.
* `m` - the window size.
* `f1` - the NFFT/NFCT flags (without the MALLOC_*/FFTW_INIT bits, which are added automatically).
* `f2` - the FFTW flags.
* `backend` - the wrapped [`NFFT`](@ref) or [`NFCT`](@ref) plan, or `nothing` before it is chosen.

# Constructor
    NFMT{D}( basis_vect::NTuple{D,String}, N::NTuple{D,Int32}, M::Int32, n::Union{Nothing,NTuple{D,Int32}}, m::Int32, f1::UInt32, f2::UInt32 ) where {D}

# Additional Constructor
    NFMT( basis_vect::NTuple{D,String}, N::NTuple{D,Integer}, M::Integer, n::NTuple{D,Integer}, m::Integer, f1::UInt32, f2::UInt32) where {D}
    NFMT( basis_vect::NTuple{D,String}, N::NTuple{D,Integer}, M::Integer) where {D}

# See also
[`NFMT`](@ref)
"""
mutable struct NFMT{D}
    basis_vect::NTuple{D,String}
    N::NTuple{D,Int32}
    M::Int32
    n::Union{Nothing,NTuple{D,Int32}}
    m::Int32
    f1::UInt32
    f2::UInt32
    backend::Union{Nothing,NFFT{D},NFCT{D}}
    x_raw::Union{Nothing,Array{Float64}}
    function NFMT{D}(
        basis_vect::NTuple{D,String},
        N::NTuple{D,Int32},
        M::Int32,
        n::Union{Nothing,NTuple{D,Int32}},
        m::Int32,
        f1::UInt32,
        f2::UInt32,
    ) where {D}
        if length(basis_vect) != D
            error("basis_vect must have one entry per dimension")
        end
        for s in basis_vect
            if !haskey(BASES, s)
                error("Unknown basis '$s'. Allowed: $(collect(keys(BASES)))")
            end
        end
        if M <= 0
            throw(DomainError(M, "argument must be a positive integer"))
        end
        new(basis_vect, N, M, n, m, f1, f2, nothing, nothing)
    end
end

function NFMT(
    basis_vect::NTuple{D,String},
    N::NTuple{D,Integer},
    M::Integer,
    n::NTuple{D,Integer},
    m::Integer,
    f1::UInt32 = (D > 1 ? NFFT3.f1_default : NFFT3.f1_default_1d),
    f2::UInt32 = NFFT3.f2_default,
) where {D}
    NFMT{D}(
        basis_vect,
        NTuple{D,Int32}(N),
        Int32(M),
        NTuple{D,Int32}(n),
        Int32(m),
        f1,
        f2,
    )
end

# additional constructor for easy use [NFMT(basis_vect,N,M) instead of NFMT{2}(basis_vect,N,M)]
function NFMT(basis_vect::NTuple{D,String}, N::NTuple{D,Integer}, M::Integer) where {D}
    f1 = D > 1 ? NFFT3.f1_default : NFFT3.f1_default_1d
    NFMT{D}(basis_vect, NTuple{D,Int32}(N), Int32(M), nothing, Int32(8), f1, NFFT3.f2_default)
end

@doc raw"""
    nfmt_finalize_plan(P::NFMT{D})

destroys the wrapped NFFT/NFCT plan structure.

# See also
[`NFMT{D}`](@ref)
"""
function nfmt_finalize_plan(P::NFMT{D}) where {D}
    if P.backend !== nothing
        if P.backend isa NFFT{D}
            nfft_finalize_plan(P.backend)
        else
            nfct_finalize_plan(P.backend)
        end
    end
end

function finalize_plan(P::NFMT{D}) where {D}
    return nfmt_finalize_plan(P)
end

# True if no dimension uses the exp basis, i.e. the transform is real-valued
# and can be realised with a plain NFCT instead of an NFFT.
_nfmt_all_cos_alg(basis_vect) = all(BASES[s] > 0 for s in basis_vect)

function _nfmt_real_backend_bandwidths(N::NTuple{D,Int32}) where {D}
    half = Int32.(N .÷ 2)
    if any(h -> h % 2 != 0, half)
        error(
            "Real-valued NFMT (NFCT backend) requires N/2 to be even in every " *
            "dimension, since NFCT itself needs an even bandwidth; got N=$(N) -> " *
            "N/2=$(half). Use an N that is a multiple of 4 in every dimension, or " *
            "pass ComplexF64 fhat/f to use the general NFFT-based path instead.",
        )
    end
    return Tuple(half)
end

function _nfmt_ensure_backend!(P::NFMT{D}, want_complex::Bool) where {D}
    if P.backend isa NFCT{D} && want_complex
        error(
            "This NFMT plan already uses the real-valued NFCT backend (its first " *
            "fhat/f assignment was real); it cannot accept complex data afterwards. " *
            "Create a new NFMT plan to use complex coefficients on it.",
        )
    end
    if P.backend !== nothing
        return
    end

    f1_bits = P.f1 | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT
    backend = if want_complex || !_nfmt_all_cos_alg(P.basis_vect)
        if P.n === nothing
            # P.n is only nothing via the 3-arg NFMT convenience constructor,
            # which already set P.f1/P.m to the same defaults NFFT's own
            # 2-arg constructor would pick, so this is equivalent to passing
            # them through explicitly.
            NFFT(Tuple(P.N), Int(P.M))
        else
            NFFT{D}(P.N, P.M, P.n, P.m, f1_bits, P.f2)
        end
    else
        half_N = _nfmt_real_backend_bandwidths(P.N)
        if P.n === nothing
            NFCT(half_N, Int(P.M))
        else
            NFCT{D}(half_N, P.M, P.n, P.m, f1_bits, P.f2)
        end
    end
    Core.setfield!(P, :backend, backend)

    if P.x_raw !== nothing
        backend.x = _nfmt_transform_x(P, P.x_raw)
    end
end

# Maps nodes into the native domain [-0.5, 0.5) (NFFT) resp. [0, 0.5] (NFCT)
# - both backends share this image exactly for cos/alg dimensions, so this
# transform is the same regardless of which backend ends up being used.
function _nfmt_node_to_backend_domain(col, basis)
    b = BASES[basis]
    if b == 0  # exp: already on the NFFT's native domain [-0.5, 0.5)
        return col
    elseif b == 1  # cos on [0, 1] -> native domain via a half-scale
        return col ./ 2.0
    elseif b == 2  # alg (Chebyshev-type) on [0, 1] via arccos substitution
        return acos.(2.0 .* col .- 1.0) ./ (2.0 * pi)
    else
        error("Unknown basis id for $basis")
    end
end

function _nfmt_transform_x(P::NFMT{D}, xh) where {D}
    xh = copy(xh)
    if D == 1
        return _nfmt_node_to_backend_domain(xh, P.basis_vect[1])
    else
        out = similar(xh)
        for i = 1:D
            out[i, :] = _nfmt_node_to_backend_domain(xh[i, :], P.basis_vect[i])
        end
        return out
    end
end

# number of dimensions with a cos/alg basis
_nfmt_a(basis_vect) = sum(getindex.([BASES], basis_vect) .> 0)

_nfmt_reduced_size(N::NTuple{D,Int32}, basis_vect) where {D} =
    prod(N) ÷ (2^_nfmt_a(basis_vect))

function _nfmt_strides(N::NTuple{D,Int32}, basis_vect; reduced::Bool) where {D}
    pv = zeros(Int64, D)
    pv[D] = 1
    for i = D-1:-1:1
        pv[i] = pv[i+1] * N[i+1]
        if reduced && BASES[basis_vect[i+1]] > 0
            pv[i] ÷= 2
        end
    end
    return pv
end

# Decode linear index i (1-based) into the FULL (NFFT-sized) fhat array,
# returning (idx_red, e, valid): the matching 1-based index into the reduced
# array, the exponent for the sqrt(2)^e normalisation factor, and whether i
# corresponds to a Nyquist frequency that has no reduced counterpart.
function _nfmt_decode_full_index(i::Int, N::NTuple{D,Int32}, basis_vect, pv_exp, pv_red) where {D}
    idx = 1
    e = 0
    for j = 1:D
        if BASES[basis_vect[j]] > 0
            k = abs(((i - 1) ÷ pv_exp[j]) % N[j] - N[j] ÷ 2)
            if k == N[j] ÷ 2
                return 1, 0, false
            end
            if k != 0
                e -= 1
            end
            idx += k * pv_red[j]
        else
            idx += (((i - 1) ÷ pv_exp[j]) % N[j]) * pv_red[j]
        end
    end
    return idx, e, true
end

function _nfmt_expand_reduced_to_full(fhat_reduced, N::NTuple{D,Int32}, basis_vect) where {D}
    p = prod(N)
    pv_exp = _nfmt_strides(N, basis_vect; reduced = false)
    pv_red = _nfmt_strides(N, basis_vect; reduced = true)
    full = zeros(ComplexF64, p)
    for i = 1:p
        idx, e, valid = _nfmt_decode_full_index(i, N, basis_vect, pv_exp, pv_red)
        if valid
            full[i] = sqrt(2.0)^e * fhat_reduced[idx]
        end
    end
    return full
end

function _nfmt_reduce_full_to_reduced(fhat_full, N::NTuple{D,Int32}, basis_vect) where {D}
    p = prod(N)
    pv_exp = _nfmt_strides(N, basis_vect; reduced = false)
    pv_red = _nfmt_strides(N, basis_vect; reduced = true)
    reduced = zeros(ComplexF64, _nfmt_reduced_size(N, basis_vect))
    for i = 1:p
        idx, e, valid = _nfmt_decode_full_index(i, N, basis_vect, pv_exp, pv_red)
        if valid
            reduced[idx] += sqrt(2.0)^e * fhat_full[i]
        end
    end
    return reduced
end

# sqrt(2)^e per reduced index for the NFCT backend (every dimension cos/alg):
# NFCT's own N_j is exactly N_j/2 (see _nfmt_real_backend_bandwidths) and its
# native linear index already equals this module's reduced index one for
# one (same strides), with no folding or Nyquist exclusion needed.
function _nfmt_cos_alg_scale(N::NTuple{D,Int32}, basis_vect) where {D}
    half_N = Int.(N) .÷ 2
    pv_red = _nfmt_strides(N, basis_vect; reduced = true)
    nf = _nfmt_reduced_size(N, basis_vect)
    scale = zeros(Float64, nf)
    for i = 1:nf
        e = 0
        for j = 1:D
            k = ((i - 1) ÷ pv_red[j]) % half_N[j]
            if k != 0
                e += 1
            end
        end
        scale[i] = sqrt(2.0)^e
    end
    return scale
end

# overwrite dot notation for struct to set values and do the transformations at once
function Base.setproperty!(P::NFMT{D}, v::Symbol, val) where {D}
    if v == :x
        P.x_raw = copy(val)
        if P.backend !== nothing
            P.backend.x = _nfmt_transform_x(P, P.x_raw)
        end
    elseif v == :fhat
        want_complex = eltype(val) <: Complex
        _nfmt_ensure_backend!(P, want_complex)
        nf = _nfmt_reduced_size(P.N, P.basis_vect)
        if size(val)[1] != nf
            error("fhat has to be a vector of size $nf.")
        end
        if P.backend isa NFFT{D}
            P.backend.fhat = _nfmt_expand_reduced_to_full(
                ComplexF64.(val),
                P.N,
                P.basis_vect,
            )
        else
            scale = _nfmt_cos_alg_scale(P.N, P.basis_vect)
            P.backend.fhat = Float64.(val) .* scale
        end
    elseif v == :f
        want_complex = eltype(val) <: Complex
        _nfmt_ensure_backend!(P, want_complex)
        if size(val)[1] != P.M
            error("f has to be a vector of size M = $(P.M).")
        end
        if P.backend isa NFFT{D}
            P.backend.f = ComplexF64.(val)
        else
            P.backend.f = Float64.(val)
        end
    elseif v == :backend
        @warn "You can't modify the backend plan directly."
    elseif v == :basis_vect
        @warn "You can't modify the bases, please create an additional plan."
    elseif v == :N || v == :M || v == :n || v == :m || v == :f1 || v == :f2
        @warn "You can't modify this parameter, please create an additional plan."
    else
        Core.setfield!(P, v, val)
    end
end

# overwrite dot notation to get values from the struct and do the transformations at once
function Base.getproperty(P::NFMT{D}, v::Symbol) where {D}
    if v == :x
        return Core.getfield(P, :x_raw)
    elseif v == :num_threads
        if P.backend === nothing
            error("NFMT not initialized.")
        end
        return P.backend.num_threads
    elseif v == :f
        if P.backend === nothing
            error("f is not set.")
        end
        return P.backend.f
    elseif v == :fhat
        if P.backend === nothing
            error("fhat is not set.")
        end
        if P.backend isa NFFT{D}
            return _nfmt_reduce_full_to_reduced(P.backend.fhat, P.N, P.basis_vect)
        else
            return P.backend.fhat .* _nfmt_cos_alg_scale(P.N, P.basis_vect)
        end
    else
        return Core.getfield(P, v)
    end
end

@doc raw"""
    nfmt_trafo(P::NFMT{D})

computes the NFMT using the fast approximate algorithm (NFFT or NFCT, depending on the backend) for the provided nodes in `P.x` and coefficients in `P.fhat`.

# See also
[`NFMT{D}`](@ref), [`nfmt_trafo_direct`](@ref)
"""
function nfmt_trafo(P::NFMT{D}) where {D}
    if P.backend === nothing
        error("fhat has not been set.")
    end
    if P.backend isa NFFT{D}
        nfft_trafo(P.backend)
    else
        nfct_trafo(P.backend)
    end
    return P.f
end

function trafo(P::NFMT{D}) where {D}
    return nfmt_trafo(P)
end

@doc raw"""
    nfmt_adjoint(P::NFMT{D})

computes the adjoint NFMT using the fast approximate algorithm (NFFT or NFCT, depending on the backend) for the provided nodes in `P.x` and coefficients in `P.f`.

# See also
[`NFMT{D}`](@ref), [`nfmt_adjoint_direct`](@ref)
"""
function nfmt_adjoint(P::NFMT{D}) where {D}
    if P.backend === nothing
        error("f has not been set.")
    end
    if P.backend isa NFFT{D}
        nfft_adjoint(P.backend)
    else
        nfct_transposed(P.backend)
    end
    return P.fhat
end

function adjoint(P::NFMT{D}) where {D}
    return nfmt_adjoint(P)
end

@doc raw"""
    nfmt_trafo_direct(P::NFMT{D})

computes the NDMT via direct (naive) transformation for the provided nodes in `P.x` and coefficients in `P.fhat`.

# See also
[`NFMT{D}`](@ref), [`nfmt_trafo`](@ref)
"""
function nfmt_trafo_direct(P::NFMT{D}) where {D}
    if P.backend === nothing
        error("fhat has not been set.")
    end
    if P.backend isa NFFT{D}
        nfft_trafo_direct(P.backend)
    else
        nfct_trafo_direct(P.backend)
    end
    return P.f
end

function trafo_direct(P::NFMT{D}) where {D}
    return nfmt_trafo_direct(P)
end

@doc raw"""
    nfmt_adjoint_direct(P::NFMT{D})

computes the adjoint NDMT via direct (naive) transformation for the provided nodes in `P.x` and coefficients in `P.f`.

# See also
[`NFMT{D}`](@ref), [`nfmt_adjoint`](@ref)
"""
function nfmt_adjoint_direct(P::NFMT{D}) where {D}
    if P.backend === nothing
        error("f has not been set.")
    end
    if P.backend isa NFFT{D}
        nfft_adjoint_direct(P.backend)
    else
        nfct_transposed_direct(P.backend)
    end
    return P.fhat
end

function adjoint_direct(P::NFMT{D}) where {D}
    return nfmt_adjoint_direct(P)
end

@doc raw"""
    nfmt_get_LinearMap(basis_vect, N::Vector{Int}, X::Array{Float64}; n, m::Integer, f1::UInt32, f2::UInt32)::LinearMap

gives an linear map which computes the NFMT (always via the complex, NFFT-based path).

# See also
[`NFMT{D}`](@ref)
"""
function nfmt_get_LinearMap(
    basis_vect::NTuple{D,String},
    N::Vector{Int},
    X::Array{Float64};
    n = undef,
    m::Integer = 5,
    f1::UInt32 = (size(X, 1) > 1 ? f1_default : f1_default_1d),
    f2::UInt32 = f2_default,
)::LinearMap where {D}
    if size(X, 1) == 1
        X = vec(X)
        d = 1
        M = length(X)
    else
        (d, M) = size(X)
    end

    if N == []
        return LinearMap{ComplexF64}(fhat -> fill(fhat[1], M), f -> [sum(f)], M, 1)
    end

    b = copy(N)
    for (idx, s) in enumerate(basis_vect)
        if (BASES[s] > 0)
            b[idx] ÷= 2
        end
    end
    N = Tuple(N)

    if n == undef
        n = Tuple(2 * collect(N))
    end

    plan = NFMT(basis_vect, N, M, n, m, f1, f2)
    plan.x = X

    function trafo(fhat::Vector{ComplexF64})::Vector{ComplexF64}
        plan.fhat = fhat
        nfmt_trafo(plan)
        return plan.f
    end

    function adjoint(f::Vector{ComplexF64})::Vector{ComplexF64}
        plan.f = f
        nfmt_adjoint(plan)
        return plan.fhat
    end

    N = prod(b)
    return LinearMap{ComplexF64}(trafo, adjoint, M, N)
end

@doc raw"""
    nfmt_get_coefficient_vector(fhat::Array{ComplexF64})::Vector{ComplexF64}

reshapes an coefficient array to an vector for multiplication with the linear map of the NFMT.

# See also
[`NFMT{D}`](@ref), [`nfmt_get_LinearMap`](@ref)
"""
function nfmt_get_coefficient_vector(fhat::Array{ComplexF64})::Vector{ComplexF64}
    N = size(fhat)
    return vec(permutedims(fhat, length(N):-1:1))
end

@doc raw"""
    nfmt_get_coefficient_array(fhat::Vector{ComplexF64},P::NFMT{D})::Array{ComplexF64} where {D}

reshapes an coefficient vector returned from a linear map of the NFMT to an array.

# See also
[`NFMT{D}`](@ref), [`nfmt_get_LinearMap`](@ref)
"""
function nfmt_get_coefficient_array(
    fhat::Vector{ComplexF64},
    P::NFMT{D},
)::Array{ComplexF64} where {D}
    b = copy([P.N...])
    for (idx, s) in enumerate(P.basis_vect)
        if (BASES[s] > 0)
            b[idx] ÷= 2
        end
    end
    N = Tuple(b)
    return permutedims(reshape(fhat, reverse(N)), length(N):-1:1)
end

@doc raw"""
    nfmt_get_coefficient_array(fhat::Vector{ComplexF64},N::Vector{Int64},basis_vect)::Array{ComplexF64}

reshapes an coefficient vector returned from a linear map of the NFMT to an array.

# See also
[`NFMT{D}`](@ref), [`nfmt_get_LinearMap`](@ref)
"""
function nfmt_get_coefficient_array(
    fhat::Vector{ComplexF64},
    N::Vector{Int64},
    basis_vect::NTuple{D,String},
)::Array{ComplexF64} where {D}
    b = copy(N)
    for (idx, s) in enumerate(basis_vect)
        if (BASES[s] > 0)
            b[idx] ÷= 2
        end
    end
    N = Tuple(b)
    return permutedims(reshape(fhat, reverse(N)), length(N):-1:1)
end

@doc raw"""
    *(plan::NFMT{D}, fhat::Array{ComplexF64})::Vector{ComplexF64}

This function defines the multiplication of an NFMT plan with an coefficient array.
"""
function Base.:*(plan::NFMT{D}, fhat::Array{ComplexF64})::Vector{ComplexF64} where {D}
    if plan.x_raw === nothing
        error("x is not set.")
    end
    plan.fhat = nfmt_get_coefficient_vector(fhat)
    nfmt_trafo(plan)
    return plan.f
end

struct Adjoint_NFMT{D}
    plan::NFMT{D}
end

Adjoint_NFMT{D}(plan::NFMT{D}) where {D} = plan

@doc raw"""
    adjoint(plan::NFMT{D})::Adjoint_NFMT{D}

This function defines the adjoint operator for an NFMT.
"""
function Base.adjoint(plan::NFMT{D})::Adjoint_NFMT{D} where {D}
    return Adjoint_NFMT(plan)
end

@doc raw"""
    *(plan::Adjoint_NFMT{D}, f::Vector{ComplexF64})::Array{ComplexF64}

This function defines the multiplication of an adjoint NFMT plan with an vector of function values.
"""
function Base.:*(plan::Adjoint_NFMT{D}, f::Vector{ComplexF64})::Array{ComplexF64} where {D}
    if plan.plan.x_raw === nothing
        error("x is not set.")
    end
    plan.plan.f = f
    nfmt_adjoint(plan.plan)
    return nfmt_get_coefficient_array(plan.plan.fhat, plan.plan)
end
