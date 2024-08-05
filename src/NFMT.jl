using NFFT3

# NFMT plan struct
@doc raw"""
    NFMT{D}

A NFMT (nonequispaced fast mixed transform) plan, where D is the dimension. 

The NFCT realizes a direct and fast computation of the discrete nonequispaced mixed transform. The aim is to compute

$$f^{\pmb{d}}(\pmb{x}_j) \coloneqq \sum_{\pmb{k} \in I_{\pmb{N},\pmb{d}}^d} \hat{f}_{\pmb{k}}^{\pmb{d}} \, \phi_{\pmb{k}}^{\pmb{d}}(\pmb{x}_j)$$

with 

$$\phi_{\pmb{k}}^{\pmb{d}}(\pmb{x})=\prod_{j=1}^d\begin{cases}1,&k_j=0\\\exp(2\pi\mathrm{i}k_jx_j),&d_j=\exp,\;k_j\neq0\\
\sqrt{2}\cos(\pi k_jx_j),&d_j=\cos,\;k_j\neq0\\
\sqrt{2}\cos(k_j\arccos(2x_j-1)),&d_j=\mathrm{alg},\;k_j\neq0\end{cases}$$

for $\pmb{d}\in\{\exp,\cos,\mathrm{alg}\}^d$ at given arbitrary knots ``\pmb{x}_j \in [0,1]^D, j = 1, \cdots, M``, for coefficients ``\hat{f}^{c}_{\pmb{k}} \in \mathbb{R}``, 

$$\pmb{k}\inI_{\pmb{N},\pmb{d}}^d \coloneqq \overset{d}{\underset{j=1}{\vphantom{\mathop{\raisebox{-.5ex}{\hbox{\huge{$\times$}}}}}⨉}}\begin{cases}\Big\{-\frac{N_j}{2},-\frac{N_j}{2}+1,\ldots,\frac{N_j}{2}\Big\},&d_j=\exp\\\Big\{0,1,\ldots,\frac{N_j}{2}\Big\},&d_j\neq\exp\end{cases},$$

and a multibandlimit vector ``\pmb{N} \in (2\mathbb{N})^{D}``. The transposed problem reads as

$$\hat{h}^{\pmb{d}}_{\pmb{k}} = \sum_{j=1}^M f^{\pmb{d}}_j \, \phi_{\pmb{k}}^{\pmb{d}}(\pmb{x}_j)$$

for the frequencies ``\pmb{k} \in I_{\pmb{N},\pmb{d}}^D`` with given coefficients ``f^{\pmb{d}}_j \in \mathbb{R}, j = 1,2,\ldots,M``.

# Fields
* `basis_vect` - vector with the bases (`exp`, `cos`, `alg`) for each dimension
* `NFFT_struct` - underlying [NFFT plan](@ref NFFT_site# Plan structure).

# Constructor
    NFMT{D}( basis_vect::NTuple{D,String}, P::NFFT{D}N::NTuple{D,Int32}) where {D}


> **WARNING**: Not for direct usage!

# Additional Constructor
    NFMT( N::NTuple{D,Int32}, M::Int32, n::NTuple{D,Int32}, m::Int32, f1::UInt32, f2::UInt32) where {D}
    NFMT( N::NTuple{D,Int32}, M::Int32) where {D}

with
* `N` - the multibandlimit ``(N_1, N_2, \ldots, N_D)`` of the trigonometric polynomial ``f^s``.
* `M` - the number of nodes.
* `n` - the oversampling ``(n_1, n_2, \ldots, n_D)`` per dimension.
* `m` - the window size. A larger m results in more accuracy but also a higher computational cost. 
* `f1` - the NFST flags.
* `f2` - the FFTW flags.

# See also
[`NFFT`](@ref)
"""
mutable struct NFMT{D}
    basis_vect::NTuple{D,String}             # dimensions with cos
    NFFT_struct::NFFT{D}             # structure for the NFFT
    function NFMT{D}(
        basis_vect::NTuple{D,String},
        P::NFFT{D}
    ) where {D}
        new(basis_vect, P)
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
    NFMT{D}(basis_vect, NFFT{D}(
        NTuple{D,Int32}(N),
        Int32(M),
        NTuple{D,Int32}(n),
        Int32(m),
        (f1 | NFFT3.MALLOC_X | NFFT3.MALLOC_F_HAT | NFFT3.MALLOC_F | NFFT3.FFTW_INIT),
        f2,
    ))
end

# additional constructor for easy use [NFMT(basis_vect,N,M) instead of NFMT{2}(basis_vect,N,M)]
function NFMT(basis_vect::NTuple{D,String},N::NTuple{D,Integer}, M::Integer) where {D}
    # convert N to vector for passing it over to C
    Nv = collect(N)

    # default oversampling
    n = Array{Int32}(2 .^ (ceil.(log.(Nv) / log(2)) .+ 1))
    n = NTuple{D,Int32}(n)

    # default NFFT flags
    f1 = UInt32(0)

    if D > 1
        f1 = NFFT3.f1_default
    else
        f1 = NFFT3.f1_default_1d
    end
    NFMT{D}(basis_vect, NFFT{D}(NTuple{D,Int32}(N), Int32(M), n, Int32(8), f1, NFFT3.f2_default))
end

@doc raw"""
    nfmt_finalize_plan(P::NFMT{D})

destroys a NFFT plan structure.

# Input
* `P` - a NFMT structure.

# See also
[`NFFT{D}`](@ref), [`nfft_init`](@ref)
"""
function nfmt_finalize_plan(P::NFMT{D}) where {D}
    nfft_finalize_plan(P.NFFT_struct)
end

function finalize_plan(P::NFMT{D}) where {D}
    return nfmt_finalize_plan(P)
end

# overwrite dot notation for struct to set values and do the transformations at once
function Base.setproperty!(P::NFMT{D}, v::Symbol, val) where {D}
    if v == :x
        xh = copy(val)
        if D==1
            if (BASES[P.basis_vect[1]]==1)
                xh ./= 2
            elseif (BASES[P.basis_vect[1]]==2)
                xh = acos.(2 .*xh.-1)./(2*pi)
            end
        else
            for i in range(1, D)
                if (BASES[P.basis_vect[i]]==1)
                    xh[i,:] ./= 2
                elseif (BASES[P.basis_vect[i]]==2)
                    xh[i,:] = acos.(2 .*xh[i,:].-1)./(2*pi)
                end
            end
        end
        P.NFFT_struct.x = xh
    elseif v == :fhat
        a = sum(getindex.([BASES],P.basis_vect).>0)
        p = prod(P.NFFT_struct.N)
        l = p÷(2^a)
        if size(val)[1] != l
            error("fhat has to be a ComplexFloat64 vector of size prod(N).")
        end
        pvExp = zeros(Int64, D)
        pvExp[D] = 1
        for i = D-1:-1:1
            pvExp[i] = pvExp[i+1] * P.NFFT_struct.N[i+1]
        end
        pv = zeros(Int64, D)
        pv[D] = 1
        for i = D-1:-1:1
            pv[i] = pv[i+1] * P.NFFT_struct.N[i+1]
            if (BASES[P.basis_vect[i+1]]>0)
                pv[i] ÷= 2;
            end
        end
        fhatexp = zeros(Complex, p)
        for i = 1:p
            idx = 1
            e = 0;
            for j = 1:D
                if (BASES[P.basis_vect[j]]>0)
                    k = abs(((i - 1) ÷ pvExp[j]) % P.NFFT_struct.N[j]-P.NFFT_struct.N[j] ÷ 2)
                    if k == P.NFFT_struct.N[j] ÷ 2
                        idx = -1
                        break
                    end
                    if k ≠ 0
                        e -= 1
                    end
                    idx += k * pv[j]
                else
                    idx += (((i - 1) ÷ pvExp[j]) % P.NFFT_struct.N[j]) * pv[j]
                end
            end
            if idx == -1
                fhatexp[i] = 0
            else
                fhatexp[i] = sqrt(2)^e * val[idx]
            end
        end
        P.NFFT_struct.fhat = fhatexp
    elseif v == :f
        P.NFFT_struct.f = val
    elseif v == :NFFT_struct
        @warn "You can't modify the pointer to the NFFT plan."
    elseif v == :basis_vect
        @warn "You can't modify the cosinus dimensions, please create an additional plan."
    else
        Base.setproperty!(P.NFFT_struct, v, val)
    end
end

# overwrite dot notation to get values from the struct and do the transformations at once
function Base.getproperty(P::NFMT{D}, v::Symbol) where {D}
    if v == :NFFT
        return P.NFFT_struct
    elseif v == :f || v == :plan || v == :num_threads || v == :init_done || v == :N || v == :M || v == :n || v == :m || v == :f1 || v == :f2
        return Base.getproperty(P.NFFT_struct, v)
    elseif v == :x
        xd = copy(P.NFFT_struct.x)
        if D==1
            if (BASES[P.basis_vect[1]]==1)
                xd .*= 2
            elseif (BASES[P.basis_vect[1]]==2)
                xd = (cos.(2 .*pi .*xd).+1)./2
            end
        else
            for i in range(1, D)
                if (BASES[P.basis_vect[i]]==1)
                    xd[i,:] .*= 2
                elseif (BASES[P.basis_vect[i]]==2)
                    xd[i,:] = (cos.(2 .*pi .*xd[i,:]).+1)./2
                end
            end
        end
        return xd
    elseif v == :fhat
        a = sum(getindex.([BASES],P.basis_vect).>0)
        p = prod(P.NFFT_struct.N)
        l = p÷(2^a)
        pvExp = zeros(Int64, D)
        pvExp[D] = 1
        for i = D-1:-1:1
            pvExp[i] = pvExp[i+1] * P.NFFT_struct.N[i+1]
        end
        pv = zeros(Int64, D)
        pv[D] = 1
        for i = D-1:-1:1
            pv[i] = pv[i+1] * P.NFFT_struct.N[i+1]
            if (BASES[P.basis_vect[i+1]]>0)
                pv[i] ÷= 2;
            end
        end
        fhat = zeros(Complex, l)
        fhat_old = P.NFFT_struct.fhat
        for i = 1:p
            idx = 1
            e = 0;
            for j = 1:D
                if (BASES[P.basis_vect[j]]>0)
                    k = abs(((i - 1) ÷ pvExp[j]) % P.NFFT_struct.N[j]-P.NFFT_struct.N[j] ÷ 2)
                    if k == P.NFFT_struct.N[j] ÷ 2
                        idx = -1
                        break
                    end
                    if k ≠ 0
                        e -= 1
                    end
                    idx += k * pv[j]
                else
                    idx += (((i - 1) ÷ pvExp[j]) % P.NFFT_struct.N[j]) * pv[j]
                end
            end
            if idx ≠ -1
                fhat[idx] += sqrt(2)^e * fhat_old[i]
            end
        end
        return fhat
    else
        return Core.getfield(P, v)
    end
end

# NFMT trafo [call with NFFT.trafo outside module]
@doc raw"""
    nfmt_trafo(P)

computes the NDMT via the fast NFMT algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}}^{\pmb{d}} \in \mathbb{R}, \pmb{k} \in I_{\pmb{N},\pmb{d}}^D,`` in `P.fhat`.

# Input
* `P` - a NFMT plan structure.

# See also
[`NFMT{D}`](@ref), [`nfmt_trafo_direct`](@ref)
"""
function nfmt_trafo(P::NFMT{D}) where {D}
    return nfft_trafo(P.NFFT_struct)
end

function trafo(P::NFMT{D}) where {D}
    return nfmt_trafo(P)
end

# adjoint trafo [call with NFFT.adjoint outside module]
@doc raw"""
    nfmt_transposed(P)

computes the transposed NDFCT via the fast transposed NFMT algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j^{\pmb{d}} \in \mathbb{R}, j =1,2,\dots,M,`` in `P.f`.

# Input
* `P` - a NFMT plan structure.

# See also
[`NFMT{D}`](@ref), [`nfmt_transposed_direct`](@ref)
"""
function nfmt_adjoint(P::NFMT{D}) where {D}
    return nfft_adjoint(P.NFFT_struct)
end

function adjoint(P::NFMT{D}) where {D}
    return nfmt_adjoint(P)
end

# NFMT trafo direct [call with NFMT.trafo_direct outside module]
@doc raw"""
    nfmt_trafo_direct(P)

computes the NDMT via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}} \in \mathbb{C}, \pmb{k} \in I_{\pmb{N}}^D,`` in `P.fhat`.

# Input
* `P` - a NFMT plan structure.

# See also
[`NFMT{D}`](@ref), [`nfmt_trafo`](@ref)
"""
function nfmt_trafo_direct(P::NFMT{D}) where {D}
    return nfft_trafo_direct(P.NFFT_struct)
end

function trafo_direct(P::NFMT{D}) where {D}
    return nfmt_trafo_direct(P)
end

# adjoint trafo direct [call with NFMT.adjoint_direct outside module]
@doc raw"""
    nfmt_adjoint_direct(P)

computes the adjoint NDMT via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j \in \mathbb{C}, j =1,2,\dots,M,`` in `P.f`.

# Input
* `P` - a NFMT plan structure.

# See also
[`NFMT{D}`](@ref), [`nfmt_adjoint`](@ref)
"""
function nfmt_adjoint_direct(P::NFMT{D}) where {D}
    return nfft_adjoint_direct(P.NFFT_struct)
end

function adjoint_direct(P::NFMT{D}) where {D}
    return nfmt_adjoint_direct(P)
end