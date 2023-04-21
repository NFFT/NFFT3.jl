using NFFT3

mutable struct NFFCT{D}
    dcos::NTuple{D,String}             # dimensions with cos
    NFFT_struct::NFFT{D}             # structure for the NFFT
    function NFFCT{D}(
        dcos::NTuple{D,String},
        P::NFFT{D}
    ) where {D}
        new(dcos, P)
    end
end

function NFFCT(
    dcos::NTuple{D,String},
    N::NTuple{D,Integer},
    M::Integer,
    n::NTuple{D,Integer},
    m::Integer,
    f1::UInt32 = (D > 1 ? NFFT3.f1_default : NFFT3.f1_default_1d),
    f2::UInt32 = NFFT3.f2_default,
) where {D}
    NFFCT{D}(dcos, NFFT{D}(
        NTuple{D,Int32}(N),
        Int32(M),
        NTuple{D,Int32}(n),
        Int32(m),
        (f1 | NFFT3.MALLOC_X | NFFT3.MALLOC_F_HAT | NFFT3.MALLOC_F | NFFT3.FFTW_INIT),
        f2,
    ))
end

# additional constructor for easy use [NFFCT(dcos,N,M) instead of NFFCT{2}(dcos,N,M)]
function NFFCT(dcos::NTuple{D,String},N::NTuple{D,Integer}, M::Integer) where {D}
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
    NFFCT{D}(dcos, NFFT{D}(NTuple{D,Int32}(N), Int32(M), n, Int32(8), f1, NFFT3.f2_default))
end

@doc raw"""
    nffct_finalize_plan(P::NFFCT{D})

destroys a NFFT plan structure.

# Input
* `P` - a NFFCT structure.

# See also
[`NFFT{D}`](@ref), [`nfft_init`](@ref)
"""
function nffct_finalize_plan(P::NFFCT{D}) where {D}
    nfft_finalize_plan(P.NFFT_struct)
end

function finalize_plan(P::NFFCT{D}) where {D}
    return nffct_finalize_plan(P)
end

# overwrite dot notation for struct to set values and do the transformations at once
function Base.setproperty!(P::NFFCT{D}, v::Symbol, val) where {D}
    if v == :x
        xh = copy(val)
        if D==1
            if (BASES[P.dcos[1]]==1)
                xh ./= 2
            elseif (BASES[P.dcos[1]]==2)
                xh = acos(2*xh-1)/(2*pi)
            end
        else
            for i in range(1, D)
                if (BASES[P.dcos[i]]==1)
                    xh[i,:] ./= 2
                elseif (BASES[P.dcos[i]]==2)
                    xh[i,:] = acos.(2 .*xh[i,:].-1)./(2*pi)
                end
            end
        end
        P.NFFT_struct.x = xh
    elseif v == :fhat
        a = sum(getindex.([BASES],P.dcos).>0)
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
            if (BASES[P.dcos[i+1]]>0)
                pv[i] ÷= 2;
            end
        end
        fhatexp = zeros(Complex, p)
        for i = 1:p
            idx = 1
            e = 0;
            for j = 1:D
                if (BASES[P.dcos[j]]>0)
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
    elseif v == :dcos
        @warn "You can't modify the cosinus dimensions, please create an additional plan."
    else
        Base.setproperty!(P.NFFT_struct, v, val)
    end
end

# overwrite dot notation to get values from the struct and do the transformations at once
function Base.getproperty(P::NFFCT{D}, v::Symbol) where {D}
    if v == :NFFT
        return P.NFFT_struct
    elseif v == :f || v == :plan || v == :num_threads || v == :init_done || v == :N || v == :M || v == :n || v == :m || v == :f1 || v == :f2
        return Base.getproperty(P.NFFT_struct, v)
    elseif v == :x
        xd = copy(P.NFFT_struct.x)
        if D==1
            if (BASES[P.dcos[1]]==1)
                xd .*= 2
            elseif (BASES[P.dcos[1]]==2)
                xd = (cos(2*pi*xd)+1)/2
            end
        else
            for i in range(1, D)
                if (BASES[P.dcos[i]]==1)
                    xd[i,:] .*= 2
                elseif (BASES[P.dcos[i]]==2)
                    xd[i,:] = (cos.(2 .*pi .*xd[i,:]).+1)./2
                end
            end
        end
        return xd
    elseif v == :fhat
        a = sum(getindex.([BASES],P.dcos).>0)
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
            if (BASES[P.dcos[i+1]]>0)
                pv[i] ÷= 2;
            end
        end
        fhat = zeros(Complex, l)
        fhat_old = P.NFFT_struct.fhat
        for i = 1:p
            idx = 1
            e = 0;
            for j = 1:D
                if (BASES[P.dcos[j]]>0)
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

# nffct trafo [call with NFFT.trafo outside module]
@doc raw"""
    nffct_trafo(P)

computes the NDFCT via the fast NFFCT algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}}^c \in \mathbb{R}, \pmb{k} \in I_{\pmb{N},\mathrm{c}}^D,`` in `P.fhat`.

# Input
* `P` - a NFFCT plan structure.

# See also
[`NFFCT{D}`](@ref), [`nffct_trafo_direct`](@ref)
"""
function nffct_trafo(P::NFFCT{D}) where {D}
    return nfft_trafo(P.NFFT_struct)
end

function trafo(P::NFFCT{D}) where {D}
    return nffct_trafo(P)
end

# adjoint trafo [call with NFFT.adjoint outside module]
@doc raw"""
    nffct_transposed(P)

computes the transposed NDFCT via the fast transposed NFFCT algorithm for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j^c \in \mathbb{R}, j =1,2,\dots,M,`` in `P.f`.

# Input
* `P` - a NFFCT plan structure.

# See also
[`NFFCT{D}`](@ref), [`nffct_transposed_direct`](@ref)
"""
function nffct_adjoint(P::NFFCT{D}) where {D}
    return nfft_adjoint(P.NFFT_struct)
end

function adjoint(P::NFFCT{D}) where {D}
    return nffct_adjoint(P)
end

# nffct trafo direct [call with NFFCT.trafo_direct outside module]
@doc raw"""
    nffct_trafo_direct(P)

computes the NDFCT via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``\hat{f}_{\pmb{k}} \in \mathbb{C}, \pmb{k} \in I_{\pmb{N}}^D,`` in `P.fhat`.

# Input
* `P` - a NFFCT plan structure.

# See also
[`NFFCT{D}`](@ref), [`nffct_trafo`](@ref)
"""
function nffct_trafo_direct(P::NFFCT{D}) where {D}
    return nfft_trafo_direct(P.NFFT_struct)
end

function trafo_direct(P::NFFCT{D}) where {D}
    return nffct_trafo_direct(P)
end

# adjoint trafo direct [call with NFFCT.adjoint_direct outside module]
@doc raw"""
    nffct_adjoint_direct(P)

computes the adjoint NDFCT via naive matrix-vector multiplication for provided nodes ``\pmb{x}_j, j =1,2,\dots,M,`` in `P.X` and coefficients ``f_j \in \mathbb{C}, j =1,2,\dots,M,`` in `P.f`.

# Input
* `P` - a NFFCT plan structure.

# See also
[`NFFCT{D}`](@ref), [`nffct_adjoint`](@ref)
"""
function nffct_adjoint_direct(P::NFFCT{D}) where {D}
    return nfft_adjoint_direct(P.NFFT_struct)
end

function adjoint_direct(P::NFFCT{D}) where {D}
    return nffct_adjoint_direct(P)
end