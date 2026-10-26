# [Fast Spherical Fourier Transform on an equispaced grid (FSFT)](@id FSFT_site)

```@meta
    CurrentModule = NFFT3
```

## FSFT

FSFT is the equispaced variant of [`NFSFT`](@ref NFSFT_site): the nodes are not user-supplied, they are the implicit equispaced grid of size ``(2N+2)(N+2)`` that the underlying C library always uses for this case - the `M` argument to the constructor is therefore ignored (`FSFT` computes the real node count itself).

## Plan structure

```@docs
    FSFT
```

## Functions

```@docs
    fsft_trafo
    fsft_trafo_direct
    fsft_adjoint
    fsft_adjoint_direct
    fsft_finalize_plan
    fsft_init
    fsft_index
```
