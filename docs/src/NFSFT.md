# [Nonequispaced Fast Spherical Fourier Transform (NFSFT)](@id NFSFT_site)

```@meta
    CurrentModule = NFFT3
```

## NDSFT and NFSFT

We consider spherical Fourier coefficients ``\hat{f}_k^n \in \mathbb{C}`` of degree ``k`` and order ``n``, ``0 \leq k \leq N``, ``-k \leq n \leq k``, and a bandwidth ``N \in \mathbb{N}``. The NDSFT is the evaluation of

```math
f(\pmb{x}_j) = \sum_{k=0}^N \sum_{n=-k}^k \hat{f}_k^n \, Y_k^n(\pmb{x}_j)
```

at given arbitrary nodes ``\pmb{x}_j`` on the unit sphere, ``j = 1, \ldots, M``, where ``Y_k^n`` are the spherical harmonics of degree ``k`` and order ``n``. The transposed NDSFT is the evaluation of

```math
\hat{h}_k^n = \sum_{j=1}^M f_j \, \overline{Y_k^n(\pmb{x}_j)}
```

for given coefficients ``f_j \in \mathbb{C}, j = 1, \ldots, M``.

The NFSFT realizes a fast algorithm for the NDSFT and the transposed NDSFT, based on the [NFFT](@ref NFFT_site) and the fast polynomial transform (FPT). A node ``\pmb{x} = (x_0, x_1)`` encodes spherical coordinates via ``x_0 = \varphi/(2\pi) \in [-0.5, 0.5)`` (azimuth) and ``x_1 = \vartheta/\pi \in [0, 1)`` (polar angle).

See [`FSFT`](@ref FSFT_site) for the variant on an equispaced grid.

## Plan structure

```@docs
    NFSFT
```

## Functions

```@docs
    nfsft_trafo
    nfsft_trafo_direct
    nfsft_adjoint
    nfsft_adjoint_direct
    nfsft_finalize_plan
    nfsft_init
    nfsft_index
```
