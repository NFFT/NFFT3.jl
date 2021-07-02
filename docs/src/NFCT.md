# [Nonequispaced Fast Cosine Transform (NFCT)](@id NFCT_site)

```@meta
    CurrentModule = NFFT3
```

## NDCT and NFCT

We consider the even, 1-periodic trigonometric polynomial

```math
    f^c(\pmb{x}) \coloneqq \sum_{\pmb{k} \in I_{\pmb{N},\mathrm{c}}^d} \hat{f}_{\pmb{k}}^c \, \cos(2\pi \, \pmb{k} \odot \pmb{x}), \quad \pmb{x} \in \mathbb{R}^d,
```

with multibandlimit ``\pmb{N} \in \mathbb{N}^d`` and index set

```math
  I_{\pmb{N},\mathrm{c}}^d \coloneqq \left\{ \pmb{k} \in \mathbb{Z}^d: 0 \leq k_i \leq N_i - 1, \, i = 1,2,\ldots,d \right\}.
```

Note that we define ``\cos(\pmb{k} \circ \pmb{x}) \coloneqq \prod_{i=1}^d \cos(k_i \cdot x_i)``. The NDCT is the evaluation of

```math
  f^c(\pmb{x}_j) = \sum_{\pmb{k} \in I_{\pmb{N},\mathrm{c}}^d} \hat{f}_{\pmb{k}}^c \, \cos(2\pi \, \pmb{k} \odot \pmb{x}_j)
```
at arbitrary nodes ``\pmb{x}_j \in [0,0.5]^d`` for given coefficients ``\hat{f}_{\pmb{k}}^c \in \mathbb{R}, \pmb{k} \in I_{\pmb{N},\mathrm{c}}^d``. Similarly to the NDFT, the transposed NDCT is the evaluation of

```math
  \hat{h}^c_{\pmb{k}} = \sum_{j=1}^M f^c_j \, \cos(2\pi \, \pmb{k} \odot \pmb{x}_j)
```

for the frequencies ``\pmb{k} \in I_{\pmb{N},\mathrm{c}}^d`` with given coefficients ``f^c_j \in \mathbb{R}, j = 1,2,\ldots,M``.

We modify the [NFFT](@ref NFFT_site) in order to derive a fast algorithm for the computation of the NDCT and transposed NDCT, obtaining the NFCT and its transposed counterpart. For details we refer to Chapter 7 in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)].

## Plan structure

```@docs
    NFCT{D}
```

## Functions

```@docs
    nfct_trafo
    nfct_adjoint
    nfct_trafo_direct
    nfct_adjoint_direct
  	nfct_finalize_plan
    nfct_init
```

## Literature

```@raw html
<ul>
<li id="PlonkaPottsSteidlTasche2018">[<a>Plonka, Potts, Steidl, Tasche, 2018</a>]
  G. Plonka, D. Potts, G. Steidl and M. Tasche. Numerical Fourier Analysis: Theory and Applications.</emph>
  Springer Nature Switzerland AG, 2018.
  doi: <a href="https://doi.org/10.1007/978-3-030-04306-3">10.1007/978-3-030-04306-3</a>.
</li>
</ul>
```