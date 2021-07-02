# [Nonequispaced Fast Sine Transform (NFST)](@id NFST_site)

```@meta
    CurrentModule = NFFT3
```

## NDST and NFST

We consider the odd, 1-periodic trigonometric polynomial

```math
    f^s(\pmb{x}) \coloneqq \sum_{\pmb{k} \in I_{\pmb{N},\mathrm{s}}^d} \hat{f}_{\pmb{k}}^s \, \sin(2\pi \, \pmb{k} \odot \pmb{x}), \quad \pmb{x} \in \mathbb{R}^d,
```

with multibandlimit ``\pmb{N} \in \mathbb{N}^d`` and index set

```math
  I_{\pmb{N},\mathrm{s}}^d \coloneqq \left\{ \pmb{k} \in \mathbb{Z}^d: 1 \leq k_i \leq N_i - 1, \, i = 1,2,\ldots,d \right\}.
```

Note that we define ``\sin(\pmb{k} \circ \pmb{x}) \coloneqq \prod_{i=1}^d \sin(k_i \cdot x_i)``. The NDST is the evaluation of

```math
  f^s(\pmb{x}_j) = \sum_{\pmb{k} \in I_{\pmb{N},\mathrm{s}}^d} \hat{f}_{\pmb{k}}^s \, \sin(2\pi \, \pmb{k} \odot \pmb{x}_j)
```
at arbitrary nodes ``\pmb{x}_j \in [0,0.5]^d`` for given coefficients ``\hat{f}_{\pmb{k}}^s \in \mathbb{R}, \pmb{k} \in I_{\pmb{N},\mathrm{s}}^d``. Similarly to the NDFT, the transposed NDST is the evaluation of

```math
  \hat{h}^s_{\pmb{k}} = \sum_{j=1}^M f^s_j \, \sin(2\pi \, \pmb{k} \odot \pmb{x}_j)
```

for the frequencies ``\pmb{k} \in I_{\pmb{N},\mathrm{s}}^d`` with given coefficients ``f^s_j \in \mathbb{R}, j = 1,2,\ldots,M``.

We modify the [NFFT](@ref NFFT_site) in order to derive a fast algorithm for the computation of the NDST and transposed NDST, obtaining the NFST and its transposed counterpart. For details we refer to Chapter 6 in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)].

## Plan structure

```@docs
    NFST{D}
```

## Functions

```@docs
    nfst_trafo
    nfst_trafo_direct
    nfst_transposed
    nfst_transposed_direct
  	nfst_finalize_plan
    nfst_init
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