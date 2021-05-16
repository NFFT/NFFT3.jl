# [Nonequispaced Fast Fourier Transform](@id NFFT)

```@meta
    CurrentModule = NFFT3
```

```@docs
    NFFT{D}
```

# NFFT

The nonequispaced fast Fourier transform [[Keiner, Kunis, Potts, 2006](#KeinerKunisPotts2006)] (NFFT or NUFFT) overcomes one of the main shortcomings of the FFT - the need for an equispaced sampling grid. Considering a d-dimensional trigonometric polynomial

## Pseudocode

**Input:** ``N, M \in \mathbb{N}, \, \sigma > 1, \, m \in \mathbb{N}, \, x_j \in \mathbb{T}^d \text{ for } j \in I^{1}_{M}, \, f_k \in \mathbb{C} \text{ for } k \in I^{d}_{N}.``

**Precomputation:** Compute the nonzero Fourier coefficients ``c_k(\tilde{\varphi}) \text{ for all } k \in I^{d}_{N}.`` Compute the values ``\tilde{\psi} (x_j -\frac{\pi \pmb{l}}{\sigma N}) \text{ for } j \in I^{1}_{M} \text{ and } \pmb{l} \in I_{\sigma N, \ m} (x_j).``

1. Set ``\hat{g}_k \coloneqq \hat{f}_k / c_k(\tilde{\varphi}) \text{ for } k \in I^{d}_{N}.`` 
2. Compute `` g_{\pmb{l}} = \frac{1}{(\sigma N)^d} \sum_{k \in I^{d}_{N}} \hat{g}_k e^{2 \mathrm{i} \ k \cdot \pmb{l} / (\sigma N)}, \quad \pmb{l} \in I^{d}_{\sigma N}`` using a d-variate FFT. 
3. Compute ``s(x_j) \coloneqq \sum_{\pmb{l} \in I_{\sigma N, \ m} (x_j)} g_{\pmb{l}} \ \tilde{\psi}(x_j -\frac{\pi \pmb{l}}{\sigma N}), \quad \in I^{1}_{M}``.

**Output:** ``s(x_j), \, j \in I^{1}_{M}``, approximate values of ``f(x_j)``.

**Computational cost:** ``\mathcal{O}(N^d \log{N} + m^d \ M)`` 

## Functions

```@docs
  	nfft_finalize_plan
    nfft_init
    nfft_trafo
    nfft_adjoint
    nfft_trafo_direct
    nfft_adjoint_direct
```

## Literature

```@raw html
<ul>
<li id="KeinerKunisPotts2006">[<a>Keiner, Kunis, Potts, 2006</a>]
  J. Keiner, S. Kunis, and D. Potts. Fast summation of radial functions on the sphere. Computing, 78:1–15, 2006.
</li>
</ul>
```