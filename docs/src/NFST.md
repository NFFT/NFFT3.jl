# [Nonequispaced Fast Sine Transform](@id NFST)

```@meta
    CurrentModule = NFFT3
```

```@docs
    NFST{D}
```

# NFST

## Pseudocode

**Input:** ``N, M \in \mathbb{N}, \, \sigma > 1, \, m \in \mathbb{N}, \, x_j \in (0, \pi) \text{ for } j = 0, \ldots, M - 1, \, f^{s}_k \in \mathbb{R} \text{ for } k = 0, \ldots, N - 1.``

**Precomputation:** Compute the nonzero Fourier coefficients ``c_k(\tilde{\varphi}) \text{ for all } k = 0, \ldots, N - 1.`` Compute the values ``\tilde{\psi} (x_j -\frac{\pi \ell}{\sigma N}) \text{ for } j = 0, \ldots, M - 1 \text{ and } \ell \in I^{\mathrm{T}}_{\sigma N, \ m} (x_j).``

1. Set ``\hat{g}_k \coloneqq \begin{cases} \frac{\hat{f}^{s}_k}{2 c_k(\tilde{\varphi})} & k = 1, \ldots, N - 1, \\ 0 & k = 0 \text{ and } k = N, \ldots, \sigma N .\end{cases}`` 
2. Compute `` g_{\ell} = \frac{1}{\sigma N} \sum_{k=1}^{\sigma N - 1} \hat{g}_k \sin ( \frac{\pi k \ell}{\sigma N} ), \quad \ell = 1, \ldots, \sigma N - 1`` using a fast algorithm of DST-I(``\sigma N - 1``), see  Table 6.1 and Remark 6.40 in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)] and set ``g_0 \coloneqq 0``. 
3. Compute `` s(x_j) \coloneqq \sum_{\ell=\lfloor 2\sigma N x_j \rfloor -m}^{\lceil 2 \sigma N x_j  \rceil +m} g_\ell\ \tilde{\psi}(x_j -\frac{\pi \ell}{\sigma N}), \quad j = 0, \ldots, M-1``.

**Output:** ``s(x_j), \, j = 0, \ldots, M-1``, approximate values of ``f^s(x_j)``.

**Computational cost:** ``\mathcal{O}(N \log{N} + m \ M)`` 

## Functions

```@docs
  	nfst_finalize_plan
    nfst_init
    nfst_trafo
    nfst_adjoint
    nfst_trafo_direct
    nfst_adjoint_direct
```