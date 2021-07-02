# [Nonequispaced Fast Sine Transform (NFST)](@id NFST_site)

```@meta
    CurrentModule = NFFT3
```

## NFST algorithm

We modify the [NFFT](@ref NFFT_site) in order to derive a fast algorithm for the evaluation of the odd, ``2 \pi``-periodic trigonometric polynomial

```math
    f^s (x) \coloneqq \sum^{N-1}_{k=1} \hat{f}_{k}^s \, \sin( k \ x), \quad x \in \mathbb{R}
```

at nonequispaced nodes ``x_j \in (0,\pi)``. To this end, we rewrite ``f^s`` as a sum of exponentials and obtain

```math
    \mathrm{i} \ f^s (x) = f(x) = \sum^{N-1}_{k=-N} \hat{f}_{k} \, \mathrm{e}^{\mathrm{i} \ k \ x} = \mathrm{i} \sum^{N-1}_{k=1} 2 \hat{f}_{k} \, \sin( k \ x), \quad x \in \mathbb{R}
```

with ``\hat{f}_0 = \hat{f}_{-N} = 0`` and ``\hat{f}_k = -\hat{f}_{-k} = \frac{1}{2}\hat{f}^{s}_k`` for ``k = 1, \ldots, N - 1``. Similarly as before, we approximate ``f(x)`` by a function ``s_1(x)`` and obtain for the coefficients ``g_{\ell}`` for ``\ell = 1, \ldots, \sigma N - 1``

```math
    - \mathrm{i} g_{\ell} = \frac{- \mathrm{i}}{2 \sigma N} \sum^{\sigma N - 1}_{k=-\sigma N} \hat{g}_{k} \, \mathrm{e}^{\pi \ \mathrm{i} \ k \ \ell / (\sigma N)} = \frac{1}{\sigma N} \sum^{\sigma N - 1}_{k=1} \hat{g}_{k} \, \sin( \frac{\pi \ k \ \ell}{\sigma N})
```

and particularly ``g_0 = g_{\sigma N} = 0``. Moreover, we observe that ``g_{2 \sigma N r - \ell} = -g_{\ell}`` for all ``r \in \mathbb{Z}``. Finally we compute the sum

```math
  	\mathrm{i} \ s(x_j) \coloneqq \sum_{\ell = \lfloor 2 \sigma N x \rfloor - m }^{\lceil 2 \sigma N x \rceil + m} \mathrm{i} \ g_{\ell} \ \tilde{\psi}(x_j - \frac{\pi \ell}{\sigma N})
```

and obtain the approximate values of ``f^s(x_j) = \mathrm{i} \ f(x_j) \approx \mathrm{i} \ s(x_j), \, j = 0, \ldots, M - 1``.

### Pseudocode

**Input:** ``N, M \in \mathbb{N}, \, \sigma > 1, \, m \in \mathbb{N}, \, x_j \in (0, \pi) \text{ for } j = 0, \ldots, M - 1, \, f^{s}_k \in \mathbb{R} \text{ for } k = 0, \ldots, N - 1.``

**Precomputation:** Compute the nonzero Fourier coefficients ``c_k(\tilde{\varphi}) \text{ for all } k = 0, \ldots, N - 1.`` Compute the values ``\tilde{\psi} (x_j -\frac{\pi \ell}{\sigma N}) \text{ for } j = 0, \ldots, M - 1 \text{ and } \ell \in I^{\mathrm{T}}_{\sigma N, \ m} (x_j).``

1. Set ``\hat{g}_k \coloneqq \begin{cases} \frac{\hat{f}^{s}_k}{2 c_k(\tilde{\varphi})} & k = 1, \ldots, N - 1, \\ 0 & k = 0 \text{ and } k = N, \ldots, \sigma N .\end{cases}`` 
2. Compute `` g_{\ell} = \frac{1}{\sigma N} \sum_{k=1}^{\sigma N - 1} \hat{g}_k \sin ( \frac{\pi k \ell}{\sigma N} ), \quad \ell = 1, \ldots, \sigma N - 1`` using a fast algorithm of DST-I(``\sigma N - 1``), see  Table 6.1 and Remark 6.40 in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)] and set ``g_0 \coloneqq 0``. 
3. Compute `` s(x_j) \coloneqq \sum_{\ell=\lfloor 2\sigma N x_j \rfloor -m}^{\lceil 2 \sigma N x_j  \rceil +m} g_\ell\ \tilde{\psi}(x_j -\frac{\pi \ell}{\sigma N}), \quad j = 0, \ldots, M-1``.

**Output:** ``s(x_j), \, j = 0, \ldots, M-1``, approximate values of ``f^s(x_j)``.

**Computational cost:** ``\mathcal{O}(N \log{N} + m \ M)`` 

### Transposed algorithm

The transposed problem reads as

```math
	h(k) \coloneqq \sum_{ j \in I_M^l} f_{j}\ \sin( k \ x_j), \quad  k \in I_{ N}^d = \{ k \in \mathbb{N}^d \colon 0 \leq k_i \leq N_i  \}
```

for given knots ``{x}_k \in [ 0,\pi ]^d, \, k=0,\ldots,M-1``, and coefficients ``f_j \in \mathbb{C}, j \in I_M^l``.
The algorithm for the fast evaluation of 

```math
	h(k) \coloneqq \sum^{M-1}_{ j = 0} h_{j} \ \sin(k \ x_j), \quad k = 0, \ldots, N-1,
```

with nonequispaced nodes ``{x}_j \in [ 0,\pi ], \, j=0,\ldots,M-1,`` can easily be derived from the duality of the two problems.

## Plan structure

```@docs
    NFST{D}
```

## Functions

```@docs
  	nfst_finalize_plan
    nfst_init
    nfst_trafo
    nfst_adjoint
    nfst_trafo_direct
    nfst_adjoint_direct
```