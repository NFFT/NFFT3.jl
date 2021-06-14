# [Nonequispaced Fast Fourier Transform](@id NFFT_site)

```@meta
    CurrentModule = NFFT3
```

```@docs
    NFFT{D}
```

# NFFT

The nonequispaced fast Fourier transform [[Keiner, Kunis, Potts, 2006](#KeinerKunisPotts2006)] (NFFT or NUFFT) overcomes one of the main shortcomings of the FFT - the need for an equispaced sampling grid. Considering the evaluation of the ``d``-dimensional trigonometric polynomial

```math
  f \colon \mathbb{T}^d \to \mathbb{C}, \; \pmb{x} \mapsto \sum_{\pmb{k} \in I_{N}} \hat{f}_{\pmb{k}} \mathrm{e}^{-2 \pi \mathrm{i} \ \pmb{k} \cdot \pmb{x}}
```
 
with multibandlimit ``\pmb{N} \in 2 \mathbb{N}^d`` and index set

```math
  I_{N} \coloneqq \{ \pmb{k} \in \mathbb{Z}^d: - \frac{N_i}{2} \leq k_i \leq \frac{N_i}{2} - 1 \; \forall i \in \{ 1,2,\ldots, d \} \}.
```

The first approximation is a linear combination of a shifted periodized window function ``\tilde{\varphi}``

```math
  s_1(\pmb{x}) = \sum_{\pmb{\ell} \in I_{\pmb{n}}} g_{\pmb{\ell}} \tilde{\varphi} \left( \pmb{x} - \frac{1}{\pmb{n}} \odot \pmb{\ell} \right)
```

where ``\frac{1}{\pmb{n}}`` is the elementwise inversion of the vector ``\pmb{n}``. We choose an oversampling vector ``\pmb{\sigma} > 1`` componentwise and obtain the index set by

```math
  \pmb{n} \coloneqq \pmb{\sigma} \odot \pmb{N}.
```

Here, ``\odot`` denotes the componentwise product. Note that one could skip the choice of ``\pmb{\sigma}`` entirely and directly choose ``n_i > N_i`` for ``i= 1, 2, \ldots,d``. The standard choice in the C library is at the moment 

```math
  n_i = 2^{\lceil \log_2 N_i \rceil + 1}, i= 1, 2, \ldots, d. 
```

If the window function ``\varphi: \R^d \to \R`` has a one-periodization 

```math
  \tilde{\varphi}(\pmb{x}) = \sum_{\pmb{r} \in \mathbb{Z}^d} \varphi(\pmb{x}+\pmb{r}) 
```

with a uniformly convergent Fourier series, we use the Poisson summation formula and obtain the Fourier coefficients 

```math
  c_{\pmb{k}}(\tilde{\varphi}) = \hat{\varphi}(\pmb{k}).
```

Here, ``\hat{\varphi}`` is the Fourier transform of ``\varphi`` defined by

```math
  \hat{\varphi}(\pmb{k}) = \int_{\mathbb{T}^d} \varphi(\pmb{x}) \mathrm{e}^{-2 \pi \mathrm{i} \ \pmb{k} \cdot \pmb{x}} \mathrm{d} \pmb{x}.
```
Replacing ``\tilde{\varphi}`` by its Fourier series and splitting the sum in ``s_1(x)`` yields


```math
\begin{aligned}
  s_1(\pmb{x}) & = \sum_{\pmb{\ell} \in I_{\pmb{n}}} g_{\pmb{\ell}} \sum_{\pmb{k} \in \mathbb{Z}^d} c_{\pmb{k}}(\tilde{\varphi}) \mathrm{e}^{-2 \pi \mathrm{i} \ \pmb{k} \cdot \left( \pmb{x} - \frac{1}{\pmb{n}} \odot \pmb{\ell}\right) } \\ 
	& = \sum_{\pmb{k} \in \mathbb{Z}^d} c_{\pmb{k}} (\tilde{\varphi}) \underbrace{ \left(\sum_{\pmb{\ell} \in I_{\pmb{n}}} g_{\pmb{\ell}} \mathrm{e}^{2 \pi \mathrm{i} \ \frac{1}{\pmb{n}} \odot (\pmb{k} \cdot \pmb{\ell})} \right)}_{ \eqqcolon \hat{g}_{\pmb{k}} \mathrm{e}^{-2 \pi \mathrm{i} \ \pmb{k} \cdot \pmb{x}}} \\
	& = \sum_{\pmb{k} \in I_{\pmb{n}}} c_{\pmb{k}}(\tilde{\varphi}) \hat{g}_{\pmb{k}} \mathrm{e}^{-2 \pi \mathrm{i} \ \pmb{k} \cdot \pmb{x} } + \sum_{\pmb{r} \in \mathbb{Z}^d \setminus \{ \pmb{0} \}} \sum_{\pmb{k} \in I_{\pmb{n}}} c_{\pmb{k}}(\tilde{\varphi}) \hat{g}_{\pmb{k}} \mathrm{e}^{-2 \pi \mathrm{i} \ (\pmb{k} + \pmb{n} \odot \pmb{r})\cdot \pmb{x} }.
\end{aligned}
```

Furthermore, a comparison of the equation above and ``f`` suggests the choice

```math
  \hat{g}_{\pmb{k}} = \begin{cases} \frac{\hat{f}_{\pmb{k}}}{c_{\pmb{k}}(\tilde{\varphi})} \quad &\text{for } \pmb{k} \in I_{\pmb{N}} \\ 0 &\text{for } \pmb{k} \in I_{\pmb{n}} \setminus I_{\pmb{N}} \end{cases}.
```

Now, we are able to compute ``g_{\pmb{\ell}}`` by an FFT of size ``n_1 \times n_2 \times \cdots \times n_d``. The error of the first approximation ``f \approx s_1`` is called aliasing error. For further analysis, we refer to Chapter 7 in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)].
\
We obtain an approximation ``\psi`` for ``\varphi`` if the window function is well localized in the spatial domain by the truncation

```math
  \psi(\pmb{x}) \coloneqq \varphi(\pmb{x}) \mathbb{1}_{\times_{i = 1}^d [-m/n_i,m/n_i]}(\pmb{x})
```

with a chosen window size ``m \ll \min_{i\in\{1,2,\dots,d\}}\{n_i\}, m \in \mathbb{N}``. Following the same scheme as above, we can use the periodization 

```math
  \tilde{\psi}(\pmb{x}) = \sum_{\pmb{r} \in \mathbb{Z}^d} \psi(\pmb{x}+\pmb{r}) 
``` 

and look at the linear combination 

```math
  s(\pmb{x}_j) \coloneqq \sum_{\pmb{\ell} \in I_{\pmb{n}}} g_{\pmb{\ell}} \tilde{\psi} \left( \pmb{x}_j - \frac{1}{\pmb{n}} \odot \pmb{\ell} \right).
```

The following calculations show that ``s \approx s_1``

```math
\begin{aligned}
  s(\pmb{x}_j) &= \sum_{\pmb{\ell} \in I_{\pmb{n}}} g_{\pmb{\ell}} \sum_{\pmb{r} \in \mathbb{Z}^d} \psi \left( \pmb{x}_j - \frac{1}{\pmb{n}} \odot \pmb{\ell} + \pmb{r} \right) \\
  &= \sum_{\pmb{\ell} \in I_{\pmb{n}}} g_{\pmb{\ell}} \sum_{\pmb{r} \in \mathbb{Z}^d} \varphi \left(\pmb{x}_j - \frac{1}{\pmb{n}} \odot \pmb{\ell} + \pmb{r} \right) \mathbb{1}_{\times_{i = 1}^d [-m/n_i,m/n_i]}(\pmb{x}_j) \\
  &= \sum_{\ell \in I_n} g_{\pmb{\ell}} \mathbb{1}_{\times_{i = 1}^d [-m/n_i,m/n_i]}(\pmb{x}_j) \tilde{\varphi} \left(\pmb{x}_j - \frac{1}{\pmb{n}} \odot \pmb{\ell} \right) \\
  &= \sum_{\pmb{\ell} \in I_{\pmb{n},m} (\pmb{x}_j)} g_{\pmb{\ell}} \tilde{\varphi} \left(\pmb{x}_j - \frac{1}{\pmb{n}} \odot \pmb{\ell} \right)
\end{aligned}
```

with the index set 

```math
  I_{\pmb{n},m}(\pmb{x}_j) = \{ \pmb{\ell} \in I_{\pmb{n}} \colon \pmb{n} \odot \pmb{x}_j - m \pmb{1} \leq \pmb{\ell} \leq \pmb{n} \odot \pmb{x}_j +m \pmb{1} \}
```

for a fixed node ``\pmb{x}_j``. This is motivated by 

```math
  -\frac{m}{n_i} \leq \left( \pmb{x}_j \right)_i \leq \frac{m}{n_i} 
```

in order to ensure that ``\pmb{x}_j`` is within the support. This second approximation error is called the truncation error. Summarizing, we have ``f \approx s_1 \approx s``.


### Pseudocode

**Input:** ``N, M \in \mathbb{N}, \, \sigma > 1, \, m \in \mathbb{N}, \, x_j \in \mathbb{T}^d \text{ for } j \in I^{1}_{M}, \, f_k \in \mathbb{C} \text{ for } k \in I^{d}_{N}.``

**Precomputation:** Compute the nonzero Fourier coefficients ``c_k(\tilde{\varphi}) \text{ for all } k \in I^{d}_{N}.`` Compute the values ``\tilde{\psi} (x_j -\frac{\pi \pmb{l}}{\sigma N}) \text{ for } j \in I^{1}_{M} \text{ and } \pmb{l} \in I_{\sigma N, \ m} (x_j).``

1. Set ``\hat{g}_k \coloneqq \hat{f}_k / c_k(\tilde{\varphi}) \text{ for } k \in I^{d}_{N}.`` 
2. Compute `` g_{\pmb{l}} = \frac{1}{(\sigma N)^d} \sum_{k \in I^{d}_{N}} \hat{g}_k e^{2 \mathrm{i} \ k \cdot \pmb{l} / (\sigma N)}, \quad \pmb{l} \in I^{d}_{\sigma N}`` using a d-variate FFT. 
3. Compute ``s(x_j) \coloneqq \sum_{\pmb{l} \in I_{\sigma N, \ m} (x_j)} g_{\pmb{l}} \ \tilde{\psi}(x_j -\frac{\pi \pmb{l}}{\sigma N}), \quad \in I^{1}_{M}``.

**Output:** ``s(x_j), \, j \in I^{1}_{M}``, approximate values of ``f(x_j)``.

**Computational cost:** ``\mathcal{O}(N^d \log{N} + m^d \ M)`` 

## Adjoint Problem

Using the transposed index set 

```math
  I_{\pmb{n},m}^\top(\pmb{\ell}) = \{ j= 0, 1, \ldots, M-1 : \pmb{\ell} - m\pmb{1} \leq \pmb{n} \odot \pmb{x}_j \leq \pmb{\ell} + m \pmb{1} \},
```

we obtain the adjoint NFFT algorithm for the fast evaluation the adjoint problem

```math
	\hat{h}_{\pmb{k}} = \sum_{j = 0}^{M-1} f_j \mathrm{e}^{2 \pi \ \mathrm{i} \ \pmb{k} \ \pmb{x}_j}, \pmb{k} \in I_{\pmb{N}}
```

for given coefficients ``f_j \in \mathbb{C}``.

# Functions

```@docs
  	nfft_finalize_plan
    nfft_init
    nfft_trafo
    nfft_adjoint
    nfft_trafo_direct
    nfft_adjoint_direct
```

### Literature

```@raw html
<ul>
<li id="KeinerKunisPotts2006">[<a>Keiner, Kunis, Potts, 2006</a>]
  J. Keiner, S. Kunis, and D. Potts. Fast summation of radial functions on the sphere. Computing, 78:1–15, 2006.
</li>
</ul>
```

```@raw html
<ul>
<li id="PlonkaPottsSteidlTasche2018">[<a>Plonka, Potts, Steidl, Tasche, 2018</a>]
  G. Plonka, D. Potts, G. Steidl and M. Tasche. Numerical Fourier Analysis: Theory and Applications.</emph>
  Springer Nature Switzerland AG, 2018.
  doi: <a href="https://doi.org/10.1007/978-3-030-04306-3">10.1007/978-3-030-04306-3</a>.
</li>
</ul>
```