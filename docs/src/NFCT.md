# [Nonequispaced Fast Cosine Transform](@id NFCT)

```@docs
    NFCT{D}
```

## Background

The NFCT (Nonequispaced fast cosine transform) ([[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)], Sec. 7.4) realizes a multivariate fast cosine transform for nonequispaced knots. The aim is to compute 

```math
    f^c (x) \coloneqq \sum_{ k \in I_{ N}} \hat{f}_{ k}^c \cos( k \cdot x), \quad x \in \mathbb{R}^d
```
at given (nonequidistant) knots `` x_{k} \in \left[ 0,\pi \right]^d, \ k=0,\ldots,M-1``, coefficients ``\hat{f_{ k}^c} \in \mathbb{R}, \  k\in I_{ N} \coloneqq \{ k\in \mathbb{N}^d:\ 0\leq k_i \leq N_i\ \forall\, i=1,\ldots,d \}`` for some multibandlimit vector `` N \in \mathbb{N}^d``. 
\
The transposed (adjoined) problem reads as

```math
	h(k) \coloneqq \sum_{ j\in I_M^l} f_{ j}\ \cos({k\,x_j}), \quad  k\in I_{ N}^d \coloneqq \{ k\in \mathbb{N}^d: 0\leq k_i\leq N_i  \}
```

for given knots ``{x}_k\in \left[ 0,\pi \right]^d, \ k=0,\ldots,M-1``, and coefficients ``f_j \in \mathbb{C}, j \in I_M^l``.

## Algorithm

For simplicity, we only consider the univariate case (``d=1``). The multivariate case, which is also implemented in the module, can be derived analogously from the NFFT algorithm (see e.g. [[Schmischke, 2018](#Schmischke2018)] or [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)], Chapter 7). Hence, we want to approximate 

```math
    f^c \colon [ 0, \pi ] \to \mathbb{C}, \quad f^c(x_j) = \sum_{k=0}^{N-1} \hat{f}^{c}_k \cos(k x_j) = \sum_{k=-N}^{N-1} \hat{f}_k \mathrm{e}^{ \mathrm{i} k x_j},
```
with ``x_j \in [0,\pi], \ j=0, \ldots, M-1`` the nodes and ``\hat{f}_k^c \in \mathbb{R}, \ k = 0, \ldots, N-1`` the coefficients. Observe that the second equality sign is only justified when ``\hat{f}_0 = \hat{f}_0^c,\ \hat{f}_{_N}=0,\ \hat{f}_k = \hat{f}_{-k}= \frac{1}{2}\hat{f}_k^c``. Similarly to the NFFT, we choose ``g_l\in \mathbb{R}`` in the following term such that ``s_1`` with a given window function ``\varphi\in L_2(\mathbb{R})\cap L_1(\mathbb{R})`` approximates ``f``

```math
    s_1(x) \coloneqq \sum_{\ell=0}^{\sigma N} g_\ell\ \tilde{\varphi}\left(x-\frac{\pi\ell}{\sigma N}\right), \quad x\in \mathbb{R}.
```

``\sigma \in \mathbb{N}`` is an oversampling factor fulfilling ``\sigma N \in \mathbb{N}``. Replacing ``\varphi`` by its Fourier series ``\varphi (x) = \sum_{k\in \mathbb{N}} \hat{g}_k \mathrm{e}^{\mathrm{i} k x}``, comparing the definition of ``s_1(\cdot)`` and ``f^c(\cdot)`` and taking the symmetry into account, we obtain 

```math
    g_{\ell} = \frac{1}{\sigma N} \sum_{k=0}^{\sigma N} (\varepsilon_{\sigma N}(k))^2 \hat{g}_k \cos ( \frac{2 \pi k \ell}{\sigma N} ), \quad \ell = 0, \ldots, \sigma N.
```

Here, ``\varepsilon_{\sigma N}(j)\coloneqq 1,\ j=1,\ldots,\sigma N-1``,  and  ``\varepsilon_{\sigma N}(0)\coloneqq \varepsilon_{\sigma N}(\sigma N)\coloneqq \frac{\sqrt{2}}{2}``. Now, we can compute ``g_l`` by an DCT-I of length ``\sigma N +1`` and obtain ``f(x) \approx s_1(x)``. 
\
Next, we want to replace ``\varphi`` by its truncation

```math
    \psi(x) = \varphi(x) \chi_Q (x) = \begin{cases} \varphi(x) & x\in Q\coloneqq\left[ -\frac{2\pi m}{\sigma N},\frac{2\pi m}{\sigma N} \right] , \\ 0, & else. \end{cases}
```

As in the definition of ``s_1(\cdot)``, we replace ``\psi`` by its periodization ``\tilde{\psi}`` and define

```math
  	s(x)\coloneqq \sum_{\ell=\lfloor 2\sigma N x\rfloor -m}^{\lceil 2\sigma N x\rceil +m} g_\ell\ \tilde{\psi}\left(x-\frac{\pi \ell}{\sigma N}\right),\quad x\in\mathbb{R}
```

Finally, we arrive at ``f(x)\approx s_1(x) \approx s(x)``. This gives rise to Algorithm \ref{alg:NFCT}.

# Pseudocode


The algorithm for the transposed problem 

```math
  	h(k) \coloneqq \sum_{j=0}^{M-1} h_j \ \cos(kx_j), \quad k=0,\ldots,N-1
```
can easily be derived from the duality of the two problems.

## Functions

```@docs
  	finalize_plan
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

```@raw html
<ul>
<li id="Schmischke2018">[<a>Schmischke, 2018</a>]
  M. Schmischke. Nonequispaced Fast Fourier Transform (NFFT) Interface for Julia.</emph>
  2018.
  arXiv: <a href="https://arxiv.org/abs/1810.09891">1512.02814</a>.
</li>
</ul>
```
