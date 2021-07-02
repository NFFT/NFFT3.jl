# [Fast Summation](@id fastsum)

```@meta
    CurrentModule = NFFT3
```

## Fastsum algorithm

The fast summation algorithm evaluates the function 

```math
    f(\pmb{y}) \coloneqq \sum_{k=1}^{N} \alpha_k \ \mathscr{K} (\pmb{y}-\pmb{x}_k)
``` 

for given arbitrary source knots ``\pmb{x}_k \in \mathbb{R}^d,\ k = 1, \ldots, N`` and a given kernel function ``\mathscr{K}(\pmb{x}) \coloneqq K(\lVert \pmb{x} \rVert_2),\ \pmb{x} \in \mathbb{R}^d``. Here, ``K`` is required to be infinitely differentiable at all points in ``\mathbb{R}^d \setminus \{ 0 \}``. If ``K`` is even infinitely differentiable at the origin, ``\mathscr{K}`` is called **nonsingular kernel function**, otherwise **singular kernel function**. 
\
The evaluation is done at ``M`` different points ``\pmb{y}_j \in \mathbb{R}^d,\ j=1,2, \ldots, M``. W.l.o.g. we assume ``\pmb{y}_j \neq \pmb{x}_k  \ (j\neq k)``, i.e., one wants to compute

```math
    f(\pmb{y}_j) \coloneqq \sum_{k=1}^{N} \alpha_k \ \mathscr{K}(\pmb{y}_j-\pmb{x}_k), \qquad j=1,\ldots,M.
```

We replace ``\mathscr{K}`` in the definition above by a periodic function. Hence, we assume ``\lVert \pmb{y}_k \rVert_2 \leq \frac{1}{2} \ (\frac{1}{2} - \varepsilon_B)`` and ``\lVert \pmb{x}_k \rVert_2 \leq \frac{1}{2} \ (\frac{1}{2} - \varepsilon_B)`` with ``\varepsilon_B \in (0,0.5)``. By the triangle inequality, we conclude 

```math
    \pmb{x}_j-\pmb{x}_k \in [-\frac{1}{2}+\varepsilon_B, \frac{1}{2}-\varepsilon_B].
```
Next, we regularize the kernel near zero and ``\pm 0.5`` to obtain the function

```math
    K_R(r) \coloneqq \begin{cases} T_I(r) & \lvert r \rvert < \varepsilon_I \\ K(r) & \varepsilon_I < \lvert r \rvert \leq 0.5 - \varepsilon_B \\ T_B(\lvert r \rvert) & 0.5 - \varepsilon_B < \lvert r \rvert \leq 0.5. \end{cases}
```

Here, ``T_I,\, T_B \in \mathcal{P}_{2p-1}`` are chosen such that the 1-periodic extension of ``K_R`` is of class ``\mathcal{C}^{p-1}``. Further details about the computation of ``T_I, T_B`` are given in Chapter 7.5 of [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)]. 
\
We define 

```math
    \mathscr{K}_R(\pmb{x}) \coloneqq \begin{cases} K_R(\lVert \pmb{x} \rVert_2) & \lVert \pmb{x} \rVert_2 < 0.5 \\ T_B(0.5) & \lVert \pmb{x} \rVert_2 \geq 0.5 \end{cases}
```

Since ``\mathscr{K}_R`` is sufficiently smooth, we replace it by its partial sum  ``\mathscr{K}_{RF}(\pmb{x}) \coloneqq \sum_{\pmb{\ell} \in I_n^d} b_{\pmb{\ell}}\, \mathrm{e}^{2\pi\mathrm{i} \, \pmb{\ell} \cdot \pmb{x}}`` and obtain ``\mathscr{K} \approx \mathscr{K} - \mathscr{K}_R + \mathscr{K}_{RF}``. Using quadrature rules, we replace ``b_{\pmb{\ell}}`` by

```math
    b_{\pmb{\ell}} \coloneqq \frac{1}{n^d} \sum_{\pmb{h} \in I_n^d} \mathscr{K}_R(\frac{\pmb{h}}{n}) \, \mathrm{e}^{-2 \pi \mathrm{i} \, \pmb{h} \cdot \pmb{\ell} / n}, \quad \pmb{\ell} \in I_n^d \coloneqq [-\frac{n}{2},\frac{n}{2}]^d \cap \mathbb{Z}^d 
```

Hence, we have to compute the **nearfield sum**

```math
    f_{NE}(\pmb{x}) \coloneqq \sum_{k=1}^{N} \alpha_k \, \mathscr{K}_{NE}(\pmb{x}-\pmb{x}_k), \quad \mathscr{K}_{NE} \coloneqq \mathscr{K} - \mathscr{K}_R
```

and the *far field sum*

```math
    f_{RF}(\pmb{x}) \coloneqq \sum_{k=1}^{N} \alpha_k \, \mathscr{K}_{RF} (\pmb{x}-\pmb{x}_k)
```

and approximate ``f`` by ``\tilde{f} \coloneqq f_{NE} + f_{RF}``.
\
From ``\pmb{y}_j-\pmb{x}_k \in [-0.5+\varepsilon_B, 0.5-\varepsilon_B]`` and the definition of ``K_R(\cdot)`` we conclude that for the evaluation of ``f_{NE}`` we only need to consider  ``\pmb{x}`` fulfilling ``\lVert \pmb{x} \rVert \leq \varepsilon_I``. Assuming a certain kind of uniform distribution of the ``\pmb{y}_j,\ \pmb{x}_k``, we obtain an arithmetic cost of ``\mathcal{O}(M)``.
\
For the evaluation of ``f_{RF}(\cdot)``, we rearrange the sums

```math
    f_{RF}(\pmb{y})= \sum_{k=1}^{N} \alpha_k \sum_{\pmb{\ell} \in I_n^d} b_{\pmb{\ell}} \ \mathrm{e}^{ 2\pi\mathrm{i} \, \pmb{\ell} \cdot ( \pmb{y} - \pmb{x}_k )} = \sum_{\pmb{\ell} \in I_n^d} b_{\pmb{\ell}}\, (\sum_{k=1}^{N} \alpha_k \ \mathrm{e}^{- 2\pi\mathrm{i} \,\pmb{\ell} \cdot \pmb{x}_k})\, \mathrm{e}^{2\pi\mathrm{i} \, \pmb{\ell} \cdot \pmb{y}}
```

where the inner sum can be computed by an adjoint NFFT which is then followed by a NFFT for the computation of the outer sum. 

### Pseudocode

**Input:** ``\alpha_k \in \mathbb{C} \text{ for } k = 1, \ldots, N, \pmb{x}_k \in \mathbb{R}^d \text{ for } k = 1, \ldots, M \text{ with } \lVert \pmb{x}_k \rVert_2 \leq \frac{1}{2} \ (0.5 - \varepsilon_B), \pmb{y}_j \in \mathbb{R}^d \text{ for } j = 1, \ldots, M \text{ with } \lVert \pmb{y}_j \rVert_2 \leq \frac{1}{2} \ (0.5 - \varepsilon_B).``

**Precomputation:** Compute the polynomials ``T_I`` and ``T_B``. Compute ``(b_{\pmb{\ell}})_{\pmb{\ell} \in I^d_n}``. Compute ``\mathscr{K}_{NE}(\pmb{y}_j - \pmb{x}_k) \text{ for } j = 1, \ldots, M \text{ and } k \in I^{NE}_{\varepsilon_I}(j), \text{ where } I^{NE}_{\varepsilon_I}(j) \coloneqq \{ k \in \{ 1, \ldots, N \} \colon \lVert \pmb{y}_j - \pmb{x}_k \rVert_2 < \varepsilon_I \}``.

1. For each ``\pmb{\ell} \in I^d_n`` compute ``a_{\pmb{\ell}} \coloneqq \sum^{N}_{k = 1} \alpha_k \, e^{ -2\pi \mathrm{i} \, \pmb{\ell} \cdot \pmb{x}_k}`` using the d-variate adjoint NFFT of size ``n \times \ldots \times n``, see Algorithm 7.3 in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)].
2. For each ``\pmb{\ell} \in I^d_n`` compute the products ``d_{\pmb{\ell}} \coloneqq a_{\pmb{\ell}} b_{\pmb{\ell}}``.
3. For ``j = 1, \ldots, M`` compute the far field sums ``f_{RF}(\pmb{y}_j) = \sum_{\pmb{\ell} \in I_n^d} d_{\pmb{\ell}} \, \mathrm{e}^{ 2\pi\mathrm{i} \, \pmb{\ell} \cdot \pmb{y}_j}`` using the d-variate NFFT of size ``n \times \ldots \times n``, see Algorithm 7.1 in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)].
4. For ``j = 1, \ldots, M`` compute the near field sums ``f_{NE}(\pmb{y}_j) \coloneqq \sum_{k \in I^{NE}_{\varepsilon_I}(j)} \alpha_k \, \mathscr{K}_{NE}(\pmb{y}_j - \pmb{x}_k)``.
5. For ``j = 1, \ldots, M`` compute the near field corrections ``\tilde{f}(\pmb{y}_j) \coloneqq f_{NE}(\pmb{y}_j) + f_{RF}(\pmb{y}_j)``.

**Output:** ``\tilde{f}(\pmb{y}_j), \, j = 1, \ldots, M``, approximate values of ``f(\pmb{y}_j)``.

**Computational cost:** ``\mathcal{O}(m^d (M_1+M_2)+ n^d \log n) `` 

## Plan structure

```@docs
    FASTSUM
```

## Functions

```@docs
  	fastsum_init
    fastsum_finalize_plan
    fastsum_trafo
    fastsum_trafo_exact
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
