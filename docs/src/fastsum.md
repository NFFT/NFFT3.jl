# [Fast Summation (fastsum)](@id fastsum)

```@meta
    CurrentModule = NFFT3
```

## fastsum algorithm

The fast summation algorithm evaluates the function 

```math
    f(y) \coloneqq \sum_{k=1}^{M_1} \alpha_k \ \mathscr{K} (y-x_k)
``` 

for given (nonequispaced) source knots ``x_k \in \mathbb{R}^d,\ k = 1, \ldots, M_1`` and a given kernel function ``\mathscr{K} (x) \coloneqq K(\lVert x \rVert_2),\ x \in \mathbb{R}^d``. Here, ``K`` is required to be infinitely differentiable at ``x \in \mathbb{R} \setminus \{ 0 \}``. If ``K`` is even infinitely differentiable at ``0``, ``\mathscr{K}`` is called *nonsingular kernel function*, otherwise *singular kernel function*. 
\
The evaluation is done at ``M_2`` different points ``y_j \in \mathbb{R}^d,\ j=0, \ldots, M_2``. W.l.o.g. we assume ``y_j \neq x_i  \ (j\neq i)``, i.e., one wants to compute

```math
    f(y_j) \coloneqq \sum_{k=1}^{M_1} \alpha_k \ \mathscr{K} (y_j-x_k), \qquad j=1,\ldots,M_2.
```

We replace ``\mathscr{K}`` in the definition above by a periodic function. Hence, we assume ``\lVert y_k \rVert_2 \leq \frac{1}{2} \ (\pi - \varepsilon_B)`` and ``\lVert x_k \rVert_2 \leq \frac{1}{2} \ (\pi - \varepsilon_B)`` with ``\varepsilon_B \in (0,\pi)``. By the triangle inequality, we conclude 

```math
    y_j-x_i \in [-\pi+\varepsilon_B, \pi-\varepsilon_B]
```
Next, we regularize the kernel near zero and ``\pm \pi``:

```math
    K_R (x) \coloneqq \begin{cases} T_I(x) & \lvert x \rvert < \varepsilon_I \\ K(x) & \varepsilon_I < \lvert x \rvert \leq \pi - \varepsilon_B \\ T_B(\lvert x \rvert) & \pi - \varepsilon_B < \lvert x \rvert \leq \pi \end{cases}
```

With ``\lvert x \rvert \coloneqq (\lvert x_k \rvert)_{k=0}^{N-1}`` we mean to componentwise absolute value and ``T_I,\, T_B \in \mathcal{P}_{2p-1}`` are chosen such that the ``2\pi``-periodic extension of ``K_R`` is of class ``\mathcal{C}^{p-1}``. Further details about the computation of ``T_I, T_B`` are given in Chapter 7.5, [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)]. 
\
We define 

```math
    \mathscr{K}_R (x) \coloneqq \begin{cases} K_R (\lVert x \rVert_2) & \lVert x \rVert < \pi \\ T_B (\pi) & \lVert x \rVert_2 \geq \pi \end{cases}
```

Since ``\mathscr{K}_R`` is sufficiently smooth, we replace it by its partial sum  ``\mathscr{K}_{RF} (x) \coloneqq \sum_{\pmb{l} \in I_n^d} b_{\pmb{l}}\mathrm{e}^{\mathrm{i} \ \pmb{l} \cdot x}`` and obtain ``\mathscr{K} \approx \mathscr{K} - \mathscr{K}_R + \mathscr{K}_{RF}``. Using quadrature rules, we replace ``b_{\pmb{l}}`` by

```math
    b_{\pmb{l}} \coloneqq \frac{1}{n^d} \sum_{j \in I_n^d} \mathscr{K}_R(\frac{2 \pi j}{n}) \, \mathrm{e}^{-2 \pi j \cdot \pmb{l} / n}, \quad \pmb{l} \in I_n^d \coloneqq [- \frac{n}{2},\frac{n}{2} ]^d \cap \mathbb{Z}^d 
```

Hence, we have to compute the *nearfield sum*

```math
    f_{NE} (x) \coloneqq \sum_{k=1}^{M_1} \alpha_k \, \mathscr{K}_{NE} (x-x_k), \quad \mathscr{K}_{NE} \coloneqq \mathscr{K} - \mathscr{K}_R
```

and the *far field sum*

```math
    f_{RF} (x) \coloneqq \sum_{k=1}^{M_1} \alpha_k \, \mathscr{K}_{RF} (x-x_k)
```

and approximate ``f`` by ``\tilde{f} \coloneqq f_{NE} + f_{RF}``.
\
From ``y_j-x_i \in [-\pi+\varepsilon_B, \pi-\varepsilon_B]`` and the definition of ``K_R (\cdot)`` we conclude that for the evaluation of ``f_{NE}`` we only need to consider  ``x`` fulfilling ``\lVert x \rVert \leq \varepsilon_I``. Assuming a certain kind of uniform distribution of the ``x_k,\ y_i``, we obtain an arithmetic cost of ``\mathcal{O}(M_2)``.
\
For the evaluation of ``f_{RF}(\cdot)``, we rearrange the sums

```math
    f_{RF}(y)= \sum_{k=1}^{M_1} \alpha_k \sum_{\pmb{l} \in I_n^d} b_{\pmb{l}} \ \mathrm{e}^{ \mathrm{i} \ \pmb{l} \cdot ( y - x_k )} = \sum_{\pmb{l} \in I_n^d} b_{\pmb{l}}\, (\sum_{k=1}^{M_1} \alpha_k \ \mathrm{e}^{- \mathrm{i} \pmb{l} x_k})\, \mathrm{e}^{\mathrm{i} \ \pmb{l} \cdot y}
```

where the inner sum can be computed by a NFFT``^\mathrm{T}`` which is then followed by a NFFT for the computation of the outer sum. 

### Pseudocode

**Input:** ``\alpha_k \in \mathbb{C} \text{ for } k = 1, \ldots, M_1, \; x_k \in \mathbb{R}^d \text{ for } k = 1, \ldots, M_1 \text{ with } \lVert x_k \rVert_2 \leq \frac{1}{2} \ (\pi - \varepsilon_B), \; y_j \in \mathbb{R}^d \text{ for } j = 1, \ldots, M_2 \text{ with } \lVert y_j \rVert_2 \leq \frac{1}{2} \ (\pi - \varepsilon_B).``

**Precomputation:** Compute the polynomials ``T_I`` and ``T_B``. Compute ``(b_{\pmb{l}})_{\pmb{l} \in I^d_n}``. Compute ``\mathscr{K}_{NE} (y_j - x_k) \text{ for } j = 1, \ldots, M_2 \text{ and } k \in I^{NE}_{\varepsilon_I}(j), \text{ where } I^{NE}_{\varepsilon_I}(j) \coloneqq \{ k \in \{ 1, \ldots, M_1 \} \colon \lVert y_j - x_k \rVert_2 < \varepsilon_I \}``.

1. For each ``\pmb{l} \in I^d_n`` compute ``\alpha_{\pmb{l}} \coloneqq \sum^{M_1}_{k = 1} \alpha_k e^{ - \mathrm{i} \ \pmb{l} \cdot x_k}`` using the d-variate NFFT``^\mathrm{T}`` of size ``n \times \ldots \times n``, see Algorithm 7.3 in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)].
2. For each ``\pmb{l} \in I^d_n`` compute the products ``d_{\pmb{l}} \coloneqq a_{\pmb{l}} b_{\pmb{l}}``.
3. For ``j = 1, \ldots, M_2`` compute the far field sums ``f_{RF}(y_j) = \sum_{\pmb{l} \in I_n^d} d_{\pmb{l}} \ \mathrm{e}^{ \mathrm{i} \ \pmb{l} \cdot y_j}`` using the d-variate NFFT of size ``n \times \ldots \times n``, see Algorithm 7.1 in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)].
4. For ``j = 1, \ldots, M_2`` compute the near field sums ``f_{NE} (y_j) \coloneqq \sum_{k \in I^{NE}_{\varepsilon_I}(j)} \alpha_k \mathscr{K}_{NE} (y_j - x_k)``.
5. For ``j = 1, \ldots, M_2`` compute the near field corrections ``\tilde{f}(y_j) \coloneqq f_{NE}(y_j) + f_{RF}(y_j)``.

**Output:** ``\tilde{f}(y_j), \, j = 1, \ldots, M_2``, approximate values of ``f(y_j)``.

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
