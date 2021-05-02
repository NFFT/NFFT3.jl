# [Fast Summation Algorithm](@id fastsum)

```@meta
    CurrentModule = NFFT3
```

```@docs
    FASTSUM
```

## Background

The fast summation algorithm evaluates the function 

```math
    f(y) \coloneqq \sum_{k=1}^{M_1} \alpha_k \ \mathcal{K} (y-x_k)
``` 

for given (nonequispaced) source knots ``x_k \in \mathbb{R}^d,\ k = 1, \ldots, M_1`` and a given kernel function ``\mathcal{K} (x) \coloneqq K(\lVert x \rVert_2),\ x \in \mathbb{R}^d``. Here, ``K`` is required to be infinitely differentiable at ``x\in \mathbb{R} \setminus \{ 0 \}``. If ``K`` is even infinitely differentiable at 0, ``\mathcal{K}`` is called \emph{nonsingular kernel function}, otherwise \emph{singular kernel function}. 

# Pseudocode

The evaluation is done at ``M_2`` different points ``y_j \in \mathbb{R}^d,\ j=0,\ldots,M_2``. W.l.o.g. we assume ``y_j \neq x_i  \ (j\neq i)``, i.e., one wants to compute

```math
    f(y_j) \coloneqq \sum_{k=1}^{M_1} \alpha_k \ \mathcal{K} (y_j-x_k), \qquad j=1,\ldots,M_2
```

We replace ``\mathcal{K}`` in the definition above by a periodic function. Hence, we assume ``\lVert y_k \rVert_2 \leq \frac{\pi}{2} - \frac{\varepsilon_B}{2}`` and ``\lVert x_k \rVert_2 \leq \frac{\pi}{2} - \frac{\varepsilon_B}{2}`` with ``\varepsilon_B \in (0,\pi)``. By the triangle inequality, we conclude 

```math
    y_j-x_i \in [-\pi+\varepsilon_B, \pi-\varepsilon_B]
```
Next, we regularize the kernel near zero and ``\pm \pi``:

```math
    K_R (x) \coloneqq \begin{cases} T_I(x) & \lvert x \rvert < \varepsilon_I \\ K(x) & \varepsilon_I < \lvert x \rvert \leq \pi - \varepsilon_B \\ T_B(\lvert x \rvert) & \pi - \varepsilon_B < \lvert x \rvert \leq \pi \end{cases}
```

With ``\lvert x \rvert \coloneqq (\lvert x_k \rvert)_{k=0}^{N-1}`` we mean to componentwise absolute value. ``T_I,\, T_B \in \mathscr{P}_{2p-1}`` are chosen such that the ``2\pi``-periodic extension of ``K_R`` is of class ``\mathcal{C}^{p-1}``. Further details about the computation of ``T_I,T_B`` is given in [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)]. 
\
We define 

```math
    \mathcal{K}_R (x) \coloneqq \begin{cases} K_R (\lVert x \rVert_2) & \lVert x \rVert < \pi \\ T_B (\pi) & \lVert x \rVert_2 \geq \pi \end{cases}
```

Since ``\mathcal{K}_R`` is sufficiently smooth, we replace it by its partial sum  ``\mathcal{K}_{RF} (x) \coloneqq \sum_{\ell \in I_n^d} b_{\ell}\mathrm{e}^{\mathrm{i} \ell\cdot x}`` and obtain `` \mathcal{K} \approx \mathcal{K}-\mathcal{K}_R + \mathcal{K}_{RF}``.  Using quadrature rules, we replace ``b_{\ell}`` by

```math
    b_{\ell} \coloneqq \frac{1}{n^d} \sum_{j \in I_n^d} \mathcal{K}_R(\frac{2 \pi j}{n}) \mathrm{e}^{-2 \pi j * \ell / n}, \quad \ell \in I_n^d \coloneqq [- \frac{n}{2},\frac{n}{2} ]^d \cap \mathbb{Z}^d 
```

Hence, we have to compute the \emph{nearfield sum}

```math
    f_{N} (x) \coloneqq \sum_{k=1}^{M_1} \alpha_k \mathcal{K}_{N} (x-x_k), \quad \mathcal{K}_N \coloneqq \mathcal{K} - \mathcal{K}_R
```

and the \emph{far field sum}

```math
  f_{RF} (x) \coloneqq \sum_{k=1}^{M_1} \alpha_k \mathcal{K}_{RF} (x-x_k)
```

and approximate ``f`` by ``\tilde{f} \coloneqq f_N + f_{RF}``.
\
From ``y_j-x_i \in [-\pi+\varepsilon_B, \pi-\varepsilon_B]`` and the definition of ``K_R (\cdot)`` we conclude that for the evaluation of ``f_N`` we only need to consider  ``x`` fulfilling ``\lVert x \rVert \leq \varepsilon_I``. Assuming a certain kind of uniform distribution of the ``x_k,\ y_i``, we obtain an arithmetic cost of ``\mathcal{O}(M_2)``.
\
For the evaluation of ``f_{RF}(\cdot)``, we rearrange the sums

```math
    f_{RF}(y)= \sum_{k=1}^{M_1} \alpha_k \sum_{\ell \in I_n^d} b_{\ell} \ \mathrm{e}^{ \mathrm{i} \ell ( y - x_k )} = \sum_{\ell \in I_n^d} b_{\ell}\, (\sum_{k=1}^{M_1} \alpha_k \ \mathrm{e}^{- \mathrm{i} \ell x_k})\, \mathrm{e}^{\mathrm{i} \ell \cdot y}
```
where the inner sum can be computed by a NFFT``^T`` which is then followed by a NFFT for the computation of the outer sum. 
\
In total, we obtain an arithmetic cost of ``\mathcal{O}(m^d(M_1+M_2)+(\rho n)^d \log(\rho n)) `` and arrive at Algorithm.

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
