# [Nonequispaced Fast Cosine Transform (NFMT)](@id NFMT_site)

```@meta
    CurrentModule = NFFT3
```

## NDMT and NFMT

We consider for $\pmb{d}\in\{\exp,\cos,\mathrm{alg}\}^d$ the trigonometric polynomial

$$f^{\pmb{d}}(\pmb{x}) \coloneqq \sum_{\pmb{k} \in I_{\pmb{N},\pmb{d}}^d} \hat{f}_{\pmb{k}}^{\pmb{d}} \, \phi_{\pmb{k}}^{\pmb{d}}(\pmb{x}), \quad \pmb{x} \in \mathbb{R}^d,$$

with 

$$\phi_{\pmb{k}}^{\pmb{d}}(\pmb{x})=\prod_{j=1}^d\begin{cases}1,&k_j=0\\\exp(2\pi\mathrm{i}k_jx_j),&d_j=\exp,\;k_j\neq0\\
\sqrt{2}\cos(\pi k_jx_j),&d_j=\cos,\;k_j\neq0\\
\sqrt{2}\cos(k_j\arccos(2x_j-1)),&d_j=\mathrm{alg},\;k_j\neq0\end{cases} $$

and multibandlimit $\pmb{N} \in (2\mathbb{N})^d$ and index set

$$I_{\pmb{N},\pmb{d}}^d \coloneqq \overset{d}{\underset{j=1}{\vphantom{\mathop{\raisebox{-.5ex}{\hbox{\huge{$\times$}}}}}⨉}}\begin{cases}\Big\{-\frac{N_j}{2},-\frac{N_j}{2}+1,\ldots,\frac{N_j}{2}\Big\},&d_j=\exp\\\Big\{0,1,\ldots,\frac{N_j}{2}\Big\},&d_j\neq\exp\end{cases}.$$

The NDMT is the evaluation of

$$f^{\pmb{d}}(\pmb{x}_j) \coloneqq \sum_{\pmb{k} \in I_{\pmb{N},\pmb{d}}^d} \hat{f}_{\pmb{k}}^{\pmb{d}} \, \phi_{\pmb{k}}^{\pmb{d}}(\pmb{x}_j)$$

at arbitrary nodes $\pmb{x}_j \in [0,1]^d$ for given coefficients $\hat{f}_{\pmb{k}}^{\pmb{d}} \in \mathbb{R}, \pmb{k} \in I_{\pmb{N},\pmb{d}}^d$. Similarly to the NDFT, the transposed NDMT is the evaluation of

$$\hat{h}^{\pmb{d}}_{\pmb{k}} = \sum_{j=1}^M f^{\pmb{d}}_j \, \phi_{\pmb{k}}^{\pmb{d}}(\pmb{x}_j)$$

for the frequencies $\pmb{k} \in I_{\pmb{N},\pmb{d}}^d$ with given coefficients $f^{\pmb{d}}_j \in \mathbb{R}, j = 1,2,\ldots,M$.

We modify the [NFFT](@ref NFFT_site) in order to derive a fast algorithm for the computation of the NDMT and transposed NDMT, obtaining the NFMT and its transposed counterpart. For details we refer to [[Potts, Schröter, 2024](#PottsSchröter2024)].

## Plan structure

```@docs
    NFMT{D}
```

## Functions

```@docs
    nfmt_trafo
    nfmt_trafo_direct
    nfmt_transposed
    nfmt_transposed_direct
  	nfmt_finalize_plan
    nfmt_init
```

## Literature

```@raw html
<ul>
<li id="PottsSchröetr2024">[<a>Potts, Schröter, 2024</a>]
  D. Potts, P.Schröter. Linear Algebra Appl.</emph>
  arXiv: <a href="https://arxiv.org/abs/2306.09174">2306.09174</a>.
</li>
</ul>
```