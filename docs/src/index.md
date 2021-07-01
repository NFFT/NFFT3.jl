# Welcome to NFFT3.jl

The nonequispaced fast Fourier transform or NFFT, see [[Keiner, Kunis, Potts, 2006](#KeinerKunisPotts2006)] and [[Plonka, Potts, Steidl, Tasche, 2018](#PlonkaPottsSteidlTasche2018)], overcomes one of the main shortcomings of the FFT - the need for an equispaced sampling grid. Considering a ``d``-dimensional trigonometric polynomial 

```math
  	f(\pmb{x}) \coloneqq \sum_{ \pmb{k} \in I_{\pmb{N}}} \hat{f}_{\pmb{k}} \mathrm{e}^{-2\pi\mathrm{i} \pmb{k}\pmb{x}}
```

with an index set ``I_{\pmb{N}} \coloneqq \{ \pmb{k} \in \mathbb{Z}^d: -\frac{N_i}{2} \leq \pmb{k}_i \leq \frac{N_i}{2}-1, i=0,\ldots,d-1 \}`` where ``\pmb{N} \in 2\mathbb{N}^d`` is the multibandlimit, the nonequispaced fast Fourier transform (NDFT) is its evaluation at ``M \in \mathbb{N}`` nonequispaced points ``\pmb{x}_j \in \mathbb{T}^d`` for ``j = 0, 1, \ldots, M``,

```math
  	f(\pmb{x}_j) =\sum_{\pmb{k} \in I_{\pmb{N}}} \hat{f}_{\pmb{k}} \mathrm{e}^{-2 \pi \mathrm{i} \pmb{k} \cdot \pmb{x}_j},
```

with given coefficients ``\hat{f}_{\pmb{k}} \in \mathbb{C}`` where we identify the smooth manifold of the torus ``\mathbb{T}`` with ``[-1/2, 1/2)``. The NFFT is an algorithm for the fast evaluation of the sums ``f(\pmb{x}_j)`` as well as the adjoint problem, the fast evaluation of

```math
	\hat{h}_{\pmb{k}} = \sum_{j = 0}^{M-1} f_j \mathrm{e}^{2 \pi \mathrm{i} \pmb{k} \pmb{x}_j}, \pmb{k} \in I_{\pmb{N}}
```

for given coefficients ``f_j \in \mathbb{C}``. The available NFFT3 library [[Keiner, Kunis, Potts, NFFT3](#KeinerKunisPottsNFFT3)] provides C routines for the NFFT, applications such as the fast evaluation of sums

```math
  	g(\pmb{y}_j) \coloneqq \sum_{k=1}^{N} \alpha_k K(\lVert \pmb{y}_j - \pmb{x}_k \rVert_2), j = 1, \ldots, M,
```

for given coefficients ``\alpha_k \in \mathbb{C}``, nodes ``\pmb{x}_k,\pmb{y}_j \in \R^d``  and a radial kernel function ``K: [0,\infty) \to [0,\infty)``, and generalizations such as the NNFFT for nonequispaced nodes in time and frequency domain. 

The NFFT3 C library has been developed at the Mathematical Institute of the University of Luebeck, at the Mathematical Institute of the University Osnabrueck and at the Faculty of Mathematics of the Chemnitz University of Technology by Jens Keiner, Stefan Kunis and Daniel Potts. Further contributions, in particular applications, are due to Dr. Markus Fenn, Steffen Klatt, Tobias Knopp and Antje Vollrath. The support for OpenMP was developed by Toni Volkmer. Many contributions to the release 3.3.* and later have been done by Toni Volkmer, Michael Quellmalz, and Michael Schmischke.

This package offers a Julia wrapper for the NFFT, NFCT, NFST, and fastsum algorithms, see [[Schmischke, 2018](#Schmischke2018)]

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

```@raw html
<ul>
<li id="KeinerKunisPotts2006">[<a>Keiner, Kunis, Potts, 2006</a>]
  J. Keiner, S. Kunis, and D. Potts. Fast summation of radial functions on the sphere. </emph>
  Computing, 78:1--15, 2006.
  doi: <a href="https://doi.org/10.1007/s00607-006-0169-z">1512.02814</a>.
</li>
</ul>
```

```@raw html
<ul>
<li id="KeinerKunisPottsNFFT3">[<a>Keiner, Kunis, Potts, NFFT3</a>]
  J. Keiner, S. Kunis, and D. Potts. NFFT 3.0, C subroutine library. </emph>
  url: <a href="http://www.tu-chemnitz.de/~potts/nfft">1512.02814</a>.
</li>
</ul>
```