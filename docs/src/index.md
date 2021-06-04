# Welcome to NFFT3.jl

## NNFT online 
Fast Fourier transforms (FFTs) belong to the '10 algorithms with the greatest influence on the development and practice of science and engineering in the 20th century'. The classic algorithm computes the discrete Fourier transform

```math
	f_j = \sum^{\frac{N}{2} - 1}_{k = -\frac{N}{2}} \hat{f}_k \mathrm{e}^{2 \pi \ \mathrm{i} \ \frac{k \ j}{N}}
```

for ``j = - \frac{N}{2}, \ldots, \frac{N}{2} - 1 `` and given complex coefficients ``\hat{f}_{k} \in \mathbb{C}``. Using a divide and conquer approach, the number of floating point operations is reduced from ``\mathcal{O}(N^2)`` for a straightforward computation to only ``\mathcal{O}(N \log N)``. In conjunction with publicly available efficient implementations the fast Fourier transform has become of great importance in scientific computing.
\
However, two shortcomings of traditional schemes are the need for equispaced sampling and the restriction to the system of complex exponential functions. The NFFT (nonequispaced fast Fourier transform or nonuniform fast Fourier transform, NUFFT) is a C subroutine library for computing the nonequispaced discrete Fourier transform (NDFT) and its generalisations in one or more dimensions, of arbitrary input size, and of complex data.
\
More precisely, we collect the possible frequencies ``\mathbf{k} \in \mathbb{Z}^d`` in the multi-index set

```math
  	I_{\mathbf{N}} := \left\{ \mathbf{k}=\left(k_t\right)_{t=0,\ldots,d-1} \in \mathbb{Z}^d: - \frac{N_t}{2} \le k_t \lt \frac{N_t}{2} ,\;t=0,\ldots,d-1\right\},
```

where ``\mathbf{N}=\left(N_t \right)_{t=0,\ldots,d-1}`` is the multibandlimit, i.e., ``N_t \in 2\mathbb{N}``. For a finite number of given Fourier coefficients ``\hat{f}_{\mathbf{k}} \in \mathbb{C}, \; \mathbf{k} \in I_{\mathbf{N}}``, we consider the fast evaluation of the trigonometric polynomial

```math
  	f\left(\mathbf{x}\right) := \sum_{ \mathbf{k}\in I_{ N}} \hat{f}_{\mathbf{ k}} {\mathrm{e}}^{-2\pi{\mathrm{i}}\mathbf{k}\mathbf{ x}}
```
 
at given nonequispaced nodes ``\mathbf{x}_j \in \mathbb{T}^d, \; j=0,\ldots, M-1``, from the ``d``-dimensional torus as well as the adjoint problem, the fast evaluation of sums of the form

```math
  	\hat h_{\mathbf{k}} := \sum_{j=0}^{M-1} {f}_{j} {\mathrm{e}}^{2\pi{\mathrm{i}}\mathbf{k}\mathbf{ x}_j}.
```

In general, the adjoint NDFT is not the inverse transform of the NDFT.

The NFFT is a C subroutine library for computing the nonequispaced discrete Fourier transform (NDFT) in one or more dimensions, of arbitrary input size, and of complex data. New: A Matlab interface is part of the NFFT3. We believe that our library, which is free software, and based on FFTW (FFTW 3.x) should become the NFFT library of choice for most applications.

The NFFT package has been developed at the Mathematical Institute of the University of L\"ubeck, at the Mathematical Institute of the University Osnabr\"uck and at the Faculty of Mathematics of the Chemnitz University of Technology by Jens Keiner, Stefan Kunis and Daniel Potts. Further contributions, in particular applications, are due to Dr. Markus Fenn, Steffen Klatt, Tobias Knopp and Antje Vollrath. The support for OpenMP was developed by Toni Volkmer. Many contributions to the release 3.3.* and later are done by Toni Volkmer and Michael Quellmalz.

## Michael

The nonequispaced fast Fourier transform [[Keiner, Kunis, Potts, 2006](#KeinerKunisPotts2006)] (NFFT or NUFFT) overcomes one of the main shortcomings of the FFT - the need for an equispaced sampling grid. Considering a ``d``-dimensional trigonometric polynomial 

```math
  	f(\pmb{x}) \coloneqq \sum_{ \pmb{k} \in I_{\pmb{N}}} \hat{f}_{\pmb{k}} \mathrm{e}^{-2\pi\mathrm{i} \pmb{k}\pmb{x}}
```

with an index set ``I_{\pmb{N}} \coloneqq \{ \pmb{k} \in \mathbb{Z}^d: -\frac{N_i}{2} \leq \pmb{k}_i \leq \frac{N_i}{2}-1, i=0,\ldots,d-1 \}`` where ``\pmb{N} \in 2\mathbb{N}^d`` is the multibandlimit, the NDFT is its evaluation at ``M \in \mathbb{N}`` nonequispaced points ``\pmb{x}_j \in \mathbb{T}^d`` for ``j = 0, 1, \ldots, M``,

```math
  	f(\pmb{x}_j) =\sum_{\pmb{k} \in I_{\pmb{N}}} \hat{f}_{\pmb{k}} \mathrm{e}^{-2 \pi \mathrm{i} \pmb{k} \cdot \pmb{x}_j},
```

with given coefficients ``\hat{f}_{\pmb{k}} \in \mathbb{C}`` where we identify the smooth manifold of the torus ``\mathbb{T}`` with ``[-1/2, 1/2)``. The NFFT is an algorithm for the fast evaluation of the sums ``f(\pmb{x}_j)`` and the adjoint problem, the fast evaluation of

```math
	\hat{h}_{\pmb{k}} = \sum_{j = 0}^{M-1} f_j \mathrm{e}^{2 \pi \mathrm{i} \pmb{k} \pmb{x}_j}, \pmb{k} \in I_{\pmb{N}}
```

for given coefficients ``f_j \in \mathbb{C}``. The available NFFT3 library [[Keiner, Kunis, Potts, NFFT3](#KeinerKunisPottsNFFT3)] provides C routines for the NFFT, applications such as the fast evaluation of sums

```math
  	g(\pmb{y}_j) \coloneqq \sum_{k=1}^{N} \alpha_k K(\lVert \pmb{y}_j - \pmb{x}_k \rVert_2), j = 1, \ldots, M,
```

for given coefficients ``\alpha_k \in \mathbb{C}``, nodes ``\pmb{x}_k,\pmb{y}_j \in \R^d``  and a radial kernel function ``K: [0,\infty) \to [0,\infty)``, and generalizations such as the NNFFT for nonequispaced nodes in time and frequency domain. The major NFFT3 release included interfaces for the numerical computing environments MATLAB and OCTAVE, expanding the user base.  

In the past years, a new dynamic programming language called Julia, see \cite{Bezanson2017}, caught interest in the field of numerical computing. The language is designed around the idea of high performance and flexibility with a lot of available packages for different purposes. The growing community around Julia has led to the idea of further expanding the NFFT3's user base by adding a corresponding interface to the library. 


### Literature

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