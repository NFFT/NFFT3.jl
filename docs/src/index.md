# Welcome to NFFT3.jl

Fast Fourier transforms (FFTs) belong to the '10 algorithms with the greatest influence on the development and practice of science and engineering in the 20th century'. The classic algorithm computes the discrete Fourier transform

```math
    f_j = \sum^{\frac{N}{2} - 1}_{k = -\frac{N}{2}} \hat{f}_k \mathrm{e}^{2 \pi \ \mathrm{i} \ \frac{k \ j}{N}}
```

for ``j = - \frac{N}{2}, \ldots, \frac{N}{2} - 1 `` and given complex coefficients ``\hat{f}_{k} \in \mathbb{C}``. Using a divide and conquer approach, the number of floating point operations is reduced from ``\mathcal{O}(N^2)`` for a straightforward computation to only ``\mathcal{O}(N \log N)``. In conjunction with publicly available efficient implementations the fast Fourier transform has become of great importance in scientific computing.
\
However, two shortcomings of traditional schemes are the need for equispaced sampling and the restriction to the system of complex exponential functions. The NFFT (nonequispaced fast Fourier transform or nonuniform fast Fourier transform, NUFFT) is a C subroutine library for computing the nonequispaced discrete Fourier transform (NDFT) and its generalisations in one or more dimensions, of arbitrary input size, and of complex data.


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