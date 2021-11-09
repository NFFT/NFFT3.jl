# NFFT3.jl

Julia Interface for the [NFFT C library](https://github.com/NFFT/nfft) 

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://nfft.github.io/NFFT3.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://nfft.github.io/NFFT3.jl/dev)
[![ci](https://github.com/NFFT/NFFT3.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/NFFT/NFFT3.jl/actions?query=workflow%3ACI+branch%3Amain)
[![codecov](https://codecov.io/gh/NFFT/NFFT3.jl/branch/main/graph/badge.svg?token=YCTMXP64FK)](https://codecov.io/gh/NFFT/NFFT3.jl)
[![Aqua QA](https://img.shields.io/badge/Aqua.jl-%F0%9F%8C%A2-aqua.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5656757.svg)](https://doi.org/10.5281/zenodo.5656757)

`NFFT3.jl` provides the following fast algorithms:
- nonequispaced fast Fourier transform (NFFT) 
- nonequispaced fast cosine transform (NFCT) 
- nonequispaced fast sine transform (NFST)
- fast summation (fastsum) 

## Getting started

In Julia you can get started by just typing

```julia
] add NFFT3
```

then checkout the [documentation](https://nfft.github.io/NFFT3.jl/stable/).
