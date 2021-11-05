using NFFT3
using Test
using LinearAlgebra
using Aqua

Aqua.test_all(NFFT3)

include( "NFFT.jl" )
include( "NFCT.jl" )
include( "NFST.jl" )
include( "fastsum.jl" )
