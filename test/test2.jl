using NFFT3

x = [0.05,0.4]
fhat = [       0.0 + 0.0im
0.7071067811865475 + 0.0im
0.7071067811865475 + 0.0im
0.7071067811865475 + 0.0im
0.7071067811865475 + 0.0im
0.7071067811865475 + 0.0im
1.0 + 0.0im
0.7071067811865475 + 0.0im
0.7071067811865475 + 0.0im
0.7071067811865475 + 0.0im
0.7071067811865475 + 0.0im
0.7071067811865475 + 0.0im]

F1 = NFFT{1}(NTuple{1,Int32}(12), Int32(2), NTuple{1,Int32}(64), Int32(8), NFFT3.f1_default_1d, NFFT3.f2_default)
F2 = NFFT{1}(NTuple{1,Int32}(12), Int32(2), NTuple{1,Int32}(64), Int32(12), NFFT3.f1_default_1d, NFFT3.f2_default)
F3 = NFFT{1}(NTuple{1,Int32}(12), Int32(2), NTuple{1,Int32}(64), Int32(8), NFFT3.f1_default_1d, NFFT3.f2_default)

F1.x = x
F2.x = x
F3.x = x

F1.fhat = fhat
F2.fhat = fhat
F3.fhat = fhat

println(F1)
#println(F1.x)
println(unsafe_wrap(Vector{Float64}, Core.getfield(F1, :x), F1.M))
#println(F1.fhat)
println(unsafe_wrap(Vector{ComplexF64}, Core.getfield(F1, :fhat), prod(F1.N)))
println(F1.N)
println(F1.M)
println(F1.n)
println(F1.m)
println(F1.f1)
println(F1.f2)
println(F1.init_done)
println(F1.finalized)
nfft_trafo(F1)
println(F1.f)


#nfft_trafo(F2)
#nfft_trafo_direct(F3)

#f1 = F1.f
#f2 = F2.f
#f3 = F3.f

#e = [1.9733023985986597,0.7283598003728021]


#println(f2)
#println(f3)

#println(F1)