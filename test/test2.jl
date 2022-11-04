using NFFT3

x = [0.1,0.8]
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

#F1.x[1] -= 0.3
#F2.x[1] -= 0.3
#F3.x[1] -= 0.3

F1.fhat = fhat
F2.fhat = fhat
F3.fhat = fhat

F1.fhat[4] += 10.0
#F2.fhat[4] += 10.0
F3.fhat[4] += 10.0

nfft_trafo(F1)
nfft_trafo(F2)
nfft_trafo_direct(F3)

f1 = F1.f
f2 = F2.f
f3 = F3.f

println(f1)
println(f2)
println(f3)

#println(F1)