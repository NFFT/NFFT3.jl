using CpuId

# file ending for OS
ending = ".so"
path = "/../src/"

if Sys.iswindows()
    ending = ".dll"
    path = "\\..\\src\\"
elseif Sys.isapple()
    ending = ".dylib"
end

if cpufeature(:AVX2)
    flag = "avx2"
elseif cpufeature(:AVX)
    flag = "avx"
else
    flag = ""
end

lib_path_nfft = string(@__DIR__, path, "libnfftjulia", flag, ending)
lib_path_nfct = string(@__DIR__, path, "libnfctjulia", flag, ending)
lib_path_nfst = string(@__DIR__, path, "libnfstjulia", flag, ending)
lib_path_fastsum = string(@__DIR__, path, "libfastsumjulia", flag, ending)

println( lib_path_nfft )

chmod(lib_path_nfft, filemode(lib_path_nfft) | 0o755)
chmod(lib_path_nfct, filemode(lib_path_nfct) | 0o755)
chmod(lib_path_nfst, filemode(lib_path_nfst) | 0o755)
chmod(lib_path_fastsum, filemode(lib_path_fastsum) | 0o755)