using CpuId

# file ending for OS
ending = ".so"
path = "/../src/lib/"

if Sys.iswindows()
    ending = ".dll"
    path = "\\..\\src\\lib\\"
elseif Sys.isapple()
    ending = ".dylib"
end

if cpufeature(:AVX2)
    flag = "AVX2/"
elseif cpufeature(:AVX)
    flag = "AVX/"
else
    flag = "SSE2/"
end

lib_path_nfft = string(@__DIR__, path, flag, "libnfftjulia", ending)
lib_path_nfct = string(@__DIR__, path, flag, "libnfctjulia", ending)
lib_path_nfst = string(@__DIR__, path, flag, "libnfstjulia", ending)
lib_path_fastsum = string(@__DIR__, path, flag, "libfastsumjulia", ending)

println( lib_path_nfft )

chmod(lib_path_nfft, filemode(lib_path_nfft) | 0o755)
chmod(lib_path_nfct, filemode(lib_path_nfct) | 0o755)
chmod(lib_path_nfst, filemode(lib_path_nfst) | 0o755)
chmod(lib_path_fastsum, filemode(lib_path_fastsum) | 0o755)