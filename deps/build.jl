# file ending for OS
ending = ".so"

if Sys.iswindows()
    ending = ".dll"
elseif Sys.isapple()
    ending = ".dylib"
end

lib_path_nfft = string(@__DIR__, "/../src", "/libnfftjulia", ending)
lib_path_nfct = string(@__DIR__, "/../src", "/libnfctjulia", ending)
lib_path_nfst = string(@__DIR__, "/../src", "/libnfstjulia", ending)
lib_path_fastsum = string(@__DIR__, "/../src", "/libfastsumjulia", ending)

println( lib_path_nfft )

chmod(lib_path_nfft, filemode(lib_path_nfft) | 0o755)
chmod(lib_path_nfct, filemode(lib_path_nfct) | 0o755)
chmod(lib_path_nfst, filemode(lib_path_nfst) | 0o755)
chmod(lib_path_fastsum, filemode(lib_path_fastsum) | 0o755)