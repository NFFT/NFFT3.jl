using CpuId

# file ending for OS
ending = ".so"
path = "/../src/lib/"
glibcversion = ""

if Sys.iswindows()
    ending = ".dll"
    path = "\\..\\src\\lib\\"
elseif Sys.isapple()
    ending = ".dylib"
else
    glibcversion = "glibc2.40/"
    glibc = VersionNumber(unsafe_string(@ccall gnu_get_libc_version()::Cstring))
    # the glibc2.40 builds need glibc 2.35 (x86_64) or glibc 2.40 (ARM64)
    if glibc < (Sys.ARCH === :aarch64 ? v"2.40" : v"2.35")
        glibcversion = "glibc2.22/"
    end
end

# ARM64 (Linux aarch64, Windows ARM64, Apple Silicon) has a single build
if Sys.ARCH === :aarch64
    flag = "ARM64/"
elseif cpufeature(:AVX2)
    flag = "AVX2/"
elseif cpufeature(:AVX)
    flag = "AVX/"
else
    flag = "SSE2/"
end

lib_path_nfft = string(@__DIR__, path, flag, glibcversion, "libnfftjulia", ending)
lib_path_nfct = string(@__DIR__, path, flag, glibcversion, "libnfctjulia", ending)
lib_path_nfst = string(@__DIR__, path, flag, glibcversion, "libnfstjulia", ending)
lib_path_fastsum = string(@__DIR__, path, flag, glibcversion, "libfastsumjulia", ending)

println(lib_path_nfft)

chmod(lib_path_nfft, filemode(lib_path_nfft) | 0o755)
chmod(lib_path_nfct, filemode(lib_path_nfct) | 0o755)
chmod(lib_path_nfst, filemode(lib_path_nfst) | 0o755)
chmod(lib_path_fastsum, filemode(lib_path_fastsum) | 0o755)