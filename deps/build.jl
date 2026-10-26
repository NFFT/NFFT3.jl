using CpuId
using Libdl

# file ending for OS
ending = ".so"
path = "/../src/lib/"
glibcversion = ""

if Sys.iswindows()
    ending = ".dll"
    path = "\\..\\src\\lib\\"
elseif Sys.isapple()
    ending = ".dylib"
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

if !Sys.iswindows() && !Sys.isapple()
    # Try the newer (glibc2.40) build first, falling back to glibc2.22 if it
    # fails to load on this system - mirrors src/NFFT3.jl's glibc selection
    # (see that file for why: guessing from the host's reported glibc
    # *version number* doesn't verify the chosen build actually loads).
    glibcversion = "glibc2.40/"
    probe_path = string(@__DIR__, path, flag, glibcversion, "libnfftjulia", ending)
    handle = Libdl.dlopen(probe_path; throw_error = false)
    if handle === nothing
        glibcversion = "glibc2.22/"
    else
        Libdl.dlclose(handle)
    end
end

lib_path_nfft = string(@__DIR__, path, flag, glibcversion, "libnfftjulia", ending)
lib_path_nfct = string(@__DIR__, path, flag, glibcversion, "libnfctjulia", ending)
lib_path_nfst = string(@__DIR__, path, flag, glibcversion, "libnfstjulia", ending)
lib_path_nfsft = string(@__DIR__, path, flag, glibcversion, "libnfsftjulia", ending)
lib_path_fastsum = string(@__DIR__, path, flag, glibcversion, "libfastsumjulia", ending)

println(lib_path_nfft)

chmod(lib_path_nfft, filemode(lib_path_nfft) | 0o755)
chmod(lib_path_nfct, filemode(lib_path_nfct) | 0o755)
chmod(lib_path_nfst, filemode(lib_path_nfst) | 0o755)
chmod(lib_path_nfsft, filemode(lib_path_nfsft) | 0o755)
chmod(lib_path_fastsum, filemode(lib_path_fastsum) | 0o755)