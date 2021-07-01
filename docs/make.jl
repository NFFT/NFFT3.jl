using Documenter, NFFT3

makedocs(
    sitename = "NFFT3.jl",
    format = Documenter.HTML(; prettyurls = false),
    modules = [NFFT3],
    pages = [
        "Home" => "index.md",
        "About" => "about.md",
        "Transformations" =>
            ["NFFT" => "NFFT.md", "NFST" => "NFST.md", "NFCT" => "NFCT.md"],
        "Fastsum" => "fastsum.md",
        "Flags" => "Flags.md",
    ],
)

deploydocs(
    repo = "github.com/NFFT/NFFT3.jl.git",
    devbranch = "main",
)
