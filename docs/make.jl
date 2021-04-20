using Documenter, NFFT3

makedocs(
    sitename = "NFFT3",
    format = Documenter.HTML(; prettyurls = false),
    modules = [NFFT3],
    pages = [
        "Home" => "index.md",
        "About" => "about.md",
        "Flags" => "Flags.md",
        "Transformations" =>
            ["NFFT" => "NFFT.md", "NFST" => "NFST.md", "NFCT" => "NFCT.md"],
        "Fast Summation Algorithm" => "fastsum.md",
        "Notation" => "notation.md",
    ],
)
