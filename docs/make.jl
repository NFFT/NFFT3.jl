using Documenter, NFFT3

makedocs(
    sitename = "NFFT3",
    format = Documenter.HTML(; prettyurls=false),
    modules = [NFFT3],
    pages=[
        "Home" => "index.md",
        "About" => "about.md",
        "NFFT3" => "NFFT3.md",
        "Transformations" => [
            "NFFT" => "NFFT.md",
            "NFST" => "NFST.md",
            "NFCT" => "NFCT.md",
        ],
        "Fast Summation Algorithm" => "fastsum.md",
        "Notation" => "notation.md",
    ],
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
