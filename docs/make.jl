using Documenter
using NFFT3

makedocs(
    sitename = "NFFT3",
    format = Documenter.HTML(),
    modules = [NFFT3]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
