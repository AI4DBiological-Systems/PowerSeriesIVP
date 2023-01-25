using Documenter
using PowerSeriesIVP

makedocs(
    sitename = "PowerSeriesIVP",
    format = Documenter.HTML(),
    modules = [PowerSeriesIVP]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
