using Documenter
using SpectralCRsUtility

makedocs(
    sitename="SpectralCRsUtility",
    format=Documenter.HTML(),
    modules=[SpectralCRsUtility],
    pages = [
            "Table of Contents" => "index.md",
            "Install"           => "install.md",
            "Basics" => "basics.md",
            "Emission" => "emission.md",
            "API reference"     => "api.md"
            ]
        )


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/LudwigBoess/SpectralCRsUtility.jl.git"
)
