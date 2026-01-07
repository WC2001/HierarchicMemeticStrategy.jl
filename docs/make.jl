using hms_jl
using Documenter

DocMeta.setdocmeta!(hms_jl, :DocTestSetup, :(using hms_jl); recursive=true)

makedocs(;
    modules=[hms_jl],
    authors="Wiktor Cieślikiewicz",
    sitename="hms_jl.jl",
    format=Documenter.HTML(;
        canonical="https://WC2001.github.io/hms_jl.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Configuration" => "hms_config.md"
    ],
)

deploydocs(;
    repo="github.com/WC2001/hms_jl.jl",
    devbranch="main",
)
