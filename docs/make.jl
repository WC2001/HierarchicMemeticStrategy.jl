using HierarchicMemeticStrategy
using Documenter

DocMeta.setdocmeta!(HierarchicMemeticStrategy, :DocTestSetup, :(using HierarchicMemeticStrategy); recursive=true)

makedocs(;
    modules=[HierarchicMemeticStrategy],
    authors="Wiktor Cieślikiewicz",
    sitename="HierarchicMemeticStrategy.jl",
    format=Documenter.HTML(;
        canonical="https://WC2001.github.io/HierarchicMemeticStrategy.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Configuration" => "hms_config.md",
        "Optimization.jl Interface" => "optimization_jl.md"
    ],
)

deploydocs(;
    repo="github.com/WC2001/HierarchicMemeticStrategy.jl",
    devbranch="main",
)
