push!(LOAD_PATH,"../src/")

using TrixiShockTube
using Documenter

DocMeta.setdocmeta!(TrixiShockTube, :DocTestSetup, :(using TrixiShockTube); recursive=true)

makedocs(;
    modules=[TrixiShockTube],
    authors="Alex Ames <aames@wisc.edu>",
    repo="https://github.com/stillyslalom/TrixiShockTube.jl/blob/{commit}{path}#{line}",
    sitename="TrixiShockTube.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stillyslalom.github.io/TrixiShockTube.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Function reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/stillyslalom/TrixiShockTube.jl",
    devbranch="main",
)
