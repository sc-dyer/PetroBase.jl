using PetroBase
using Documenter

DocMeta.setdocmeta!(PetroBase, :DocTestSetup, :(using PetroBase); recursive=true)

makedocs(;
    modules=[PetroBase],
    authors="Sabastien Dyer <scdyer@uwaterloo.ca> and contributors",
    repo="https://github.com/sc-dyer/PetroBase.jl/blob/{commit}{path}#{line}",
    sitename="PetroBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://sc-dyer.github.io/PetroBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sc-dyer/PetroBase.jl",
    devbranch="main",
)
