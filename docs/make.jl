using VSFlow
using Documenter

DocMeta.setdocmeta!(VSFlow, :DocTestSetup, :(using VSFlow); recursive=true)

makedocs(;
    modules=[VSFlow],
    authors="Denis Dumoulin",
    repo="https://github.com/yosinlpet/VSFlow.jl/blob/{commit}{path}#{line}",
    sitename="VSFlow.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yosinlpet.github.io/VSFlow.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Library" => "library.md"
    ],
)

deploydocs(;
    repo="github.com/yosinlpet/VSFlow.jl",
    devbranch="main"
)
