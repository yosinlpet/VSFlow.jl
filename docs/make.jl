#!/usr/bin/env julia
# File              : docs/make.jl
# Author            : Denis Dumoulin <denis.dumoulin@uclouvain.be>
# Date              : 30.08.2021
# Last Modified Date: 30.08.2021
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
        "Manual" => "manual.md",
        "Library" => "library.md"
    ],
)

deploydocs(;
    repo="github.com/yosinlpet/VSFlow.jl.git",
    devbranch="main"
)
