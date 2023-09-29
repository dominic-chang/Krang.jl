using Kang
using Documenter

DocMeta.setdocmeta!(Kang, :DocTestSetup, :(using Kang); recursive=true)

makedocs(;
    modules=[Kang],
    authors="Dominic <dchang3419@hotmail.com> and contributors",
    repo="https://github.com/dchang10/Kang.jl/blob/{commit}{path}#{line}",
    sitename="Kang.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dchang10.github.io/Kang.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dchang10/Kang.jl",
    devbranch="main",
)
