using Kang
using Documenter
using Glob
using Literate

# Make the examples using Literate
GENERATED = joinpath(@__DIR__, "..", "examples")
OUTDIR = joinpath(@__DIR__, "src", "examples")

SOURCE_FILES = Glob.glob("*.jl", GENERATED)
foreach(fn -> Literate.markdown(fn, OUTDIR, documenter=true), SOURCE_FILES)
MD_FILES = [joinpath("examples", file) for file in readdir(OUTDIR)]

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
        "Getting Started" => "getting_started.md",
        "Examples" => MD_FILES,
        "Theory" => ["kerr_geodesic_summary.md",
            "polarization.md",
            "time_regularization.md",
        ],
        "api.md",
    ],
)

deploydocs(;
    repo="github.com/dchang10/Kang.jl",
    devbranch="main",
)
