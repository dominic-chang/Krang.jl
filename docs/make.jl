using Krang
using Documenter
using DocumenterVitepress
using Glob
using Literate

# Make the examples using Literate
GENERATED = joinpath(@__DIR__, "..", "examples")
OUTDIR = joinpath(@__DIR__, "src", "examples")

blacklist = ["gpu", "JuKeBOX", "Krang_logo", "enzyme"]
SOURCE_FILES = filter(x-> all(i->!occursin(i, x), blacklist), Glob.glob("*.jl", GENERATED))[1:2]
foreach(fn -> Literate.markdown(fn, OUTDIR, documenter=true), SOURCE_FILES)
MD_FILES = [joinpath("examples", file) for file in readdir(OUTDIR)]

Documenter.DocMeta.setdocmeta!(Krang, :DocTestSetup, :(using Krang); recursive=true)

format = Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dchang10.github.io/Krang.jl",
        edit_link="main",
        assets=String[],
    )
format=DocumenterVitepress.MarkdownVitepress(
    repo = "https://dchang10.github.io/Krang.jl", # this must be the full URL!
    devbranch = "main",
    devurl = "dev",
    #clean_md_output = true
    ;
)

makedocs(;
    sitename="Krang.jl",
    authors="Dominic <dchang3419@hotmail.com> and contributors",
    modules=[Krang],
    repo="https://github.com/dchang10/Krang.jl/blob/{commit}{path}#{line}",
    format=format,
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Examples" => MD_FILES,
        "Theory" => [
            "Raytracing" => [
                "kerr_geodesic_summary.md",
                "time_regularization.md",
            ],
            "Polarization" => [
                "newmann_penrose.md",
                "polarization.md",
            ]
        ],
        "api.md",
    ],
)

deploydocs(;
    repo="github.com/dchang10/Krang.jl",
    devbranch="main",
)
