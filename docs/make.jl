using Krang
using Documenter
using DocumenterVitepress
using Glob
using Literate

# Make the examples using Literate
GENERATED = joinpath(@__DIR__, "..", "examples")
OUTDIR = joinpath(@__DIR__, "src", "examples")

blacklist = ["gpu", "JuKeBOX", "Krang_logo", "enzyme"]
SOURCE_FILES = filter(x-> all(i->!occursin(i, x), blacklist), Glob.glob("*.jl", GENERATED))
foreach(fn -> Literate.markdown(fn, OUTDIR, documenter=true), SOURCE_FILES)
MD_FILES = [joinpath("examples", file) for file in readdir(OUTDIR)]

Documenter.DocMeta.setdocmeta!(Krang, :DocTestSetup, :(using Krang); recursive=true)

format=DocumenterVitepress.MarkdownVitepress(
    repo = "https://github.com/dominic-chang/Krang.jl", # this must be the full URL!
    devbranch = "main",
    devurl = "dev",
    #clean_md_output = true
    ;
)

makedocs(;
    sitename="Krang.jl",
    authors="Dominic <dochang@g.harvard.edu> and contributors",
    modules=[Krang],
    repo="https://github.com/dominic-chang/Krang.jl/blob/{commit}{path}#{line}",
    format=format,
    draft = false,
    source = "src",
    build = "build",
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
    repo="github.com/dominic-chang/Krang.jl",
    target = "build", # this is where Vitepress stores its output
    devbranch="main",
    push_preview=true,
)
