using Documenter
using ThreeBodyProblem

DocMeta.setdocmeta!(ThreeBodyProblem, :DocTestSetup, :(using ThreeBodyProblem); recursive = true)

makedocs(
    sitename = "ThreeBodyProblem.jl",
    modules = [ThreeBodyProblem],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", edit_link = "master"),
    pages = [
        "Introduction" => "index.md",
        "Examples" => "example.md",
        "Functions" => [
            "constants.md",
            "frames.md",
            "orbitalelements.md",
            "parameters.md",
            "util.md",
            "dynamics.md",
        ],
        "Plotting" => "plot.md",
    ],
)

deploydocs(
    repo = "github.com/jared711/ThreeBodyProblem.jl.git",
    devbranch = "master",
    push_preview = true,
)
