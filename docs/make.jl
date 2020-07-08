# make.jl
push!(LOAD_PATH,"../src/")
using Documenter
using ThreeBodyProblem

makedocs(
    sitename = "ThreeBodyProblem",
    # modules = [ThreeBodyProblem],
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "index.md"
        "Overview" => "overview.md"
        "Functions" => [
             "orbitalelements.md",
             "util.md"
        ]
        #     "costfunctions.md",
        #     "constraints.md",
        #     "creating_problems.md",
        #     "solving.md"
        # ],
        # "Interfaces" => [
        #     "costfunction_interface.md",
        #     "constraint_interface.md",
        #     "solver_interface.md"
        # ],
        # "Documentation" => [
        #     "model_types.md",
        #     "discretization.md",
        #     "cost_api.md",
        #     "constraint_api.md",
        #     "problem.md",
        #     "solvers.md",
        #     "rotations.md"
        # ]
    ]
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/jared711/ThreeBodyProblem.jl.git",
    deploy_config=Documenter.Travis(),
)
