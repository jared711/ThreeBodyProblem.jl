# make.jl
push!(LOAD_PATH,"../src/")
using Documenter
using ThreeBodyProblem

makedocs(
    sitename = "ThreeBodyProblem.jl",
    # modules = [ThreeBodyProblem],
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "index.md"
        "Examples" => "example.md"
        "Functions" => [
             "constants.md",
             "frames.md",
             "orbitalelements.md",
             "parameters.md",
             "util.md",
             "dynamics.md"
        ]
        "Plotting" => "plot.md"
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
    target = "build"
    deploy_config=Documenter.Travis(),
    push_preview = true,
)

# Let's see if this will work
