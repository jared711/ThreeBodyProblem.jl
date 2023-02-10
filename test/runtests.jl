using SafeTestsets

using Test, Documenter, ThreeBodyProblem
# doctest(ThreeBodyProblem)

@safetestset "dynamics.jl tests" begin
    include("dynamics_tests.jl")
end

@safetestset "frames.jl tests" begin
    include("frames_tests.jl")
end

@safetestset "orbitalelements.jl tests" begin
    include("orbitalelements_tests.jl")
end

@safetestset "parameters.jl tests" begin
    include("parameters_tests.jl")
end

@safetestset "plot.jl tests" begin
    include("plot_tests.jl")
end

# @safetestset "special.jl tests" begin
#     include("special_tests.jl")
# end

@safetestset "util.jl tests" begin
    include("util_tests.jl")
end
