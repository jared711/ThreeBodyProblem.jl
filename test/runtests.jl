using ThreeBodyProblem
using Test

@testset "ThreeBodyProblem.jl" begin
    # Write your own tests here.
    @test AP2a(1e6,1e5) == 5.5e5
    @test AP2a(1e6,1e6) == 1e6
    @test ThreeBodyProblem.greet() == "Hello World"
end
