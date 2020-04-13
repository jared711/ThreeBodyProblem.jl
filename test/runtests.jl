<<<<<<< HEAD
using SafeTestsets
@safetestset "orbitalelements.jl tests" begin
    include("orbitalelements_tests.jl")
=======
using ThreeBodyProblem
using Test

@testset "ThreeBodyProblem.jl" begin
    # Write your own tests here.
    @test AP2a(1e6,1e5) == 5.5e5
    @test AP2a(1e6,1e6) == 1e6
    @test ThreeBodyProblem.greet() == "Hello World"
>>>>>>> 5f1048c1ca06fbc227676f02dcf6f63025b8fe2f
end
