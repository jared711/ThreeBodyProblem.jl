using ThreeBodyProblem
using Test

μ = 0.9
rv = [1;0;0;0;0;0]
@test rot2inert!(rv,0.) == nothing
@test rv == [1;0;0;0;0;0]
