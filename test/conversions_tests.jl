using ThreeBodyProblem
using Test

μ = 0.9
rv = [1;2;3;4;5;6]
@test S2I!(rv,0.) == nothing
@test rv == [1;2;3;2;6;6]
