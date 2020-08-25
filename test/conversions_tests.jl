using ThreeBodyProblem
using Test

μ = 0.9
rv = [1;2;3;4;5;6]
@test rv == [1;2;3;4;5;6]


r = [0,0,1]
@test_broken r2latlon(r) == (90,0)
