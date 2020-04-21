using ThreeBodyProblem
using Test

@test AP2a(1e6,1e5) == 5.5e5
@test AP2a(1e6,1e6) == 1e6

p = [1;0;0]
rv = [1;3;4]
@test findrP(rv,p) == 5
rv = [rv,rv]
@test findrP(rv,p) == 5
rv = [1;3;4;5;5;7]
@test findrP(rv,p) == 5
rv = [rv,rv]
@test findrP(rv,p) == 5
findrP(rv,p)
