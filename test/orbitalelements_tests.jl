using ThreeBodyProblem
using Test

@test AP2a(1e6,1e5) == 5.5e5
@test AP2a(1e6,1e6) == 1e6
