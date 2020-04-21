using ThreeBodyProblem
using Test

rvdot = zeros(6)
@test R2BPdynamics!(rvdot,[1e4;0;0;0;1;0],398600,0) == nothing
@test rvdot == [0.0;1.0;0.0;-0.003986;-0.0;-0.0]

@test CR3BPdynamics!(rvdot,[0.5;0;0;0;0.5;0],1e-2,0) == nothing
@test rvdot ==  [0.0;0.5;0.0;-2.2645790609160836;0.0;-0.0]
