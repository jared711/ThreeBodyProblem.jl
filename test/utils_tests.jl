using ThreeBodyProblem
using Test

@test findL1(.9) == [-0.6090351100232029;0;0]
@test findL1([1,2,3]) == [-0.7122547145555811;0;0]
