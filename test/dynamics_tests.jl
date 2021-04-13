using ThreeBodyProblem
using Test

sys = earth_moon()
p = [sys.μ₁, sys.μ₂, sys.d]
rv = [0.5, 0, 0, 0, 0.5, 0]
rv_dim = [1e4, 0, 0, 0, 1, 0]
μ = sys.μ
t = 0

rvdot = zeros(6)
@test CR3BPdynamics!(rvdot, rv, sys, t) == nothing
@test rvdot == CR3BPdynamics(rv, sys, t)

rvdot = zeros(6)
@test CR3BPdynamics!(rvdot, rv_dim, p, t) == nothing
@test rvdot == CR3BPdynamics(rv_dim, p, t)

rvdot = zeros(6)
@test CR3BPinert!(rvdot, rv, sys, t) == nothing
@test rvdot == CR3BPinert(rv, sys, t)

rvdot = zeros(6)
@test CR3BPinert!(rvdot, rv_dim, p, t) == nothing
@test rvdot == CR3BPinert(rv_dim, p, t)

rvdot = zeros(6)
@test R2BPdynamics!(rvdot, rv_dim, sys.prim, t) == nothing
@test rvdot == R2BPdynamics(rv_dim, sys.prim, t)

rvdot = zeros(6)
n = 1e-3
@test CWdynamics!(rvdot, rv, n, t) == nothing
@test rvdot == CWdynamics(rv, n, t)
