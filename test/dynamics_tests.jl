using ThreeBodyProblem
using Test

PRIM, SEC, SYS = earth_moon()
p = [SYS.μ₁, SYS.μ₂, SYS.d]
rv = [0.5, 0, 0, 0, 0.5, 0]
rv_dim = [1e4, 0, 0, 0, 1, 0]
μ = SYS.μ
t = 0

rvdot = zeros(6)
@test CR3BPdynamics!(rvdot, rv, SYS, t) == nothing
@test rvdot == CR3BPdynamics(rv, SYS, t)

rvdot = zeros(6)
@test CR3BPdynamics!(rvdot, rv_dim, p, t) == nothing
@test rvdot == CR3BPdynamics(rv_dim, p, t)

rvdot = zeros(6)
@test CR3BPinert!(rvdot, rv, SYS, t) == nothing
@test rvdot == CR3BPinert(rv, SYS, t)

rvdot = zeros(6)
@test CR3BPinert!(rvdot, rv_dim, p, t) == nothing
@test rvdot == CR3BPinert(rv_dim, p, t)

rvdot = zeros(6)
@test R2BPdynamics!(rvdot, rv_dim, PRIM, t) == nothing
@test rvdot == R2BPdynamics(rv_dim, PRIM, t)

rvdot = zeros(6)
n = 1e-3
@test CWdynamics!(rvdot, rv, n, t) == nothing
@test rvdot == CWdynamics(rv, n, t)
