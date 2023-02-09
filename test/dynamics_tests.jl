using ThreeBodyProblem
using LinearAlgebra
using Test

# set up a system and some initial conditions
sys = earth_moon()
p = [sys.μ₁, sys.μ₂, sys.d]
rv = [0.5, 0, 0, 0, 0.5, 0]
rv_dim = [1e4, 0, 0, 0, 1, 0]
μ = sys.μ
t = 0

rvdot = zeros(6)

### CR3BPdynamics ###
# System input
@test isnothing(CR3BPdynamics!(rvdot, rv, sys, t))
@test rvdot == CR3BPdynamics(rv, sys, t)
# μ input
@test isnothing(CR3BPdynamics!(rvdot, rv, sys.μ, t))
@test rvdot == CR3BPdynamics(rv, sys.μ, t)
# parameters input
@test isnothing(CR3BPdynamics!(rvdot, rv_dim, p, t))
@test rvdot == CR3BPdynamics(rv_dim, p, t)

### CR3BPinert ###
# System input
@test isnothing(CR3BPinert!(rvdot, rv, sys, t))
@test rvdot == CR3BPinert(rv, sys, t)
# μ input
@test isnothing(CR3BPinert!(rvdot, rv, sys.μ, t))
@test rvdot == CR3BPinert(rv, sys, t)
# parameters input
@test isnothing(CR3BPinert!(rvdot, rv_dim, p, t))
@test rvdot == CR3BPinert(rv_dim, p, t)

### CR3BPstm ###
wdot = zeros(42)
w = [rv; reshape(I(6),36,1)]
w_dim = [rv_dim; reshape(I(6),36,1)]
# System input
@test isnothing(CR3BPstm!(wdot, w, sys, t))
@test wdot == CR3BPstm(w, sys, t)
# μ input
@test isnothing(CR3BPstm!(wdot, w, sys.μ, t))
@test wdot == CR3BPstm(w, sys.μ, t)
# parameters input
@test isnothing(CR3BPstm!(wdot, w_dim, sys.μ, t))
@test wdot == CR3BPstm(w_dim, sys.μ, t)

### R2BPdynamics ###
@test isnothing(R2BPdynamics!(rvdot, rv_dim, sys.prim, t))
@test rvdot == R2BPdynamics(rv_dim, sys.prim, t)

@test isnothing(R2BPdynamics!(rvdot, rv_dim, sys.μ₁, t))
@test rvdot == R2BPdynamics(rv_dim, sys.μ₁, t)

### CWdynamics ###
n = 1e-3
@test isnothing(CWdynamics!(rvdot, rv, n, t))
@test rvdot == CWdynamics(rv, n, t)

### BCPdynamics ###
bisys = earth_moon_sun()
@test isnothing(BCPdynamics!(rvdot, rv, bisys, t))
@test isnothing(BCPdynamics!(rvdot, rv, bisys.μ, bisys.m3, bisys.n3, t))
@test rvdot == BCPdynamics(rv, bisys, t)

### BCPstm ###
@test isnothing(BCPstm!(wdot, w, bisys.μ, bisys.m3, bisys.n3, t))
@test isnothing(BCPstm!(wdot, w, bisys, t))
@test wdot == BCPstm(w, bisys, t)

### BCPdynamics2 ###
@test isnothing(BCPdynamics2!(rvdot, rv, bisys, t))
@test isnothing(BCPdynamics2!(rvdot, rv, bisys.μ, bisys.m3, bisys.n3, t))
@test rvdot == BCPdynamics2(rv, bisys, t)
