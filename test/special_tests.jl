using ThreeBodyProblem
using LinearAlgebra
using Test

sys = sun_earth()

Az = 0.001
Lpt = 2 # Which libration point will I center about
NS = 1
npts = 100
t, rvs, T, Ax = rich3(sys, Az, Lpt, NS, npts)

rv₀, T = differential_corrector(sys, rvs[1,:], myconst=3, tf=T)

Wsp, Wsn, Wup, Wun, Φₜ, ww = invariant_manifolds(sys, rv₀, T, tf=5*T, nPts=20)
