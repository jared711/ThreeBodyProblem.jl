using ThreeBodyProblem
using Test

# circular orbit
rv = [1,0,0,0,1,0]
μ = 1
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:rad)
@test a == 1.0
@test e == 0.0
@test i == 0.0
@test isnan(Ω)
@test isnan(ω)
@test isnan(ν)
@test isnan(Π)
@test isnan(u)
@test l == 0.0
@test ℰ == -0.5
@test ℰ == computeME(rv, μ)

# polar orbit
rv = [1,0,0,0,0,1]
μ = 1
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:rad)
@test a == 1.0
@test e == 0.0
@test i == π/2
@test Ω == 0.0
@test isnan(ω)
@test isnan(ν)
@test isnan(Π)
@test u == 0.0
@test isnan(l)
@test ℰ == -0.5

# Parabolic orbit
rv = [1,0,0,0,√2,0]
μ = 1
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:rad)
@test a == Inf
@test e ≈ 1
@test i == 0.0
@test isnan(Ω)
@test isnan(ω)
@test ν == 0.0
@test Π == 0.0
@test isnan(u)
@test isnan(l)
@test ℰ == 0

# hyperbolic orbit
rv = [1,0,0,0,2,0]
μ = 1
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:rad)
@test a == -0.5
@test e == 3.0
@test i == 0.0
@test isnan(Ω)
@test isnan(ω)
@test ν == 0.0
@test Π == 0.0
@test isnan(u)
@test isnan(l)
@test ℰ == 1.0

# errors
@test_throws ErrorException cart2oe(rv, μ, ang_unit=:foo)

# inverse
rv = [1,0,0,0,√2,√2] # don't use an equatorial or circular orbit (oe2cart can't handle the special cases).
μ = 1
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:deg)
@test oe2cart(a, e, i, Ω, ω, ν, μ, ang_unit=:deg) == rv
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:rad)
@test oe2cart(a, e, i, Ω, ω, ν, μ, ang_unit=:rad) == rv
@test_throws ErrorException oe2cart(a, e, i, Ω, ω, ν, μ, ang_unit=:foo)

### nu2E and E2nu ###
rv = [0,1,2,0.1,1//2,0] # random elliptical orbit
μ = 1
# inverse
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:deg)
@test E2nu(nu2E(ν, e, ang_unit=:deg), e, ang_unit=:deg) ≈ ν 
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:rad)
@test E2nu(nu2E(ν, e, ang_unit=:rad), e, ang_unit=:rad) ≈ ν 
# edge cases
ν = -10
@test 0 <= nu2E(ν, e, ang_unit=:deg) <= 360
ν = 370
@test 0 <= nu2E(ν, e, ang_unit=:deg) <= 360
E = -10
@test 0 <= E2nu(E, e, ang_unit=:deg) <= 360
E = 370
@test 0 <= E2nu(E, e, ang_unit=:deg) <= 360
# errors
@test_throws ErrorException E2nu(E, e, ang_unit=:hey)
@test_throws TypeError E2nu(E, e, ang_unit=1)
@test_throws ErrorException E2nu(E, 1.0)
@test_throws ErrorException nu2E(ν, e, ang_unit=:hey)
@test_throws TypeError nu2E(ν, e, ang_unit=1)
@test_throws ErrorException nu2E(ν, 1.0)


### E2M and M2E ###
rv = [0,1,2,0.1,1//2,0] # random elliptical orbit
μ = 1
# inverse
E, e = 45, 0.5
@test M2E(E2M(E, e, ang_unit=:deg), e, ang_unit=:deg) ≈ E 
E, e = π/4, 0.5
@test M2E(E2M(E, e, ang_unit=:rad), e, ang_unit=:rad) ≈ E 
# edge cases
M = -10
@test 0 <= M2E(M, e, ang_unit=:deg) <= 360
M = 370
@test 0 <= M2E(M, e, ang_unit=:deg) <= 360
E = -10
@test 0 <= E2M(E, e, ang_unit=:deg) <= 360
E = 370
@test 0 <= E2M(E, e, ang_unit=:deg) <= 360
# errors
@test_throws ErrorException E2M(E, e, ang_unit=:hey)
@test_throws TypeError E2M(E, e, ang_unit=1)
@test_throws ErrorException M2E(M, e, ang_unit=:hey)
@test_throws TypeError M2E(M, e, ang_unit=1)
