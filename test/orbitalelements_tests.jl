using ThreeBodyProblem
using Test

@test AP2a(1e6,1e5) == 5.5e5
@test AP2a(1e6,1e6) == 1e6

p = [1;0;0]
rv = [1;3;4]
@test_broken findrP(rv,p) == 5
rv = [rv,rv]
@test_broken findrP(rv,p) == 5
rv = [1;3;4;5;5;7]
@test_broken findrP(rv,p) == 5
rv = [rv,rv]
@test_broken findrP(rv,p) == 5

# @test_broken M2e

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

@test_throws ErrorException("ang_unit should be :rad or :deg") cart2oe(rv, μ, ang_unit=:foo)

rv = [1,0,0,0,√2,√2]
μ = 1
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:deg)
@test oe2cart(a, e, i, Ω, ω, ν, μ, ang_unit=:deg) == rv
a, e, i, Ω, ω, ν, Π, u, l, ℰ = cart2oe(rv, μ, ang_unit=:rad)
@test oe2cart(a, e, i, Ω, ω, ν, μ, ang_unit=:rad) == rv
@test_throws ErrorException("ang_unit should be :rad or :deg") oe2cart(a, e, i, Ω, ω, ν, μ, ang_unit=:foo)
