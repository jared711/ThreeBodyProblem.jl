using ThreeBodyProblem
using LinearAlgebra
using Test

μ₁ = 398600     # Earth gravitational parameter
μ₂ = 4902       # Moon gravitational paramter
d = 384400      # Avg distance between Earth and Moon
p = [μ₁; μ₂; d]
μ = μ₂/(μ₁ + μ₂)
@test computeR1R2(μ) == (-μ, 1-μ)
@test computeL1(μ) == [0.8369247061700067;0;0]
@test computeL1(p) == [321713.8570517506;0;0]
@test computeL2(μ) == [1.1556746768583537;0;0]
@test computeL2(p) == [444241.3457843512;0;0]
@test computeL3(μ) == [-1.005061834632099;0;0]
@test computeL3(p) == [-386345.76923257887;0;0]
@test computeL4(μ) == [0.48785136133154233;0.8660254037844386;0]
@test computeL4(p) == [187530.06329584488;332900.16521473817;0]
@test computeL5(μ) == [0.48785136133154233;-0.8660254037844386;0]
@test computeL5(p) == [187530.06329584488;-332900.16521473817;0]

@test computeLpts(μ,tol=1e-3) != computeLpts(μ)
@test computeLpts(p,tol=1e-3) != computeLpts(p)

rv = [1,0,0,0,1,0]
@test computeUeff(rv,p) < 0
@test computeUeff([-μ;0;0],μ) == -Inf

@test_broken computeC(rv,p) > 0
@test_broken computeC([-μ;0;0],μ) == Inf

Φ = ones(6,6)
@test_broken stability_index(ϕ) > 1

sys = earth_moon()
Lpts = computeLpts(sys)
rv = vcat(Lpts[:,1],zeros(3))
C₁ = computeC(rv, sys)

A = [-1. 0.  0.;
      0. 1.  0.;
      0. 0. -1.]
@test det(roty(π) - A) == 0
@test det(rotyd(180) - A) == 0


mytime = [2000, 11, 17.625065277778]
mytime2 = [2000,11,17,15,0,5.63999999]
mytime3 = [1000, 1, 1]
mytime4 = [1582, 11, 1]
mytime5 = [1582, 10, 1]
mytime6 = [1582, 10, 11]
mytime7 = [1582, 10, 5]
@test date2mjd(mytime2) == date2mjd(mytime)
@test date2mjd(mytime3) == -314577.0
@test date2mjd(mytime4) == -100823.0
@test date2mjd(mytime5) == -101728.0
@test date2mjd(mytime6) == -100844.0
@test isnan(date2mjd(mytime7))

@test_throws ErrorException("ut1_date should be [Y,M,D] or [Y,M,D,h,m,s]") date2mjd(1)

@test wrapto360(314) == 314
@test wrapto360(372) == 12
@test wrapto360(721) == 1
@test wrapto360(-1) == 359
@test wrapto360(-361) == 359

@test wrapto180(181) == -179
@test wrapto180(179) == 179
@test wrapto180(-181) == 179

@test wrapto2pi(2π) == 2π
@test wrapto2pi(3π) - π == 0.
@test wrapto2pi(-π) - π == 0.

@test wraptopi(π) == π
@test wraptopi(2π) == 0.
@test wraptopi(-2π) == 0.
