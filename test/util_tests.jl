using ThreeBodyProblem
using LinearAlgebra
using Test

μ₁ = 398600     # Earth gravitational parameter
μ₂ = 4902       # Moon gravitational paramter
d = 384400      # Avg distance between Earth and Moon
p = [μ₁; μ₂; d]
μ = μ₂/(μ₁ + μ₂)

@test computed1d2(μ) == (μ, 1-μ)
sys = earth_moon()
@test computed1d2(sys) == (sys.μ, 1-sys.μ)
r₁, r₂ = computer1r2([1-sys.μ,0,0], sys)
@test r₁ ≈ 1
@test r₂ < eps()

@test computeL1(sys) != computeL1(sys, tol=1e-3)
@test computeL1(p) != computeL1(p, tol=1e-3)
@test computeL2(sys) != computeL2(sys, tol=1e-3)
@test computeL2(p) != computeL2(p, tol=1e-3)
@test computeL3(sys) != computeL3(sys, tol=1e-3)
@test computeL3(p) != computeL3(p, tol=1e-3)
@test norm(computeL4(sys) - [-sys.μ,0,0]) ≈ 1
@test norm(computeL4(p) - d*[-μ,0,0]) ≈ d
@test norm(computeL5(sys) - [-sys.μ,0,0]) ≈ 1
@test norm(computeL5(p) - d*[-μ,0,0]) ≈ d

@test computeLpts(sys,tol=1e-3) != computeLpts(sys)
@test computeLpts(p,tol=1e-3) != computeLpts(p)

### Test computeUeff and computeΩ ###
# aliases
rv = [1,0,0,0,1,0]
@test computeUeff(rv,p) == computeΩ(rv,p)
@test computeUeff(rv,sys) == computeΩ(rv,sys)
@test computeUeff(rv,μ) == computeΩ(rv,μ)
rv = [1,0,0,1] # planar
@test computeUeff(rv,p) == computeΩ(rv,p)
@test computeUeff(rv,sys) == computeΩ(rv,sys)
@test computeUeff(rv,μ) == computeΩ(rv,μ)
# errors
@test_throws ErrorException computeUeff([1,0,0,0,1],p) # wrong size rv
@test_throws ErrorException computeUeff([1,0,0,0,1],μ) # wrong size rv

### Test computeC ###
# planar and spatial cases
@test computeC([1,0,0,0,1,0],sys) == computeC([1,0,0,1],sys)
@test computeC([1e5,0,0,0,1,0],p) == computeC([1e5,0,0,1],p)
# no velocity
@test computeC([1,0,0],sys) == computeC([1,0],sys)
@test computeC([1e5,0,0],p) == computeC([1e5,0],p)
@test computeC([-μ;0;0],μ) == Inf
# trajectory
rvs = [rv for _ in 1:10]
@test length(computeC(rvs,sys)) == 10
# Lagrange points
CLpts = computeCLpts(sys)
@test CLpts[1] > CLpts[2] > CLpts[3] > CLpts[4] == CLpts[5]
CLpts = computeCLpts(p)
@test CLpts[1] > CLpts[2] > CLpts[3] > CLpts[4] == CLpts[5]
# errors
@test_throws ErrorException computeC([1,0,0,0,1],p) # wrong size rv
@test_throws ErrorException computeC([1,0,0,0,1],μ) # wrong size rv

### test computeT ###
a, e, i = 1, 0.5, 10
computeT(a, e, i) == computeT(a, e, deg2rad(i), ang_unit=:rad)
# errors
@test_throws TypeError computeT(a, e, deg2rad(i), ang_unit=1)
@test_throws ErrorException computeT(a, e, deg2rad(i), ang_unit=:abc)

### stability_index ###
@test stability_index(I(6)) == 1

### rotation matrices ###
@test det(roty(π) - I(3)) == 0
@test det(rotyd(180) - I(3)) == 0
@test det(rotz(π) - I(3)) == 0
@test det(rotzd(180) - I(3)) == 0
@test det(rotx(π) - I(3)) == 0
@test det(rotxd(180) - I(3)) == 0

### test date2mjd ###
mytime = [2000, 11, 17.625065277778]
mytime2 = [2000,11,17,15,0,5.63999999]
mytime3 = [1000, 1, 1]
mytime4 = [1582, 11, 1]
mytime5 = [1582, 10, 1]
mytime6 = [1582, 10, 11]
mytime7 = [1582, 10, 5]
mytime8 = [1582, 9, 5]
mytime9 = [1993, 4, 19]
@test date2mjd(mytime2) == date2mjd(mytime)
@test date2mjd(mytime3) == -314577.0
@test date2mjd(mytime4) == -100823.0
@test date2mjd(mytime5) == -101728.0
@test date2mjd(mytime6) == -100844.0
@test isnan(date2mjd(mytime7))
@test date2mjd(mytime8) == -101754.0
@test date2mjd(mytime9) == 49096.0

@test_throws ErrorException date2mjd(1) # ("ut1_date should be [Y,M,D] or [Y,M,D,h,m,s]")

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

@test rotlatlon(π/2,-π/2,ang_unit=:rad) == I
@test_throws ErrorException rotlatlon(0,0,ang_unit=:foo)

@test mjd2gmst(51544.5,ang_unit=:deg) == 280.4606
@test_throws ErrorException mjd2gmst(0,ang_unit=:foo)

### test deserno_hemisphere ###
xyz_sphere,N = deserno_hemisphere(100,[1,0,0])
@test abs(minimum([norm(xyz_sphere[:,i]) for i = 1:N]) - 1) < 1e-10
@test abs(maximum([norm(xyz_sphere[:,i]) for i = 1:N]) - 1) < 1e-10

### test spherical_ring ###
c = zeros(3)
r = [1,0,0]
α = 10 # deg
r_ring = spherical_ring(c, r, α)
@test norm.(r_ring) ≈ ones(100)