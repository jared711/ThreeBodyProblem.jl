using ThreeBodyProblem
using Test
using LinearAlgebra
# using SPICE
# furnsh("../src/kernels/de440s.bsp")
# furnsh("../src/kernels/naif0012.tls")

sys = sun_earth()
p = [sys.μ₁, sys.μ₂, sys.d]
ωₛ = sqrt((sys.μ₁ + sys.μ₂)/sys.d^3)

### dimensionalize and nondimensionalize ###
rv = [0.5, 0, 0, 0, 0.5, 0] # random initial condition
# in-place
@test isnothing(dimensionalize!(rv, sys))
@test isnothing(nondimensionalize!(rv, sys))
# inverse
rv_dim = dimensionalize(rv, sys)
@test nondimensionalize(rv_dim, sys) == rv

#Lpts = computeLpts(sys)
#tt = 0:π/10:2π
#xx = [[Lpts[1];zeros(3)] for t in tt]
#xx_inert = [rot2inert(xx[i],tt[i],sys) for i in eachindex(tt)]
#@test norm(xx - [inert2rot(xx_inert[i],tt[i],sys) for i = eachindex(tt)]) < 1e-12

### intert2rot and rot2inert ###
# in-place
@test isnothing(inert2rot!(rv, 0, p))
@test isnothing(inert2rot!(rv, 0, sys))
@test isnothing(inert2rot!(rv, 0, sys.μ))
@test isnothing(rot2inert!(rv, 0, p))
@test isnothing(rot2inert!(rv, 0, sys))
@test isnothing(rot2inert!(rv, 0, sys.μ))
# inverse
for θ in LinRange(0,2π,5)
    for origin in [:barycenter, :prim, :sec]
        @test norm(rv_dim - inert2rot(rot2inert(rv_dim, θ/ωₛ, p, origin=origin), θ/ωₛ, p, origin=origin)) < eps()
        @test norm(rv - inert2rot(rot2inert(rv, θ, sys, origin=origin), θ, sys, origin=origin)) < eps()
    end
end
# errors
@test_throws ErrorException inert2rot(rv, 0, sys, origin=1) # ("origin should be :barycenter, :prim, or :sec")
@test_throws ErrorException rot2inert(rv, 0, sys, origin=:hey) # ("origin should be :barycenter, :prim, or :sec")
@test_throws ErrorException inert2rot(rv, 0, p, origin=[]) # ("origin should be :barycenter, :prim, or :sec")
@test_throws ErrorException rot2inert(rv, 0, p, origin="") # ("origin should be :barycenter, :prim, or :sec")


### ecef2eci and eci2ecef ###
# inverse
@test norm(rv - ecef2eci(eci2ecef(rv, π/4), π/4)) < 1e-12
@test norm(rv - ecef2eci(eci2ecef(rv, 45, ang_unit=:deg), 45, ang_unit=:deg)) < 1e-12
# errors
@test_throws TypeError ecef2eci(rv, 0, ang_unit=1)
@test_throws ErrorException ecef2eci(rv, 0, ang_unit=:abc) # ("ang_unit should be :rad or :deg")
@test_throws TypeError eci2ecef(rv, 0, ang_unit=1)
@test_throws ErrorException eci2ecef(rv, 0, ang_unit=:abc) # ("ang_unit should be :rad or :deg")


### eci2sci and sci2eci ###
rv_sun = [10,0,0,0,10,0]
# inverse
@test norm(rv - sci2eci(eci2sci(rv, rv_sun), rv_sun)) < 1e-12

### ecef2enu and enu2ecef ###
ϕ, λ, h = π/4, -π/2, 100
@test norm(rv - ecef2enu(enu2ecef(rv, ϕ, λ, h), ϕ, λ, h)) < 1e-12
@test norm(rv - ecef2enu(enu2ecef(rv, ϕ, λ, h, ang_unit=:rad), ϕ, λ, h, ang_unit=:rad)) < 1e-12
@test norm(rv - ecef2enu(enu2ecef(rv, ϕ, λ, h, geodetic=false), ϕ, λ, h, geodetic=false)) < 1e-12
@test_throws TypeError enu2ecef(rv, ϕ, λ, h, ang_unit=1)
@test_throws ErrorException enu2ecef(rv, ϕ, λ, h, ang_unit=:abc) # ("ang_unit should be :deg or :rad")
@test_throws TypeError ecef2enu(rv, ϕ, λ, h, ang_unit=1)
@test_throws ErrorException ecef2enu(rv, ϕ, λ, h, ang_unit=:abc) # ("ang_unit should be :deg or :rad")

### cart2latlon and latlon2cart ###
# inverse
r = latlon2cart(ϕ, λ, h, ang_unit=:rad)
ϕ2, λ2, h2 = cart2latlon(r, ang_unit=:rad)
@test norm([ϕ2, λ2, h2] - [ϕ, λ, h]) < 1e-12
r = latlon2cart(ϕ, λ, h, ang_unit=:deg)
ϕ2, λ2, h2 = cart2latlon(r, ang_unit=:deg)
@test norm([ϕ2, λ2, h2] - [ϕ, λ, h]) < 1e-12
# errors
@test_throws TypeError latlon2cart(ϕ, λ, h, ang_unit=1)
@test_throws ErrorException latlon2cart(ϕ, λ, h, ang_unit=:abc) # ("ang_unit should be :deg or :rad")
@test_throws TypeError cart2latlon(r, ang_unit=1)
@test_throws ErrorException cart2latlon(r, ang_unit=:abc) # ("ang_unit should be :deg or :rad")

### cart2azel and azel2cart ###
# inverse
az, el, d = cart2azel(r)
r2 = azel2cart(az, el, d)
@test norm(r - r2) < 1e-12
az, el, d = cart2azel(r, ang_unit=:rad)
r2 = azel2cart(az, el, d, ang_unit=:rad)
@test norm(r - r2) < 1e-12
# errors
@test_throws TypeError azel2cart(az, el, d, ang_unit=1)
@test_throws ErrorException azel2cart(az, el, d, ang_unit=:abc) # ("ang_unit should be :deg or :rad")
@test_throws TypeError cart2azel(r, ang_unit=1)
@test_throws ErrorException cart2azel(r, ang_unit=:abc) # ("ang_unit should be :deg or :rad")
