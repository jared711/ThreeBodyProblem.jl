using ThreeBodyProblem
using LinearAlgebra
# This particular state is from a meteor we sensed
r_enu = [191.0973115173407; 78.7745585024403; 57.3291938159787] # [km]
v_enu = [-52.199447002937305; -21.517740337533645; -15.659828241680069] #[km/s]

rv_enu = [r_enu; v_enu]

# coordinates of observer (ALTAIR radar)
ϕ = 9.3987 # [deg] latitude
λ = 167.482 # [deg] longitude
h = 0.083532 # [km] altitude

# verified against matlab stuff 4/10/21
rv_ecef = enu2ecef(rv_enu, ϕ, λ, h, geodetic=true, ang_unit=:deg)

norm(rv_enu - ecef2enu(rv_ecef, ϕ, λ, h))

# UT1 November 17, 2000 15:00:05.634
date = [2000,11,17,15,0,5.64]
mjd = date2mjd(date)
θ = mjd2gmst(mjd)

rv_eci = ecef2eci(rv_ecef, θ, ang_unit=:rad)
norm(rv_ecef - eci2ecef(rv_eci, θ))

prim,sec,sys = sun_earth()
a,e,i,Ω,ω,ν = cart2oe(rv_eci, 3.986004328969392e+05)
norm(rv_eci - oe2cart(a,e,i,Ω,ω,ν,G*sec.m))

# Convert the calendar date to ephemeris seconds past J2000
et = utc2et(date2str(date))
# Get the position of Mars at `et` w.r.t. Earth
rv_sun_eci = spkezr("sun", et, "J2000", "none", "earth")[1]
rv_sci = eci2sci(rv_eci, rv_sun_eci)
norm(rv_eci - sci2eci(rv_sci, rv_sun_eci))
