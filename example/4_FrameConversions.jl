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
ϕ, λ, h

rv_ecef = enu2ecef(rv_enu, ϕ, λ, h, geodetic=true, ang_unit=:deg)

# Let's check that we get the same r_enu vector when working backwards
error_enu_ecef = norm(rv_enu - ecef2enu(rv_ecef, ϕ, λ, h))

# UT1 November 17, 2000 15:00:05.634
date = [2000,11,17,15,0,5.64] #[Y,M,D,h,m,s] observation date
mjd = date2mjd(date) # [days] modified julian date
θ = mjd2gmst(mjd) # [rad] Greenwich Mean Sidereal time

rv_eci = ecef2eci(rv_ecef, θ, ang_unit=:rad)
# Let's check that we get the same r_ecef vector when working backwards
error_ecef_eci = norm(rv_ecef - eci2ecef(rv_eci, θ))


sys = sun_earth()
a,e,i,Ω,ω,ν,Π,u,l,ℰ = cart2oe(rv_eci, 3.986004328969392e+05)

# Let's check that we get the same r_eci vector when working backwards
error_eci_oe = norm(rv_eci - oe2cart(a,e,i,Ω,ω,ν,G*sys.sec.m))

