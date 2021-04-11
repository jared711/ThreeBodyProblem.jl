using ThreeBodyProblem
using Test

# Altair Meteor data
r_enu = [191.0973115173407; 78.7745585024403; 57.3291938159787] # [km]
v_enu = [-52.199447002937305; -21.517740337533645; -15.659828241680069] #[km/s]

rv_enu = [r_enu; v_enu]

# coordinates of observer (ALTAIR radar)
ϕ = 9.3987 # [deg] latitude
λ = 167.482 # [deg] longitude
h = 0.083532 # [km] altitude

# verified against matlab stuff 4/10/21
rv_ecef = enu2ecef(rv_enu, ϕ, λ, h, geodetic=true, ang_unit=:deg)

@test norm(rv_enu - ecef2enu(rv_ecef, ϕ, λ, h)) < 1e-6
