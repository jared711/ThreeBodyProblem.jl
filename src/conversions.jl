"""
    rot2inert!(rv,θ)
Synodic (rotating) frame to inertial frame
"""
# function rot2inert!(rv, θ, μ; origin="barycenter")
function rot2inert!(rv, θ)
    cθ, sθ = cos(θ), sin(θ)
    A = [cθ -sθ  0;
         sθ  cθ  0;
          0   0  1]
    B = [-sθ -cθ  0;
          cθ -sθ  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]
    rv[1:6] = C*rv
    return nothing
end

function rot2inert!(rv, θ, p::Array; origin=0)
# function S2I!(rv,t,p::Array)
    μ₁,μ₂,d = p # parameters
    ωₛ = sqrt((μ₁ + μ₂)/d^3);
    cθ, sθ = cos(ωₛ*t), sin(ωₛ*t)
    A = [cθ -sθ  0;
         sθ  cθ  0;
          0   0  1]
    B = [-sθ -cθ  0;
          cθ -sθ  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]
    rv[1:6] = C*rv
    return nothing
end

function r2latlon(r::Array)
    x,y,z = r
    xy = sqrt(x^2 + y^2)
    lat = atand(z,xy)
    lon = atand(y,x)
    return lat, lon
end
# y   r
# |  /
# | /  lon
# |/------- x

"""
"""
function ECEF2ECI(rv_ECEF, gmst)
end

"""
"""
function ECI2ECEF(rv_ECI, gmst)
end

"""
"""
function ECI2SCI(rv_ECI, rv_sun_ECI)
end

"""
"""
function SCI2ECI(rv_SCI, rv_earth_SCI)
end
