# furnsh("src/kernels/de440s.bsp")
# furnsh("src/kernels/naif0012.tls")

"""
    rot2inert(rv, θ, μ; origin=:barycenter)

Convert state vector `rv = [r; v]` {NON; NON} from rotating (synodic) frame to inertial
frame in the normalized CR3BP given time `θ` {NON} and mass parameter `μ` {NON}.
"""
function rot2inert(rv, θ, μ; origin=:barycenter)
    if origin == :barycenter
    elseif origin == :prim
        rv -= [-μ, 0, 0, 0, 0, 0]
    elseif origin == :sec
        rv -= [1-μ, 0, 0, 0, 0, 0]
    end
    cθ, sθ = cos(θ), sin(θ)
    A =  [cθ -sθ  0;
          sθ  cθ  0;
           0   0  1]
    B = [-sθ -cθ  0;
          cθ -sθ  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]
    return C*rv
end

"""
    rot2inert(rv, θ, sys::System; origin=:barycenter)

Convert state vector `rv = [r; v]` {NON; NON} from rotating (synodic) frame to inertial
frame in the normalized CR3BP given time `θ` {NON} and CR3BP system `sys`.
"""
rot2inert(rv, θ, sys::System; origin=:barycenter) = rot2inert(rv, θ, sys.μ, origin=origin)

"""
    rot2inert(rv, θ, p::Array; origin=:barycenter)

Convert state vector `rv = [r; v]` {NON; NON} from rotating (synodic) frame to inertial
frame in the non-normalized CR3BP given time `θ` {NON} and p = [μ₁;μ₂;d]`
{km³/s²; km³/s²; km}, which contains the gravitational parameters of the first and second
primary bodies as well as the distance between them.
"""
function rot2inert(rv, θ, p::Array; origin=:barycenter)
    μ₁,μ₂,d = p # parameters

    if origin == :barycenter
    elseif origin == :prim
        rv -= [-d*μ₂/(μ₁ + μ₂), 0, 0, 0, 0, 0]
    elseif origin == :sec
        rv -= [d*(1-μ₂/(μ₁ + μ₂)), 0, 0, 0, 0, 0]
    end

    ωₛ = sqrt((μ₁ + μ₂)/d^3);
    cθ, sθ = cos(ωₛ*θ), sin(ωₛ*θ)
    A = [cθ -sθ  0;
         sθ  cθ  0;
          0   0  1]
    B = [-sθ -cθ  0;
          cθ -sθ  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]
    return C*rv
end

"""
    rot2inert!(rv, θ, μ)

In-place version of rot2inert(rv, θ, μ)
"""
function rot2inert!(rv, θ, μ; origin=:barycenter)
    rv[1:6] = rot2inert(rv, θ, μ, origin=origin)
    return nothing
end

"""
    rot2inert!(rv, θ, sys::System)

In-place version of rot2inert(rv, θ, sys::System)
"""
function rot2inert!(rv, θ, sys::System; origin=:barycenter)
    rv[1:6] = rot2inert(rv, θ, sys, origin=origin)
    return nothing
end

"""
    rot2inert!(rv, θ, p::Array)

In-place version of rot2inert(rv, θ, p::Array)
"""
function rot2inert!(rv, θ, p::Array; origin=0)
    rv[1:6] = rot2inert(rv, θ, p, origin=origin)
    return nothing
end

"""
    inert2rot(rv, θ, μ)

Inertial frame to rotating (synodic) frame
"""
function inert2rot(rv, t, μ; origin=:barycenter)
    ct, st = cos(t), sin(t)
    A =  [ct  st  0;
         -st  ct  0;
           0   0  1]
    # A = rotz(-t)
    B = [-st  ct  0;
         -ct -st  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]
    if origin == :barycenter
    elseif origin == :prim
        rv += [-μ, 0, 0, 0, 0, 0]
    elseif origin == :sec
        rv += [1-μ, 0, 0, 0, 0, 0]
    end
    return C*rv
end

"""
    inert2rot(rv, θ, sys::System)

Inertial frame to rotating (synodic) frame
"""
inert2rot(rv, t, sys::System; origin=:barycenter) = inert2rot(rv, t, sys.μ, origin=origin)

"""
    inert2rot(rv, θ, p::Array)

Inertial frame to rotating (synodic) frame
"""
function inert2rot(rv, t, p::Array; origin=:barycenter)
# function S2I!(rv,t,p::Array)
    μ₁,μ₂,d = p # parameters

    ωₛ = sqrt((μ₁ + μ₂)/d^3);
    ct, st = cos(ωₛ*t), sin(ωₛ*t)
    A = [ct  st  0;
        -st  ct  0;
          0   0  1]
    # A = rotz(ωₛ*t)
    B = [-st  ct  0;
         -ct -st  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]

    if origin == :barycenter
    elseif origin == :prim
        rv += [-d*μ₂/(μ₁ + μ₂), 0, 0, 0, 0, 0]
    elseif origin == :sec
        rv += [d*(1-μ₂/(μ₁ + μ₂)), 0, 0, 0, 0, 0]
    end

    return C*rv
end

"""
    inert2rot!(rv, t, μ)

In-place version of `inert2rot(rv, t, μ)`.
"""
function inert2rot!(rv, t, μ; origin=:barycenter)
    # function inert2rot!(rv, t)
    rv[1:6] = inert2rot(rv, t, μ, origin=origin)
    return nothing
end

"""
    inert2rot!(rv, t, sys::System)

In-place version of `inert2rot(rv, t, sys::System)`.
"""
function inert2rot!(rv, t, sys::System; origin=:barycenter)
    rv[1:6] = inert2rot(rv, t, sys, origin=origin)
    return nothing
end

"""
    inert2rot!(rv, t, p::Array)

In-place version of `inert2rot(rv, t, p::Array)`.
"""
function inert2rot!(rv, t, p::Array; origin=0)
    rv[1:6] = inert2rot(rv, t, p, origin=origin)
    return nothing
end


"""
    ecef2eci(rv_ecef, θ; ω=7.292115373194e-5, ang_unit::Symbol=:rad)

Convert state vector `rv = [r; v]` {km; km/s} from Earth-Centered Earth-Fixed (ECEF) frame
to Earth-Centered Inertial (ECI) frame given Greenwich Mean Sidereal Time (GMST) `θ` {rad}.
"""
function ecef2eci(rv_ecef, θ; ω=7.292115373194e-5, ang_unit::Symbol=:rad)
    if ang_unit == :deg
        θ = deg2rad(θ)
    elseif ang_unit == :rad
    else
        error("ang_unit should be :rad or :deg")
    end


    R = rotz(θ)
    O = zeros(3,3)
    R_dot = ω*[-sin(θ) -cos(θ) 0;
                cos(θ) -sin(θ) 0;
                     0       0 0]

    rv_eci = [R O; R_dot R]*rv_ecef

    return rv_eci
end

"""
    ecef2eci(rv_ecef, θ; ω=7.292115373194e-5, ang_unit::Symbol=:rad)

Convert state vector `rv = [r; v]` {km; km/s} from Earth-Centered Inertial (ECI) frame to
Earth-Centered Earth-Fixed (ECEF) frame given Greenwich Mean Sidereal Time (GMST) `θ` {rad}.
"""
function eci2ecef(rv_eci, θ; ω=7.292115373194e-5, ang_unit::Symbol=:rad)
    if ang_unit == :deg
        θ = deg2rad(θ)
    elseif ang_unit == :rad
    else
        error("ang_unit should be :rad or :deg")
    end

    R = rotz(-θ)
    O = zeros(3,3)
    R_dot = ω*[sin(-θ) cos(-θ) 0;
              -cos(-θ) sin(-θ) 0;
                     0       0 0]

    rv_ecef = [R O; R_dot R]*rv_eci

    return rv_ecef
end

"""
    eci2sci(rv_eci, rv_sun_eci; ε=23.43929)
"""
function eci2sci(rv_eci, rv_sun_eci; ε=23.439292)
    R = rotx(-ε)
    O = zeros(3,3)
    rv_sci = [R O; O R]*(rv_eci-rv_sun_eci)
    return rv_sci
end

"""
    sci2eci(rv_sci, rv_sun_eci; ε=23.439292)
"""
function sci2eci(rv_sci, rv_sun_eci; ε=23.439292)
    R = rotx(ε)
    O = zeros(3,3)
    rv_eci = [R O; O R]*rv_sci + rv_sun_eci
    return rv_eci
end

"""
    enu2ecef(rv_enu, ϕ, λ, h=0; geodetic=true, ang_unit=:deg)
Converts state from ENU frame to inertial frame
"""
function enu2ecef(rv_enu, ϕ, λ, h=0; geodetic=true, ang_unit::Symbol=:deg, e_earth=0.0818, r_earth=6.378136e3)
    d = h + r_earth # [km] radius of the earth plus altitude
    if geodetic # geodetic latitude
        if ang_unit == :rad
            N = r_earth/sqrt(1-(e_earth^2)*sin(λ)^2)
        elseif ang_unit == :deg
            N = r_earth/sqrt(1-(e_earth^2)*sind(λ)^2)
        else
            error("ang_unit should be :deg or :rad")
        end
        r_radar_ecef = latlon2cart(ϕ, λ, 1, ang_unit=ang_unit)
        geodetic_factors = [N + h; N + h; N*(1-e_earth^2) + h]
        r_radar_ecef = r_radar_ecef.*geodetic_factors
    else # geocentric latitude and longitude
        r_radar_ecef = latlon2cart(ϕ, λ, d, ang_unit=ang_unit)
    end

    R = rotlatlon(ϕ, λ, ang_unit=ang_unit)

    r_ecef = R*rv_enu[1:3] + r_radar_ecef
    v_ecef = R*rv_enu[4:6]
    rv_ecef = [r_ecef; v_ecef]
    return rv_ecef
end

"""
    ecef2enu(rv_ecef, ϕ, λ, h=0; geodetic=true, ang_unit::Symbol=:deg, e_earth=0.0818, r_earth=6.378136e3)
Converts state from ENU frame to inertial frame
"""
function ecef2enu(rv_ecef, ϕ, λ, h=0; geodetic=true, ang_unit::Symbol=:deg, e_earth=0.0818, r_earth=6.378136e3)
    d = h + r_earth # [km] radius of the earth plus altitude

    R = rotlatlon(ϕ, λ, ang_unit=ang_unit)

    if geodetic # geodetic latitude
        if ang_unit == :rad
            N = r_earth/sqrt(1-(e_earth^2)*sin(λ)^2)
        elseif ang_unit == :deg
            N = r_earth/sqrt(1-(e_earth^2)*sind(λ)^2)
        else
            error("ang_unit should be :deg or :rad")
        end
        r_radar_ecef = latlon2cart(ϕ, λ, 1, ang_unit=ang_unit)
        geodetic_factors = [N + h; N + h; N*(1-e_earth^2) + h]
        r_radar_ecef = r_radar_ecef.*geodetic_factors
    else # geocentric latitude and longitude
        r_radar_ecef = latlon2cart(ϕ, λ, d, ang_unit=ang_unit)
    end
    r_enu = R'*(rv_ecef[1:3] - r_radar_ecef)
    v_enu = R'*rv_ecef[4:6]
    rv_enu = [r_enu; v_enu]
    return rv_enu
end

"""
    latlon2cart(ϕ, λ, d; ang_unit::Symbol=:deg)
"""
function latlon2cart(ϕ, λ, d; ang_unit::Symbol=:deg)
    if ang_unit == :rad
        r = d*[cos(ϕ)*cos(λ), cos(ϕ)*sin(λ), sin(ϕ)]
    elseif ang_unit == :deg
        r = d*[cosd(ϕ)*cosd(λ), cosd(ϕ)*sind(λ), sind(ϕ)]
    else
        error("ang_unit should be :deg or :rad")
    end
    return r
end

"""
    cart2latlon(r::Array)
"""
function cart2latlon(r::Array, ang_unit::Symbol=:deg)
    x,y,z = r
    xy = sqrt(x^2 + y^2)
    d = norm(r)
    if ang_unit == :rad
        ϕ = asin(z/d)
        λ = atan(y,x)
    elseif ang_unit == :deg
        ϕ = asind(z/d)
        λ = atand(y,x)
    else
        error("ang_unit should be :deg or :rad")
    end
    return ϕ, λ, d
end
# y   r
# |  /
# | /  lon
# |/------- x

function eci2CR3BP()
end

function CR3BP2eci()
end

function CR3BP2cart()
end

function sci2CR3BP()
end

function CR3BP2sci()
end

"""
    dimensionalize(rv, sys::System)

Add dimensions to normalized state vector `rv`
"""
function dimensionalize(rv, sys::System)
    return vcat(rv[1:3]*sys.RUNIT, rv[4:6]*sys.VUNIT)
end

"""
    dimensionalize!(rv, sys::System)

In-place version of dimensionalize(rv, sys::System)
"""
function dimensionalize!(rv, sys::System)
    rv[1:6] = dimensionalize(rv,sys)
    return nothing
end

"""
    nondimensionalize(rv, sys::System)

Normalize state vector `rv`
"""
function nondimensionalize(rv, sys::System)
    return vcat(rv[1:3]/sys.RUNIT, rv[4:6]/sys.VUNIT)
end

"""
    nondimensionalize!(rv, sys::System)

In-place version of nondimensionalize(rv, sys::System)
"""
function nondimensionalize!(rv, sys::System)
    rv[1:6] = nondimensionalize(rv,sys)
    return nothing
end
