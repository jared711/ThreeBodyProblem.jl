"""
    AP2a(rA, rP)

Compute the semimajor axis given apoapsis and periapsis distances rA and rP.
"""
function AP2a(Q, q)
    a = (Q + q)/2
end

"""
    computeq(rv::Array,p::Array)
computeq
"""
function findrP(rv::Array, p::Array=[0,0,0])
    if size(rv[1]) == (6,) || size(rv[1]) == (3,)
        r = [rv[i][1:3] - p for i = 1:length(rv)]
        rP = minimum(norm,r)
        idx = findall(x -> norm(x) == rP, r)
    else
        rP = norm(rv[1:3] - p)
        idx = 1
    end
    return rP, idx
end

function computeME(rv::Array, μ)
    r_ECI = rv[1:3]
    v_ECI = rv[4:6]
    r = norm(r_ECI)
    v = norm(v_ECI)
    return 0.5*v^2 - μ/r
end

"""
    cart2oe(rv::Array, μ; ang_unit::Symbol=:deg)

Convert cartesian coordinates [km;km/s] in the ECI frame to orbital elements
(a, e, i, Ω, ω, ν, Π, u, l) given gravitational parameter μ.

See also: [`oe2cart`](@ref)

# Examples
```jldoctest
julia> cart2oe([1,0,0,0,1,0], 1, ang_unit=:deg)
1
```

"""
function cart2oe(rv::Array, μ; ang_unit::Symbol=:deg)
    # pull apart input vector
    r_ECI = rv[1:3]
    v_ECI = rv[4:6]
    r = norm(r_ECI)
    v = norm(v_ECI)

    # vectors describing orbit orientation
    h_ECI = cross(r_ECI, v_ECI)
    h = norm(h_ECI)     # [km²/s] specific angular momentum
    n_ECI = cross([0,0,1], h_ECI)
    n = norm(n_ECI)
    e_ECI = (1/μ)*((v^2 - μ/r)*r_ECI - (r_ECI'*v_ECI)*v_ECI)
    e = norm(e_ECI)     # [] eccentricity

    # size of orbit
    ℰ = 0.5*v^2 - μ/r   # [km²/s²] specific mechanical energy
    # ℰ = computeME(rv, μ)

    if e == 1
        a = Inf         # [km] semimajor axis
        p = h^2/μ       # [km] semilatus rectum
    else
        a = -μ/(2*ℰ)    # [km] semimajor axis
        p = a*(1 - e^2) # [km] semilatus rectum
    end

    # orientation of orbit
    i = acosd(h_ECI[3]/h)           # [deg] inclination
    Ω = acosd(n_ECI[1]/n)           # [deg] right ascension of the ascending node
    ω = acosd((n_ECI'*e_ECI)/(e*n)) # [deg] argument of periapsis
    ν = acosd((e_ECI'*r_ECI)/(e*r)) # [deg] true anomaly

    # put angles in correct domain
    if n_ECI[2] < 0;        Ω = 360 - Ω;    end
    if e_ECI[3] < 0;        ω = 360 - ω;    end
    if r_ECI'*v_ECI < 0;    ν = 360 - ν;    end

    # special cases
    Π = NaN;    u = NaN;    l = NaN;
    if i == 0 && e != 0     # elliptical equatorial
        Π = acosd(e_ECI[1]/e)           # [deg] longitude of periapsis (Π = Ω + ω)
        if e_ECI[2] < 0;    Π = 360 - Π;    end
    elseif i != 0 && e == 0 # circular inclined
        u = acosd((n_ECI'*r_ECI)/(n*r)) # [deg] argument of latitude (u = ω + ν)
        if r_ECI[3] < 0;    u = 360 - u;    end
    elseif i == 0 && e == 0 # circular equatorial
        l = acosd(r_ECI[1]/r)           # [deg] true latitude (l = Ω + ω + ν)
        if r_ECI[2] < 0;    l = 360 - l;    end
    end

    if ang_unit == :deg
    elseif ang_unit == :rad
        i = deg2rad(i)
        Ω = deg2rad(Ω)
        ω = deg2rad(ω)
        ν = deg2rad(ν)
        Π = deg2rad(Π)
        u = deg2rad(u)
        l = deg2rad(l)
    else
        error("ang_unit should be :rad or :deg")
    end

    return a, e, i, Ω, ω, ν, Π, u, l, ℰ
end

# function cart2oe(, sys::System, body::Symbol; ang_unit::Symbol=:deg)
#     if body == :prim
#         return cart2oe(rot2inert(rv, θ, μ; origin=:barycenter), sys.μ1, ang_unit=ang_unit)
#     elseif body == :sec
#         return cart2oe(rot2inert(rv, θ, μ; origin=:barycenter), sys.μ2, ang_unit=ang_unit)
#     else
#         error("body should be :prim or :sec")
#     end
# end

"""
    oe2cart(a, e, i, Ω, ω, ν, μ; ang_unit::Symbol=:deg)

Convert orbital elements to cartesian coordinates [km;km/s] in the ECI frame given
gravitational parameter μ.

See also: [`cart2oe`](@ref)

# Examples
```jldoctest
julia> oe2cart(1, 0, 0, 0, 0, 0, ang_unit=:deg)
[1,0,0,0,1,0]
```
"""
function oe2cart(a,e,i,Ω,ω,ν,μ; ang_unit::Symbol=:deg)
    if ang_unit == :deg
    elseif ang_unit == :rad
        i = rad2deg(i)
        Ω = rad2deg(Ω)
        ω = rad2deg(ω)
        ν = rad2deg(ν)
    else
        error("ang_unit should be :rad or :deg")
    end

    p = a*(1 - e^2) # [km] semilatus rectum
    r = p/(1 + e*cosd(ν))   # [km] magnitude of r_ECI
    r_PQW = [r*cosd(ν), r*sind(ν), 0];
    v_PQW = sqrt(μ/p)*[-sind(ν), e + cosd(ν), 0];
    rv_PQW = [r_PQW; v_PQW]

    r_ECI = rotzd(Ω)*rotxd(i)*rotzd(ω)*r_PQW;
    v_ECI = rotzd(Ω)*rotxd(i)*rotzd(ω)*v_PQW;
    rv_ECI = [r_ECI; v_ECI]
    return rv_ECI
end

"""
    nu2E(E, e; ang_unit::Symbol=:deg)

Convert true anomaly to eccentric anomaly given orbit eccentricity, e.

See also: [`E2nu`](@ref), [`M2E`](@ref), [`E2M`](@ref)

# Examples
```jldoctest
julia> nu2E(180, 0.5, ang_unit=:deg)
180
```
"""
function nu2E(ν, e; ang_unit::Symbol=:deg)
    if ang_unit == :deg
    elseif ang_unit == :rad
        ν = rad2deg(ν)
    else
        error("ang_unit should be :rad or :deg")
    end

    E = 2*atand(tand(nu/2)*sqrt((1 - e)/(1 + e)));

    while E >= 360; E -= 360;   end
    while E <= 0;   E += 360;   end

    if ang_unit == :deg
    elseif ang_unit == :rad
        E = deg2rad(E)
    end

    return E
end

"""
    E2nu(E, e; ang_unit::Symbol=:deg)

Convert eccentric anomaly to true anomaly given orbit eccentricity, e.

See also: [`nu2E`](@ref), [`M2E`](@ref), [`E2M`](@ref)

# Examples
```jldoctest
julia> E2nu(180, 0.5, ang_unit=:deg)
180
```
"""
function E2nu(E, e; ang_unit::Symbol=:deg)
    if ang_unit == :deg
    elseif ang_unit == :rad
        E = rad2deg(E)
    else
        error("ang_unit should be :rad or :deg")
    end

    ν = 2*atand(sqrt((1 + e)/(1 - e))*tand(E/2));
    while ν >= 360; ν -= 360;   end
    while ν <= 0;   ν += 360;   end

    if ang_unit == :deg
    elseif ang_unit == :rad
        ν = deg2rad(ν)
    end

    return ν
end

"""
    E2M(E, e; ang_unit::Symbol=:deg)

Convert eccentric anomaly to mean anomaly given orbit eccentricity, e.

See also: [`M2E`](@ref)

# Examples
```jldoctest
julia> E2M(180, 0.5, ang_unit=:deg)
180
```
"""
function E2M(E, e; ang_unit::Symbol=:deg)
    if ang_unit == :deg
    elseif ang_unit == :rad
        E = rad2deg(E)
    else
        error("ang_unit should be :rad or :deg")
    end

    M = E - e*sind(E)

    if ang_unit == :deg
    elseif ang_unit == :rad
        M = deg2rad(M)
    end

    return M
end

"""
    M2E(M, e; ang_unit::Symbol=:deg, tol=1e-10, i_max=100)

Convert eccentric anomaly to mean anomaly given orbit eccentricity, e.

See also: [`E2M`](@ref)

# Examples
```jldoctest
julia> M2E(180, 0.5, ang_unit=:deg, tol=1e-10)
180
```
"""
function M2E(E, e; ang_unit::Symbol=:deg, tol=1e-10, i_max=100)
    if ang_unit == :deg
    elseif ang_unit == :rad
        M = rad2deg(M)
    else
        error("ang_unit should be :rad or :deg")
    end

    E = M # initial guess
    δ = -(E - e*sind(E) - M)/(1 - e*cosd(E))
    iter = 0;
    while abs(δ) > abs(tol) && iter < i_max
        E = E + δ
        δ = -(E - e*sind(E) - M)/(1 - e*cosd(E))
        iter += 1
    end

    if ang_unit == :deg
    elseif ang_unit == :rad
        E = deg2rad(E)
    end

    return E
end

"""
    azel2cart(az, el, d=1; ang_unit::Symbol=:deg)

Convert azimuth and elevation angles to cartesian coordinates

See also: [`cart2azel`](@ref)

# Examples
```jldoctest
julia> azel2cart(90, 90, ang_unit=:deg)
[0,0,1]
```
"""
function azel2cart(az, el, d=1; ang_unit::Symbol=:deg)
    if ang_unit == :deg
    elseif ang_unit == :rad
        az = rad2deg(az)
        el = rad2deg(el)
    else
        error("ang_unit should be :rad or :deg")
    end

    # vector in East/North/Up coordinates
    r_ENU = [d*cosd(el)*sind(az);
             d*cosd(el)*cosd(az);
             d*sind(el)];

    return r_ENU
end

"""
    cart2azel(r_ENU; ang_unit::Symbol=:deg)

Convert cartesian coordinates from ENU to azimuth, elevation, and magnitude

See also: [`azel2cart2`](@ref)

# Examples
```jldoctest
julia> cart2azel([1,0,0], ang_unit=:deg)
90
```
"""
function cart2azel(r_ENU; ang_unit::Symbol=:deg)
    d = norm(r_ENU)
    el = asind(r_ENU[3]/d)
    az = atan2d(r_ENU[1], r_ENU[2])

    if ang_unit == :deg
    elseif ang_unit == :rad
        az = deg2rad(az)
        el = deg2rad(el)
    else
        error("ang_unit should be :rad or :deg")
    end

    return az, el, d
end


"""
    delaunay()
    input orbital elements? or rv?
    output delaunay variables
"""
function delaunay()
end
