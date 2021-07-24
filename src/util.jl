"""
     computeR1R2(μ)

Compute the non-dimensional distances of each body from the barycenter given the mass
parameter, μ
"""
function computeR1R2(μ::Number)
    R₁ = -μ
    R₂ = 1-μ
    return R₁, R₂
end

function computeR1R2(sys::System)
    return computeR1R2(sys.μ)
end

"""
    computeR1R2(p::Array)

Compute the dimensional distances of each body from the barycenter p = [μ₁,μ₂,d], which are
the gravitational parameters of the first and second primary bodies [km³/s²]and the distance
between them [km].
"""
function computeR1R2(p::Array)
    μ₁,μ₂,d = p
    R₁ = d*μ₂/(μ₁+μ₂)
    R₂ = d*μ₁/(μ₁+μ₂)
    return R₁, R₂
end

"""
    computer1r2(rv,μ)

Compute the non-dimensional distances of each body from the barycenter given the mass
parameter, μ
"""
function computer1r2(rv,μ)
    x,y,z = rv[1:3]
    r₁ = (x + μ)^2      + y^2 + z^2
    r₂ = (x - 1 + μ)^2  + y^2 + z^2
    return r₁, r₂
end

function computer1r2(rv,sys::System)
    return computer1r2(rv,sys.μ)
end

"""
    computer1r2(rv,p::Array)

Compute the dimensional position vector of the object from each body given p = [μ₁,μ₂,d],
which are the gravitational parameters of the first and second primary bodies [km³/s²]and
the distance between them [km]
"""
function computer1r2(rv,p::Array)
    x,y,z = rv[1:3]
    R₁,R₂ = computeR1R2(p)
    r₁ = (x + R₁)^2 + y^2 + z^2
    r₂ = (x - R₂)^2 + y^2 + z^2
    return r₁, r₂
end


"""
    computeL1(μ;tol=1e-15)

Compute position vector to L1 in a normalized CR3BP given the mass parameter, μ
"""
function computeL1(μ;tol=1e-15)
    α = (μ/3 * (1 - μ))^(1/3)
    dα = 1
    count = 0
    while abs(dα) > tol
        α₀ = α
        α = (μ*(1 - α)^2 / (3 - 2*μ - α*(3 - μ - α)))^(1/3)
        dα = α - α₀
        count += 1
        count > 100 ? break : nothing
    end
    L1 = [1 - μ - α; 0; 0]
    return L1
end
computeL1(sys::System; tol=1e-15) = computeL1(sys.μ, tol=tol)


"""
    computeL1(p::Array;tol=1e-15)

Compute 3D L1 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function computeL1(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L1 = computeL1(μ;tol=tol)*d
    return L1
end

"""
    computeL2(μ;tol=1e-15)

Compute 3D L2 in a normalized CR3BP given μ, the system mass parameter.
"""
function computeL2(μ;tol=1e-15)
    β = (μ/3 * (1 - μ))^(1/3)
    dβ = 1
    count = 0
    while abs(dβ) > tol
        β₀ = β
        β = (μ*(1 + β)^2 / (3 - 2*μ + β*(3 - μ + β)))^(1/3)
        dβ = β - β₀
        count += 1
        count > 100 ? break : nothing
    end
    L2 = [1 - μ + β; 0; 0]
    return L2
end
computeL2(sys::System; tol=1e-15) = computeL2(sys.μ, tol=tol)


"""
    computeL2(p::Array;tol=1e-15)

Compute 3D L2 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function computeL2(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L2 = computeL2(μ;tol=tol)*d
    return L2
end

"""
    computeL3(μ;tol=1e-15)

Compute 3D L3 in a normalized CR3BP given μ, the system mass parameter.
"""
function computeL3(μ;tol=1e-15)
    γ = -7*μ/12 + 1
    dγ = 1
    count = 0
    while abs(dγ) > tol
        γ₀ = γ
        γ = ((1 - μ)*(1 + γ)^2 / (1 + 2*μ + γ*(2 + μ + γ)))^(1/3)
        dγ = γ - γ₀
        count += 1
        count > 100 ? break : nothing
    end
    L3 = [-μ - γ; 0; 0]
    return L3
end

computeL3(sys::System; tol=1e-15) = computeL3(sys.μ, tol=tol)


"""
    computeL3(p::Array;tol=1e-15)

Compute 3D L3 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function computeL3(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L3 = computeL3(μ;tol=tol)*d
    return L3
end

"""
    computeL4(μ;tol=1e-15)

ComputeL4 3D L4 in a normalized CR3BP given μ, the system mass parameter.
"""
computeL4(μ;tol=1e-15) = [0.5-μ; √3/2; 0]
computeL4(sys::System; tol=1e-15) = computeL4(sys.μ, tol=tol)
"""
    computeL4(p::Array;tol=1e-15)

Compute 3D L4 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function computeL4(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L4 = computeL4(μ;tol=tol)*d
    return L4
end

"""
    computeL5(μ;tol=1e-15)

Compute 3D L5 in a normalized CR3BP given μ, the system mass parameter.
"""
computeL5(μ;tol=1e-15) = [0.5-μ; -√3/2; 0]
computeL5(sys::System; tol=1e-15) = computeL5(sys.μ, tol=tol)

"""
    computeL5(p::Array;tol=1e-15)

Compute 3D L5 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function computeL5(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L5 = computeL5(μ;tol=tol)*d
    return L5
end

"""
    computeLpts(μ;tol=1e-15)

Compute 3D Lagrange points in a normalized CR3BP given μ, the system mass parameter.
"""
function computeLpts(μ; tol=1e-15)
    Lpts = [
        computeL1(μ,tol=tol),
        computeL2(μ,tol=tol),
        computeL3(μ,tol=tol),
        computeL4(μ,tol=tol),
        computeL5(μ,tol=tol),
        ]
    return Lpts
end
"""
    computeLpts(sys::System;tol=1e-15)

Compute 3D Lagrange points in a non-normalized CR3BP given system S.
"""
computeLpts(sys::System; tol=1e-15) = computeLpts(sys.μ, tol=tol)

"""
    computeLpts(p::Array;tol=1e-15)

Compute 3D Lagrange points in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function computeLpts(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    Lpts = computeLpts(μ, tol=tol)*d
    return Lpts
end

"""
    computeUeff(rv,μ)

Compute effective potential given normalized state rv {NON} and mass parameter {NON}
"""
function computeUeff(rv,μ)
    x,y,z = rv[1:3]
    r₁,r₂ = computer1r2(rv,μ)
    Ueff = -(x^2 + y^2)/2 - (1-μ)/r₁ - μ/r₂;
    return Ueff
end
computeUeff(rv,sys::System) = computeUeff(rv,sys.μ)

"""
    computeUeff(rv,p::Array)

Compute effective potential given state rv = [r; v] {km; km/s} and p = [μ₁,μ₂,d],
which are the gravitational parameters of the first and second primary bodies
[km³/s²] and the distance between them [km].
"""
function computeUeff(rv,p::Array)
    x,y,z = rv[1:3]
    R₁,R₂ = computeR1R2(p)
    r₁,r₂ = computer1r2(rv,p)
    μ₁,μ₂,d = p
    ωₛ = sqrt((μ₁ + μ₂)/d^3)
    Ueff = -(x^2 + y^2)*ωₛ^2/2 - μ₁/r₁ - μ₂/r₂;
    return Ueff
end

"""
    computeC(rv,μ)

Compute Jacobi constant given normalized state rv {NON} and mass parameter {NON}
"""
function computeC(rv,μ)
    v = norm(rv[4:6])
    Ueff = computeUeff(rv,μ)
    C = -2*Ueff - v^2
    return C
end
"""
    computeC(rv,sys::System)

Compute Jacobi constant given normalized state rv {NON} and system
"""
computeC(rv,sys::System) = computeC(rv, sys.μ)

"""
    computeC(rv,p::Array)

Compute Jacobi constant given state rv = [r; v] {km; km/s} and p = [μ₁,μ₂,d],
which are the gravitational parameters of the first and second primary bodies
[km³/s²] and the distance between them [km].
"""
function computeC(rv,p::Array)
    μ₁,μ₂,d = p
    v = norm(rv[4:6])
    Ueff = computeUeff(rv,p)
    C = -2*Ueff - v^2
    return C
end



"""
    computeT(rv,sys::System) UNFINISHED (I don't like this name, "tisserand" is better)

Compute Tisserand parameter given orbital elements
"""
function computeT(a,e,i; aⱼ=7.783561990635208e8, ang_unit::Symbol=:deg)
    if ang_unit == :deg
    elseif ang_unit == :rad
        i = rad2deg(i)
    else
        error("ang_unit should be :rad or :deg")
    end

    T = aⱼ/a + s*sqrt(a/aⱼ*(1 - e^2)*cosd(i))
    return T
end

"""
    stability_index(Φ) UNFINISHED

Compute the stability index for a trajectory given its state transition matrix Φ
"""
function stability_index(Φ)
    return 1
end

"""
    rotx(θ)
"""
function rotx(θ)
    R = [1      0       0;
         0 cos(θ) -sin(θ);
         0 sin(θ)  cos(θ)]
    return R
end

"""
    rotxd(θ)
"""
function rotxd(θ)
    R = [1      0       0;
         0 cosd(θ) -sind(θ);
         0 sind(θ)  cosd(θ)]
    return R
end

"""
    roty(θ)
"""
function roty(θ)
    R = [cos(θ)  0  sin(θ);
              0  1       0;
        -sin(θ)  0  cos(θ)]
    return R
end

"""
    rotyd(θ)
"""
function rotyd(θ)
    R = [cosd(θ)  0  sind(θ);
               0  1        0;
        -sind(θ)  0  cosd(θ)]
    return R
end

"""
    rotz(θ)
"""
function rotz(θ)
    R = [cos(θ) -sin(θ) 0;
         sin(θ)  cos(θ) 0;
              0       0 1]
    return R
end

"""
    rotzd(θ)
"""
function rotzd(θ)
    R = [cosd(θ) -sind(θ) 0;
         sind(θ)  cosd(θ) 0;
               0        0 1]
    return R
end

"""
    rotlatlon(ϕ, λ; ang_unit=:deg)
"""
function rotlatlon(ϕ, λ; ang_unit=:deg)
    if ang_unit == :rad
        Φ = rotx(π/2 - ϕ)
        Λ = rotz(π/2 + λ)
    elseif ang_unit == :deg
        Φ = rotxd(90 - ϕ)
        Λ = rotzd(90 + λ)
    else
        error("ang_unit should be :deg or :rad")
    end
    R = Λ*Φ
    return R
end

"""
    date2mjd(ut1_date)
"""
function date2mjd(ut1_date)
    if size(ut1_date) == (3,)
        Y = ut1_date[1]
        M = ut1_date[2]
        D = ut1_date[3]
    elseif size(ut1_date) == (6,)
        Y = ut1_date[1]
        M = ut1_date[2]
        D = ut1_date[3] + ut1_date[4]/24 + ut1_date[5]/24/60 + ut1_date[6]/24/60/60
    else
        error("ut1_date should be [Y,M,D] or [Y,M,D,h,m,s]")
    end

    if M <= 2
        y = Y-1;
        m = M+12;
    else
        y = Y;
        m = M;
    end

    if Y < 1582
        B = -2 + ((y+4716/4)/4) - 1179;
    elseif Y > 1582
        B = y/400 - y/100 + y/4;
    else
        if M < 10
            B = -2 + ((y+4716/4)/4) - 1179;
        elseif M > 10
            B = y/400 - y/100 + y/4;
        else
            if D <= 4
                B = -2 + ((y+4716/4)/4) - 1179;
            elseif D >= 10
                B = y/400 - y/100 + y/4;
            else
                B = NaN;
            end
        end
    end

    ut1_mjd = 365*y - 679004 + floor(B) + floor(30.6001*(m+1)) + D;

end


"""
    mjd2gmst(ut1_mjd; ang_unit::Symbol=:rad)
"""
function mjd2gmst(ut1_mjd; ang_unit::Symbol=:rad)
    θ₀ = 280.4606 # [deg] Greenwich Mean Sidereal time at J2000 Epoch
    ωₑ = 360.9856473 # [deg/day] rotating rate of the Earth
    d = ut1_mjd - 51544.5 # Normalize by epoch (Jan. 1, 2000 12:00h)

    if ang_unit == :rad
        θ₀ = deg2rad(θ₀)
        ωₑ = deg2rad(ωₑ)
        θ = θ₀ + ωₑ*d
        θ = wrapto2pi(θ)
    elseif ang_unit == :deg
        θ = θ₀ + ωₑ*d
        θ = wrapto360(θ)
    else
        error("ang_unit should be :deg or :rad")
    end
    return θ
end

"""
    wrapto360(θ)
"""
function wrapto360(θ)
    while θ < 0;    θ += 360;   end
    while θ > 360;  θ -= 360;   end
    return θ
end

"""
    wrapto180(θ)
"""
function wrapto180(θ)
    while θ < -180;  θ += 360;   end
    while θ > 180; θ -= 360;   end
    return θ
end

"""
    wrapto2pi(θ)
"""
function wrapto2pi(θ)
    while θ < 0;    θ += 2π;    end
    while θ > 2π;   θ -= 2π;    end
    return θ
end

"""
    wraptopi(θ)
"""
function wraptopi(θ)
    while θ < -π;   θ += 2π;    end
    while θ > π;    θ -= 2π;    end
    return θ
end

"""
    date2str(date)
"""
function date2str(date)
    Y = Int(date[1])
    M = Int(date[2])
    D = Int(date[3])
    h = Int(date[4])
    m = Int(date[5])
    s = date[6]
    return string(Y, "-", M, "-", D, " ", h, ":", m, ":", s)
end

"""
    deserno_sphere(N_desired)

    Adapted from Lucas Bury's code implementing the algorithm in Deserno, M. (2004) How to
    Generate Equidistributed points on the Surface of a Sphere
"""
function deserno_sphere(N_desired)

area = 4*pi/N_desired
distance = sqrt(area)

M_theta = round(pi/distance);

d_theta = pi/M_theta

d_phi = area/d_theta;

N_new = 0;
xs = zeros(0)
ys = zeros(0)
zs = zeros(0)
for m in 0:(M_theta-1)

    theta = pi*(m+0.5)/M_theta;
    M_phi = round(2*pi*sin(theta)/d_phi); # not exact

    for n in 0:(M_phi-1)
        Phi = 2*pi*n/M_phi;

        N_new = N_new + 1;

        append!(xs, sin(theta)*cos(Phi))
        append!(ys, sin(theta)*sin(Phi))
        append!(zs, cos(theta))
    end
end

xyz_Sphere = copy(transpose(hcat(xs, ys, zs)))

return xyz_Sphere, N_new
end

"""
deserno_hemisphere(N_desired, normal)
"""
function deserno_hemisphere(N_desired, normal)
# -------------------------------------------------
### Obtaining spherical point cloud
# -------------------------------------------------
xyz_Sphere, N_sphere = deserno_sphere(N_desired*2);

# -------------------------------------------------
### Grab the desired hemisphere
# -------------------------------------------------
### Indices points with a 'z' value >= 0
indices = xyz_Sphere[3,:] .>= 0

### Use logical indexing to grab all points with a 'z' >= 0
xyz_hem = xyz_Sphere[:,indices]
N_new = size(xyz_hem, 2)

# -------------------------------------------------
### If new center is in the -z axis, just flip the signs on the existing hemisphere
# -------------------------------------------------
### In getSphere, the pattern is centered around the 'z' axis
oldNormal = [0, 0, 1]

### Turn input vector to unit vector
normal = normal ./ norm(normal)

### In case the new center is in the opposite direction of the old center, can
### save computational time with this (but also the other algorithm breaks)
if normal == -oldNormal
    xyz_unitHemisphere = zeros(size(xyz_hem))
    for kk in 1:N_new
        xyz_unitHemisphere[:,kk] = [xyz_hem[1,kk], xyz_hem[2,kk],-xyz_hem[3,kk]]
    end
    return xyz_unitHemisphere, N_new
end

# -------------------------------------------------
### Compute the rotation matrix to properly align the new hemisphere
# -------------------------------------------------
### Algorithm to compute a rotation matrix between two vectors
v = cross(oldNormal, normal)
skewSymmetricMat = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
R_old2new = [1 0 0; 0 1 0; 0 0 1] + skewSymmetricMat + skewSymmetricMat*skewSymmetricMat*(1 / (1 + dot(oldNormal,normal)))

### Rotate old hemisphere to new position, centered about normal
xyz_unitHemisphere = zeros(size(xyz_hem))
for kk in 1:N_new
    xyz_unitHemisphere[:,kk] = R_old2new * xyz_hem[:,kk]
end

return xyz_unitHemisphere, N_new
end
"""
    Generates a ring of N vectors with origin c offset by α degrees about central vector r
"""
function spherical_ring(c::Array, r::Array, α; N=100)
    ϕ = 90 - α
    λ = LinRange(360/N,360,N)
    r_ring = [latlon2cart(ϕ,λ[i],norm(r),ang_unit=:deg) for i = 1:N]
    n = [1,0,0]
    r = r/norm(r)
    if r'*n > 0.9
        n[1:3] = [0,1,0]
    end
    p = cross(r,n)
    p = p/norm(p)
    q = cross(r,p)
    q = q/norm(q)
    A = [p q r]
    r_ring = [A*r_ring[i] + c for i = 1:N]
    return r_ring
end

"""
Calculate the ratio of Lagrange point distance from closest primary to distance between two primaries
"""
function gammaL(sys::System, Lpt::Int)
    μ = sys.μ

    # poly1 = [1, -1*(3-μ), (3-2*μ),     -μ,      2*μ,     -μ];
    # poly2 = [1,    (3-μ), (3-2*μ),     -μ,     -2*μ,     -μ];
    # poly3 = [1,    (2+μ), (1+2*μ), -(1-μ), -2*(1-μ), -(1-μ)];

    poly1 = [    -μ,      2*μ,     -μ, (3-2*μ), -1*(3-μ), 1];
    poly2 = [    -μ,     -2*μ,     -μ, (3-2*μ),    (3-μ), 1];
    poly3 = [-(1-μ), -2*(1-μ), -(1-μ), (1+2*μ),    (2+μ), 1];

    rt1 = roots(poly1)
    rt2 = roots(poly2)
    rt3 = roots(poly3)

    Γ = zeros(3)
    for i=1:5
            if isreal(rt1[i]) Γ[1]=rt1[i]; end
            if isreal(rt2[i]) Γ[2]=rt2[i]; end
            if isreal(rt3[i]) Γ[3]=rt3[i]; end
    end
    γ = Γ[Lpt];
end
