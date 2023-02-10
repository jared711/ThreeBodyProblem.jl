"""
     computed1d2(μ)

Compute the non-dimensional, directional distances of each body from the barycenter given the mass
parameter `μ` {NON}.
"""
function computed1d2(μ)
    d₁ = μ
    d₂ = 1-μ
    return d₁, d₂
end

"""
     computed1d2(sys::System)

Compute the non-dimensional distances of each body from the barycenter given CR3BP system
`sys`.
"""
computed1d2(sys::System) = computed1d2(sys.μ)

"""
    computed1d2(p::Array)

Compute the dimensional distances of each body from the barycenter given `p = [μ₁,μ₂,d]`
{km³/s², km³/s², km} which are the gravitational parameters of the first and second primary
bodies and the distance between them.
"""
function computed1d2(p::Array)
    μ₁,μ₂,d = p
    d₁ = d*μ₂/(μ₁+μ₂)
    d₂ = d*μ₁/(μ₁+μ₂)
    return d₁, d₂
end

"""
    computer1r2(r, μ)

Compute the non-dimensional position vectors of the particle from each body given the position
vector `r` {NON} and the mass parameter `μ` {NON}.
"""
function computer1r2(r,μ)
    x,y,z = r[1:3]
    r₁ = sqrt((x + μ)^2      + y^2 + z^2)
    r₂ = sqrt((x - 1 + μ)^2  + y^2 + z^2)
    return r₁, r₂
end

"""
    computer1r2(r, sys::System)

Compute the non-dimensional position vectors of the particle from each body given the position
vector `r` {NON} and the CR3BP system `sys`.
"""
computer1r2(r,sys::System) = computer1r2(r,sys.μ)

"""
    computer1r2(r, p::Array)

Compute the dimensional position vectors of the particle from each body given the position vector 'r' {km}
and `p = [μ₁,μ₂,d]` {km³/s², km³/s², km} which are the gravitational parameters of the first and
second primary bodies and the distance between them.
"""
function computer1r2(r,p::Array)
    x,y,z = r[1:3]
    d₁,d₂ = computed1d2(p)
    r₁ = sqrt((x + d₁)^2 + y^2 + z^2)
    r₂ = sqrt((x - d₂)^2 + y^2 + z^2)
    return r₁, r₂
end

"""
    computeL1(μ; tol=1e-15)

Compute position vector of L1 in a normalized CR3BP given the mass parameter `μ` {NON}.
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

"""
    computeL1(μ; tol=1e-15)

Compute position vector of L1 in a normalized CR3BP given the CR3BP system `sys`.
"""
computeL1(sys::System; tol=1e-15) = computeL1(sys.μ, tol=tol)


"""
    computeL1(p::Array;tol=1e-15)

Compute position vector of L1 in a non-normalized CR3BP given `p = [μ₁,μ₂,d]`
{km³/s², km³/s², km} which are the gravitational parameters of the first and second primary
bodies and the distance between them.
"""
function computeL1(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L1 = computeL1(μ;tol=tol)*d
    return L1
end

"""
    computeL2(μ; tol=1e-15)

Compute position vector of L2 in a normalized CR3BP given the mass parameter `μ` {NON}.
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

Compute position vector of L3 in a normalized CR3BP given the mass parameter `μ` {NON}.
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
    computeL4(μ; tol=1e-15)

Compute position vector of L4 in a normalized CR3BP given the mass parameter `μ` {NON}.
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

Compute position vector of L5 in a normalized CR3BP given the mass parameter `μ` {NON}.
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
    computeΩ(r,μ)

Compute effective potential given normalized position `r` {NON} and mass parameter `μ` {NON}.
"""
function computeΩ(r,μ)
    if length(r) == 2 || length(r) == 4
        x,y = r[1:2]
        z = 0
    elseif length(r) == 3 || length(r) == 6
        x,y,z = r[1:3]
    else
        error("Invalid position vector length.")
    end
    x,y,z = r[1:3]
    r₁,r₂ = computer1r2(r,μ)
    Ω = (x^2 + y^2)/2 + (1-μ)/r₁ + μ/r₂;
    return Ω
end

"""
    computeΩ(r,sys::System)

Compute effective potential given position `r` {NON} and CR3BP system `sys`.
"""
computeΩ(r,sys::System) = computeΩ(r,sys.μ)

"""
    computeΩ(r,p::Array)

Compute effective potential given position `r` {km} and `p = [μ₁,μ₂,d]`
{km³/s²; km³/s²; km}, which are the gravitational parameters of the first and second primary
bodies and the distance between them.
"""
function computeΩ(r,p::Array)
    if length(r) == 2 || length(r) == 4
        x,y = r[1:2]
        z = 0
    elseif length(r) == 3 || length(r) == 6
        x,y,z = r[1:3]
    else
        error("Invalid position vector length.")
    end
    r₁,r₂ = computer1r2(r,p)
    μ₁,μ₂,d = p
    n = sqrt((μ₁ + μ₂)/d^3)
    Ω = ((x^2 + y^2)*n^2)/2 + μ₁/r₁ + μ₂/r₂;
    return Ω
end

"""
    computeUeff(r,μ)

Compute effective potential given normalized position `r` {NON} and mass parameter `μ` {NON}.
"""
computeUeff(r,μ) = computeΩ(r,μ)


"""
    computeUeff(r,sys::System)

Compute effective potential given position `r` {NON} and CR3BP system `sys`.
"""
computeUeff(r,sys::System) = computeΩ(r,sys)

"""
    computeΩ(r,p::Array)

Compute effective potential given position `r` {km} and `p = [μ₁,μ₂,d]`
{km³/s²; km³/s²; km}, which are the gravitational parameters of the first and second primary
bodies and the distance between them.
"""
computeUeff(r,p::Array) = computeΩ(r,p)

"""
    computeC(rv,μ)

Compute Jacobi constant given normalized state rv {NON} and mass parameter μ {NON}
"""
function computeC(rv,μ)
    if length(rv) == 2
        r = [rv[1], rv[2], 0]
        v = 0
    elseif length(rv) == 3
        r = rv[1:3]
        v = 0
    elseif length(rv) == 4
        r = [rv[1], rv[2], 0]
        v = [rv[3], rv[4], 0]
    elseif length(rv) == 6
        r = rv[1:3]
        v = rv[4:6]
    else
        error("rv must be of length 2, 3, 4, or 6")
    end
    v = norm(v)
    Ω = computeΩ(r,μ)
    C = 2*Ω - v^2
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
    if length(rv) == 2
        r = [rv[1], rv[2], 0]
        v = 0
    elseif length(rv) == 3
        r = rv[1:3]
        v = 0
    elseif length(rv) == 4
        r = [rv[1], rv[2], 0]
        v = [rv[3], rv[4], 0]
    elseif length(rv) == 6
        r = rv[1:3]
        v = rv[4:6]
    else
        error("rv must be of length 2, 3, 4, or 6")
    end
    v = norm(v)
    Ω = computeΩ(r,p)
    C = 2*Ω - v^2
    return C
end

"""
    computeC(traj::Vector{Vector{T}} where T<:Real, sys::System)

Compute Jacobi constants of each state rv = [r; v] {NON; NON} on a 
trajectory traj given a system sys
"""
function computeC(traj::Vector{Vector{T}} where T<:Real, sys::System)
    C = [computeC(traj[i],sys) for i ∈ 1:lastindex(traj)]
    return C
end

"""
    computeCLpts(μ)

Compute Jacobi constants of Lagrange Points given mass parameter μ {NON}.
"""
function computeCLpts(μ)
    Lpts = computeLpts(μ)
    CLpts = [computeC(Lpts[i],μ) for i ∈ eachindex(Lpts)]
    return CLpts
end

"""
    computeCLpts(sys::System)

Compute Jacobi constants of Lagrange Points given system.
"""
computeCLpts(sys::System) = computeCLpts(sys.μ)

"""
    computeCLpts(p::Array)

Compute Jacobi constants of Lagrange Points given p = [μ₁,μ₂,d],
which are the gravitational parameters of the first and second primary bodies
[km³/s²] and the distance between them [km].
"""
function computeCLpts(p::Array)
    Lpts = computeLpts(p)
    CLpts = [computeC(Lpts[i],p) for i ∈ eachindex(Lpts)]
    return CLpts
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

    T = aⱼ/a + 2*cosd(i)*sqrt(a/aⱼ*(1 - e^2))
    return T
end

"""
    stability_index(Φ) UNFINISHED

Compute the stability index for a trajectory given its state transition matrix Φ
"""
function stability_index(Φ)
    λ,_ = eigen(Φ) # λ is a vector of eigenvalues and V is a matrix of eigenvectors
    λmax =  maximum(real(λ))
    ν = 1/2*(abs(λmax) + 1/abs(λmax))
    return ν
end

"""
    rotx(θ)

returns the rotation matrix about the x-axis by angle θ {rad}
"""
function rotx(θ)
    R = [1      0       0;
         0 cos(θ) -sin(θ);
         0 sin(θ)  cos(θ)]
    return R
end

"""
    rotxd(θ)

returns the rotation matrix about the x-axis by angle θ {deg}
"""
function rotxd(θ)
    R = [1      0       0;
         0 cosd(θ) -sind(θ);
         0 sind(θ)  cosd(θ)]
    return R
end

"""
    roty(θ)

returns the rotation matrix about the y-axis by angle θ {rad}
"""
function roty(θ)
    R = [cos(θ)  0  sin(θ);
              0  1       0;
        -sin(θ)  0  cos(θ)]
    return R
end

"""
    rotyd(θ)

returns the rotation matrix about the y-axis by angle θ {deg}
"""
function rotyd(θ)
    R = [cosd(θ)  0  sind(θ);
               0  1        0;
        -sind(θ)  0  cosd(θ)]
    return R
end

"""
    rotz(θ)

# returns the rotation matrix about the z-axis by angle θ {rad}
"""
function rotz(θ)
    R = [cos(θ) -sin(θ) 0;
         sin(θ)  cos(θ) 0;
              0       0 1]
    return R
end

"""
    rotzd(θ)

# returns the rotation matrix about the z-axis by angle θ {deg}
"""
function rotzd(θ)
    R = [cosd(θ) -sind(θ) 0;
         sind(θ)  cosd(θ) 0;
               0        0 1]
    return R
end

"""
    rotlatlon(ϕ, λ; ang_unit=:deg)

# returns the rotation matrix to convert from ECEF to local ENU coordinates given latitude and 
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

Converts a date to Modified Julian Date (MJD) given a date in the form [Y,M,D] or [Y,M,D,h,m,s].
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

# returns the Greenwich Mean Sidereal Time (GMST) at a given UT1 time in Modified Julian Date (MJD)
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

# returns the angle θ in the range [0, 360]
"""
function wrapto360(θ)
    while θ < 0;    θ += 360;   end
    while θ > 360;  θ -= 360;   end
    return θ
end

"""
    wrapto180(θ)

# returns the angle θ in the range [-180, 180]
"""
function wrapto180(θ)
    while θ < -180;  θ += 360;   end
    while θ > 180; θ -= 360;   end
    return θ
end

"""
    wrapto2pi(θ)

# returns the angle θ in the range [0, 2π]
"""
function wrapto2pi(θ)
    while θ < 0;    θ += 2π;    end
    while θ > 2π;   θ -= 2π;    end
    return θ
end

"""
    wraptopi(θ)

# returns the angle θ in the range [-π, π]
"""
function wraptopi(θ)
    while θ < -π;   θ += 2π;    end
    while θ > π;    θ -= 2π;    end
    return θ
end

"""
    date2str(date)

# returns a string with the date in the format "YYYY-MM-DD hh:mm:ss"
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

# returns a point cloud of N_desired points in a hemisphere centered around the normal vector
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
    spherical_ring(c::Array, r::Array, α; N=100)

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

# """
# Calculate the ratio of Lagrange point distance from closest primary to distance between two primaries
# """
# function gammaL(sys::System, Lpt::Int)
#     μ = sys.μ

#     # poly1 = [1, -1*(3-μ), (3-2*μ),     -μ,      2*μ,     -μ];
#     # poly2 = [1,    (3-μ), (3-2*μ),     -μ,     -2*μ,     -μ];
#     # poly3 = [1,    (2+μ), (1+2*μ), -(1-μ), -2*(1-μ), -(1-μ)];

#     poly1 = [    -μ,      2*μ,     -μ, (3-2*μ), -1*(3-μ), 1];
#     poly2 = [    -μ,     -2*μ,     -μ, (3-2*μ),    (3-μ), 1];
#     poly3 = [-(1-μ), -2*(1-μ), -(1-μ), (1+2*μ),    (2+μ), 1];

#     rt1 = roots(poly1)
#     rt2 = roots(poly2)
#     rt3 = roots(poly3)

#     Γ = zeros(3)
#     for i=1:5
#         if isreal(rt1[i]) Γ[1]=rt1[i]; end
#         if isreal(rt2[i]) Γ[2]=rt2[i]; end
#         if isreal(rt3[i]) Γ[3]=rt3[i]; end
#     end
#     γ = Γ[Lpt];
# end
