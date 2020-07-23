using LinearAlgebra
include("util.jl")

# """
#     R2BPdynamics!(rvdot,rv,μ,t)
#
# Compute change of state vector in restricted two-body system. rv is the state in
# {km} and μ is the gravitational parameter in {km³/s²}
# """
function R2BPdynamics!(rvdot,rv,μ,t)  #make sure rv and μ are in km and km³/s²
    r,v = rv[1:3], rv[4:6]
    rvdot[1:3] = v
    rvdot[4:6] = -μ/norm(r)^3 * r   #acceleration [km/s]
    return nothing
end

# """
#     CR3BPdynamics(rvdot,rv,μ,t)
#
# Compute change of state vector in normalized CR3BP. rv is the normalized state
# {NON} and μ is the gravitational parameter {NON}.
# """
function CR3BPdynamics(rvdot,rv,μ,t) #Three body dynamics in Earth/Moon System
    x,y,z,vx,vy,vz = rv
    r₁³= ((x + μ)^2     + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    r₂³= ((x - 1 + μ)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4] = -((1 - μ)*(x + μ)/r₁³) - (μ*(x - 1 + μ)/r₂³) + 2*vy + x;
    rvdot[5] = -((1 - μ)*y      /r₁³) - (μ*y          /r₂³) - 2*vx + y;
    rvdot[6] = -((1 - μ)*z      /r₁³) - (μ*z          /r₂³);
    return rvdot
end

# """
#     CR3BPdynamics!(rvdot,rv,μ,t)
#
# Compute change of state vector in normalized CR3BP. rv is the normalized state
# {NON} and μ is the gravitational parameter {NON}.
# """
function CR3BPdynamics!(rvdot,rv,μ,t) #Three body dynamics in Earth/Moon System
    x,y,z,vx,vy,vz = rv
    r₁³= ((x + μ)^2     + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    r₂³= ((x - 1 + μ)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4] = -((1 - μ)*(x + μ)/r₁³) - (μ*(x - 1 + μ)/r₂³) + 2*vy + x;
    rvdot[5] = -((1 - μ)*y      /r₁³) - (μ*y          /r₂³) - 2*vx + y;
    rvdot[6] = -((1 - μ)*z      /r₁³) - (μ*z          /r₂³);
    return nothing
end

# """
#     CR3BPdynamics!(rvdot,rv,p::Array,t)
#
# Compute change of state vector in non-normalized restricted three-body system.
# rv is the state [r; v] {km; km/s} and p = [μ₁;μ₂;d] {km³/s²; km³/s²; km}
# contains the gravitational parameters of the first and second primary bodies and
# the distance between them.
# """
function CR3BPdynamics!(rvdot,rv,p::Array,t) #Three body dynamics in Earth/Moon System
    x,y,z,vx,vy,vz = rv
    μ₁,μ₂,d = p # parameters
    R₁ = d*μ₂/(μ₁+μ₂)
    R₂ = d*μ₁/(μ₁+μ₂)
    ωₛ = sqrt((μ₁ + μ₂)/d^3) #rotation rate of system
    r₁³= ((x+R₁)^2 + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    r₂³= ((x-R₂)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4]   = -(μ₁*(x+R₁)/r₁³) - (μ₂*(x-R₂)/r₂³) + 2*ωₛ*vy + ωₛ^2*x;
    rvdot[5]   = -(μ₁*y     /r₁³) - (μ₂*y     /r₂³) - 2*ωₛ*vx + ωₛ^2*y;
    rvdot[6]   = -(μ₁*z     /r₁³) - (μ₂*z     /r₂³);
    return nothing
end
#
# """
#     CR3BPinert!(rvdot,rv,μ,t)
# """
function CR3BPinert!(rvdot,rv,μ,t)
    x,y,z,vx,vy,vz = rv
    r₁ = [x +       μ*cos(t); y +       μ*sin(t); 0]
    r₂ = [x - (1 - μ)*cos(t); y - (1 - μ)*cos(t); 0]
    r₁³ = norm(r₁)^3
    r₂³ = norm(r₂)^3
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4:6] = -(1 - μ)*r₁/r₁³ - μ*r₂/r₂³
    return nothing
end

# """
#     CR3BPinert!(rvdot,rv,p::Array,t)
# """
# function CR3BPinert!(rvdot,rv,p::Array,t)
#     x,y,z,vx,vy,vz = rv
#     μ₁,μ₂,d = p # parameters
#     R₁,R₂ = findR1R2(p)
#     ωₛ = sqrt((μ₁ + μ₂)/d^3);
#     r₁ = [x + R₁*cos(ωₛ*t); y + R₁*sin(ωₛ*t); 0] # distance to m1, LARGER MASS
#     r₂ = [x - R₂*cos(ωₛ*t); y - R₂*sin(ωₛ*t); 0] # distance to m2, smaller mass
#     r₁³= norm(r₁)^3;
#     r₂³= norm(r₂)^3;
#     rvdot[1:3] = rv[4:6];
#     rvdot[4:6] = -μ₁*r₁/r₁³ - μ₂*r₂/r₂³;
#     return nothing
# end

# """
#     CWdynamics!(rvdot,rv,t)
#
# Clohessy-Wiltshire equations
# """
function CWdynamics!(rvdot,rv,n,t)
    x,y,z,vx,vy,vz = rv
    rvdot[1:3] = rv[4:6]
    rvdot[4] = 2*n*vy + 3*n^2*x
    rvdot[5] = -2*n*vx
    rvdot[6] = -n^2*z
    return nothing
end
