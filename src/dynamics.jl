using LinearAlgebra

"""
    R2BPdynamics!(rvdot,rv,μ,t)

Compute change of state vector in restricted two-body system. rv is the state in
{km} and μ is the gravitational parameter in {km³/s²}

# Example
```jldoctest
julia> rvdot = zeros(6)
julia> R2BPdynamics!(rvdot,[1e4;0;0;0;1;0],398600,0)
```
"""
function R2BPdynamics!(rvdot,rv,μ,t)  #make sure rv and μ are in km and km³/s²
    r,v = rv[1:3], rv[4:6]
    rvdot[1:3] = v
    rvdot[4:6] = -μ/norm(r)^3 * r   #acceleration [km/s]
    return nothing
end

"""
    CR3BPdynamics!(rvdot,rv,μ,t)

Compute change of state vector in normalized CR3BP. rv is the normalized state
{NON} and μ is the gravitational parameter {NON}.

# Example
```jldoctest
julia> rvdot = zeros(6)
julia> CR3BPdynamics!(rvdot,[0.5;0;0;0;0.5;0],1e-2,0)
```
"""
function CR3BPdynamics!(rvdot,rv,μ,t) #Three body dynamics in Earth/Moon System
    x,y,z,vx,vy,vz = rv
    μ2 = 1-μ;
    r₁³= ((x+μ )^2 + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    r₂³= ((x-μ2)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4]   = x-(μ2*(x+μ)/r₁³) -(μ*(x-μ2)/r₂³) + 2*vy;
    rvdot[5]   = y-(μ2* y   /r₁³) -(μ* y    /r₂³) - 2*vx;
    rvdot[6]   =  -(μ2* z   /r₁³) -(μ* z    /r₂³);
    return nothing
end

"""
    CR3BPdynamics!(rvdot,rv,p::Array,t)

Compute change of state vector in non-normalized restricted three-body system.
rv is the state [r; v] {km; km/s} and p = [μ₁;μ₂;d] {km³/s²; km³/s²; km}
contains the gravitational parameters of the first and second primary bodies and
the distance between them.

# Example
```jldoctest
julia> rvdot = zeros(6)
julia> CR3BPdynamics!(rvdot,[1e4;0;0;0;1;0],[398600;4970;384400],0)
```
"""
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
