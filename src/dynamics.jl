using LinearAlgebra

"""
R2BPdynamics(rv,μ,t)

Compute change of state vector in restricted two-body system. rv is the state in
{km} and μ is the gravitational parameter in {km³/s²}
"""
function R2BPdynamics(rv, μ, t)  #make sure rv and μ are in km and km³/s²
    r,v = rv[1:3], rv[4:6]
    rvdot = zeros(6)
    rvdot[1:3] = v
    rvdot[4:6] = -μ/norm(r)^3 * r   #acceleration [km/s]
    return rvdot
end

"""
    R2BPdynamics(rv,prim::Body,t)

Compute change of state vector in restricted two-body system. rv is the state in
{km} and μ is the gravitational parameter in {km³/s²}
"""
function R2BPdynamics(rv, prim::Body, t)  #make sure rv and μ are in km and km³/s²
    return R2BPdynamics(rv, G*prim.m, t)
end


"""
R2BPdynamics!(rvdot,rv,μ,t)

Compute change of state vector in restricted two-body system. rv is the state in
{km} and μ is the gravitational parameter in {km³/s²}
"""
function R2BPdynamics!(rvdot, rv, μ, t)  #make sure rv and μ are in km and km³/s²
    rvdot[:] = R2BPdynamics(rv,μ,t)
    return nothing
end

"""
    R2BPdynamics!(rvdot,rv,prim::Body,t)

Compute change of state vector in restricted two-body system. rv is the state in
{km} and μ is the gravitational parameter in {km³/s²}
"""
function R2BPdynamics!(rvdot, rv, prim::Body, t)  #make sure rv and μ are in km and km³/s²
    rvdot[:] = R2BPdynamics(rv, prim, t)
    return nothing
end

"""
CR3BPdynamics(rv,p::Array,t)

Compute change of state vector in non-normalized restricted three-body system.
rv is the state [r; v] {km; km/s} and p = [μ₁;μ₂;d] {km³/s²; km³/s²; km}
contains the gravitational parameters of the first and second primary bodies and
the distance between them.
"""
function CR3BPdynamics(rv,μ,t) #Three body dynamics in Earth/Moon System
    x,y,z,vx,vy,vz = rv
    r₁³= ((x + μ)^2     + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    r₂³= ((x - 1 + μ)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass

    rvdot = zeros(6)
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4] = -((1 - μ)*(x + μ)/r₁³) - (μ*(x - 1 + μ)/r₂³) + 2*vy + x;
    rvdot[5] = -((1 - μ)*y      /r₁³) - (μ*y          /r₂³) - 2*vx + y;
    rvdot[6] = -((1 - μ)*z      /r₁³) - (μ*z          /r₂³);
    return rvdot
end

function CR3BPdynamics(rv,sys::System,t) #Three body dynamics in Earth/Moon System
    return CR3BPdynamics(rv,sys.μ,t)
end

function CR3BPdynamics(rv,p::Array,t) #Three body dynamics in Earth/Moon System
    x,y,z,vx,vy,vz = rv
    μ₁,μ₂,d = p # parameters
    R₁ = d*μ₂/(μ₁+μ₂)
    R₂ = d*μ₁/(μ₁+μ₂)
    ωₛ = sqrt((μ₁ + μ₂)/d^3) #rotation rate of system
    r₁³= ((x+R₁)^2 + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    r₂³= ((x-R₂)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass

    rvdot = zeros(6)
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4]   = -(μ₁*(x+R₁)/r₁³) - (μ₂*(x-R₂)/r₂³) + 2*ωₛ*vy + ωₛ^2*x;
    rvdot[5]   = -(μ₁*y     /r₁³) - (μ₂*y     /r₂³) - 2*ωₛ*vx + ωₛ^2*y;
    rvdot[6]   = -(μ₁*z     /r₁³) - (μ₂*z     /r₂³);
    return rvdot
end

"""
CR3BPdynamics!(rvdot,rv,p::Array,t)

Compute change of state vector in non-normalized restricted three-body system.
rv is the state [r; v] {km; km/s} and p = [μ₁;μ₂;d] {km³/s²; km³/s²; km}
contains the gravitational parameters of the first and second primary bodies and
the distance between them.
"""
function CR3BPdynamics!(rvdot,rv,μ,t) #Three body dynamics in Earth/Moon System
    rvdot[:] = CR3BPdynamics(rv,μ,t)
    return nothing
end

function CR3BPdynamics!(rvdot,rv,sys::System,t) #Three body dynamics in Earth/Moon System
    rvdot[:] = CR3BPdynamics(rv,sys,t)
    return nothing
end

function CR3BPdynamics!(rvdot,rv,p::Array,t) #Three body dynamics in Earth/Moon System
    rvdot[:] = CR3BPdynamics(rv,p,t)
    return nothing
end


"""
    CR3BPstm(w,μ,t)

Compute change of state vector in normalized CR3BP. w is the concatenation of rv, the
normalized state {NON}, and vec(Φ), the vectorized state transition matrix {NON}, while μ
is the gravitational parameter {NON}.
"""
function CR3BPstm(w,μ,t) #Three body dynamics in Earth/Moon System
    rv = w[1:6]
    Φ = reshape(w[7:42],6,6)
    x,y,z,vx,vy,vz = rv

    r₁³= ((x + μ)^2     + y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    r₂³= ((x - 1 + μ)^2 + y^2 + z^2)^1.5; # distance to m2, smaller mass
    r₁⁵= ((x + μ)^2     + y^2 + z^2)^2.5; # distance to m1, LARGER MASS
    r₂⁵= ((x - 1 + μ)^2 + y^2 + z^2)^2.5; # distance to m2, smaller mass

    omgxx = 1 - (1-μ)*(1/r₁³ - 3*(x + μ)^2/r₁⁵) - μ*(1/r₂³ - 3*(x - 1 + μ)^2/r₂⁵);
    omgxy = 3*(1-μ)*(x + μ)*y/r₁⁵ + 3*μ*(x - 1 + μ)*y/r₂⁵;
    omgxz = 3*(1-μ)*(x + μ)*z/r₁⁵ + 3*μ*(x - 1 + μ)*z/r₂⁵;
    omgyy = 1 - (1-μ)*(1/r₁³ - 3*y^2/r₁⁵) - μ*(1/r₂³ - 3*y^2/r₂⁵);
    omgyz = 3*(1-μ)*y*z/r₁⁵ + 3*μ*y*z/r₂⁵;
    omgzz = - (1-μ)*(1/r₁³ - 3*z^2/r₁⁵) - μ*(1/r₂³ - 3*z^2/r₂⁵);


    F = [   0     0     0     1     0	 0 ;
            0     0     0     0 	1 	 0 ;
            0	  0     0     0     0    1 ;
        omgxx omgxy omgxz     0     2 	 0 ;
        omgxy omgyy omgyz    -2     0 	 0 ;
        omgxz omgyz omgzz     0	    0	 0 ];

    Φdot = F*Φ;
    wdot = zeros(42)
    wdot[1:6] = CR3BPdynamics(rv,μ,t)
    wdot[7:42] = reshape(Φdot, 36, 1);
    return wdot
end

function CR3BPstm(w,sys::System,t)
    return CR3BPstm(w,sys.μ,t)
end

function CR3BPstm!(wdot,w,μ,t)
    wdot[:] = CR3BPstm(w,μ,t)
    return nothing
end

function CR3BPstm!(wdot,w,sys::System,t) #Three body dynamics in Earth/Moon System
    wdot[:] = CR3BPstm(w,sys,t)
    return nothing
end

"""
    CR3BPinert(rvdot,rv,μ,t)
"""
function CR3BPinert(rv,μ,t)
    x,y,z,vx,vy,vz = rv
    r₁ = [x +       μ*cos(t); y +       μ*sin(t); 0]
    r₂ = [x - (1 - μ)*cos(t); y - (1 - μ)*cos(t); 0]
    r₁³ = norm(r₁)^3
    r₂³ = norm(r₂)^3
    rvdot = zeros(6)
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4:6] = -(1 - μ)*r₁/r₁³ - μ*r₂/r₂³
    return rvdot
end

"""
    CR3BPinert(rvdot,rv,sys::System,t)
"""
function CR3BPinert(rv,sys::System,t)
    return CR3BPinert(rv,sys.μ,t)
end

"""
CR3BPinert(rvdot,rv,p::Array,t)
"""
function CR3BPinert(rv,p::Array,t)
    x,y,z,vx,vy,vz = rv
    μ₁,μ₂,d = p # parameters
    R₁,R₂ = computeR1R2(p)
    ωₛ = sqrt((μ₁ + μ₂)/d^3);
    r₁ = [x + R₁*cos(ωₛ*t); y + R₁*sin(ωₛ*t); 0] # distance to m1, LARGER MASS
    r₂ = [x - R₂*cos(ωₛ*t); y - R₂*sin(ωₛ*t); 0] # distance to m2, smaller mass
    r₁³= norm(r₁)^3;
    r₂³= norm(r₂)^3;
    rvdot = zeros(6)
    rvdot[1:3] = rv[4:6];
    rvdot[4:6] = -μ₁*r₁/r₁³ - μ₂*r₂/r₂³;
    return rvdot
end

"""
CR3BPinert!(rvdot,rv,μ,t)
"""
function CR3BPinert!(rvdot,rv,μ,t)
    rvdot[:] = CR3BPinert(rv,μ,t)
    return nothing
end

"""
    CR3BPinert!(rvdot,rv,sys::System,t)
"""
function CR3BPinert!(rvdot,rv,sys::System,t)
    rvdot[:] = CR3BPinert(rv,sys,t)
    return nothing
end

"""
    CR3BPinert!(rvdot,rv,p::Array,t)
"""
function CR3BPinert!(rvdot,rv,p::Array,t)
    rvdot[:] = CR3BPinert(rv,p,t)
    return nothing
end

"""
    CWdynamics!(rvdot,rv,n,t)

Clohessy-Wiltshire equations

Inputs: n (scalar) mean motion
"""
function CWdynamics(rv,n,t)
    x,y,z,vx,vy,vz = rv
    rvdot = zeros(6)
    rvdot[1:3] = rv[4:6]
    rvdot[4] = 2*n*vy + 3*n^2*x
    rvdot[5] = -2*n*vx
    rvdot[6] = -n^2*z
    return rvdot
end

"""
    CWdynamics!(rvdot,rv,n,t)

Clohessy-Wiltshire equations

Inputs: n (scalar) mean motion
"""
function CWdynamics!(rvdot,rv,n,t)
    rvdot[:] = CWdynamics(rv,n,t)
    return nothing
end


"""
    bicircular problem dynamics
    See G. Gómez, C. Simó, J. Llibre, and R. Martínez, Dynamics and mission design near libration points. Vol. II, vol. 3. 2001.
"""
function BCPdynamics(rv, μ, m₃, n₃, t)
    x,y,z,vx,vy,vz = rv

    # Xₑ = μ*cos(t);      Yₑ = μ*sin(t);      Zₑ = 0
    # Xₘ = (μ-1)*cos(t);  Yₘ = (μ-1)*sin(t);  Zₘ = 0
    # Xₛ = a₃*cos(n₃*t);  Yₛ = a₃*sin(n₃*t);  Zₛ = 0

    a₃ = (1+m₃)^(1/3)/n₃^(2/3)
    θ = (1-n₃)*t
    x₃ =  a₃*cos(θ)
    y₃ = -a₃*sin(θ)

    r₁³ = (  (x+μ)^2 +      y^2 + z^2)^1.5; # distance to m1, LARGER MASS
    r₂³ = ((x-1+μ)^2 +      y^2 + z^2)^1.5; # distance to m2, smaller mass
    r₃³ = ( (x-x₃)^2 + (y-y₃)^2 + z^2)^1.5;

    rvdot = zeros(6)
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4] = -(1-μ)*(x+μ)/r₁³ - μ*(x-1+μ)/r₂³ - m₃*(x-x₃)/r₃³ - m₃*cos(θ)/a₃^2 + 2*vy + x;
    rvdot[5] = -(1-μ)  *  y/r₁³ - μ   *   y/r₂³ - m₃*(y-y₃)/r₃³ + m₃*sin(θ)/a₃^2 - 2*vx + y;
    rvdot[6] = -(1-μ)  *  z/r₁³ - μ   *   z/r₂³ - m₃  *   z/r₃³;
    return rvdot
end

function BCPdynamics(rv, sys::BicircularSystem, t)
    return BCPdynamics(rv, sys.μ, sys.m₃, sys.n₃, t)
end

function BCPdynamics!(rvdot, rv, μ, m₃, n₃, t)
    rvdot[:] = BCPdynamics(rv, μ, m₃, n₃, t)
    return nothing
end

function BCPdynamics!(rvdot, rv, sys::BicircularSystem, t)
    rvdot[:] = BCPdynamics(rv, sys, t)
    return nothing
end


"""
    BCPstm(wdot, w, μ, m₃, n₃, t)

Compute change of state vector in normalized Bicircular. w is the concatenation of rv, the
normalized state {NON}, and vec(Φ), the vectorized state transition matrix {NON}, while μ
is the gravitational parameter {NON}.
"""
function BCPstm(wdot,w, μ, m₃, n₃, t) #Three body dynamics in Earth/Moon System
    rv = w[1:6]
    Φ = reshape(w[7:42],6,6)
    x,y,z,vx,vy,vz = rv

    a₃ = (1+m₃)^(1/3)/n₃^(2/3)
    θ = (1-n₃)*t
    x₃ =  a₃*cos(θ)
    y₃ = -a₃*sin(θ)

    r₁ = sqrt(  (x+μ)^2 +      y^2 + z^2) # distance to m1, Larger Mass
    r₂ = sqrt((x-1+μ)^2 +      y^2 + z^2) # distance to m2, smaller mass
    r₃ = sqrt( (x-x₃)^2 + (y-y₃)^2 + z^2) # distance to m3, LARGEST MASS
    r₁³ = r₁^3; r₂³ = r₂^3; r₃³ = r₃^3
    r₁⁵ = r₁^5; r₂⁵ = r₂^5; r₃⁵ = r₃^5

    g11 = 1 - (1-μ)*(1/r₁³ - 3*(x + μ)^2/r₁⁵) - μ*(1/r₂³ - 3*(x - 1 + μ)^2/r₂⁵) - m₃*(1/r₃³ - 3*(x-x₃)^2/r₃⁵)
    g12 = 3*(1-μ)*(x + μ)*y/r₁⁵ + 3*μ*(x - 1 + μ)*y/r₂⁵ + 3*m₃*(x-x₃)*(y-y₃)/r₃⁵
    # g13 = 0 # planar assumption
    g13 = 3*(1-μ)*(x + μ)*z/r₁⁵ + 3*μ*(x - 1 + μ)*z/r₂⁵ + 3*m₃*(x-x₃)*z/r₃⁵
    g22 = 1 - (1-μ)*(1/r₁³ - 3*y^2/r₁⁵) - μ*(1/r₂³ - 3*y^2/r₂⁵) - m₃*(1/r₃³ - 3*(y-y₃)^2/r₃⁵)
    # g23 = 0 # planar assumption
    g23 = 3*(1-μ)*y*z/r₁⁵ + 3*μ*y*z/r₂⁵ + 3*m₃*(y-y₃)*z/r₃⁵
    # g33 = -(1-μ)/r₁³ - μ/r₂³ - m₃/r₃³ # planar assumption
    g33 = -(1-μ)*(1/r₁³ - 3*z^2/r₁⁵) - μ*(1/r₂³ - 3*z^2/r₂⁵) - m₃*(1/r₃³ - 3*z^2/r₃⁵)

    F = [  0   0   0  1   0   0 ;
           0   0   0  0   1   0 ;
           0   0   0  0   0   1 ;
         g11 g12 g13  0   2   0 ;
         g12 g22 g23 -2   0   0 ;
         g13 g23 g33  0	  0   0 ]

    Φdot = F*Φ
    wdot[1:6] = BCPdynamics(rv,μ,m₃,n₃,t)
    wdot[7:42] = reshape(Φdot, 36, 1);
    return wdot
end

function BCPstm(w,sys::BicircularSystem,t) #Three body dynamics in Earth/Moon System
    return BCPstm(w, sys.μ, sys.m₃, sys.n₃ ,t)
end

function BCPstm!(wdot,w,μ,m₃,n₃,t) #Three body dynamics in Earth/Moon System
    wdot[:] = BCPstm(w,μ,m₃,n₃,t)
    return nothing
end

function BCPstm!(wdot,w,sys::BicircularSystem,t) #Three body dynamics in Earth/Moon System
    wdot[:] = BCPstm(w,sys,t)
    return nothing
end
