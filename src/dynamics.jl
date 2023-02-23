"""
    R2BPdynamics(rv, μ, t)

Compute time derivative of state vector in the restricted two-body system. `rv` is the state
vector `[r; v]` {km; km/s}, `μ` is the gravitational parameter {km³/s²}, and `t` is time
{s}.
"""
function R2BPdynamics(rv, μ, t)  #make sure rv and μ are in km and km³/s²
    r,v = rv[1:3], rv[4:6]
    rvdot = zeros(6)
    rvdot[1:3] = v
    rvdot[4:6] = -μ/norm(r)^3 * r   #acceleration [km/s]
    return rvdot
end

"""
    R2BPdynamics(rv, prim::Body, t)

Compute time derivative of state vector in the restricted two-body system. `rv` is the state
vector `[r; v]` {km; km/s}, `prim` is the central body, and `t` is time {s}.
"""
function R2BPdynamics(rv, prim::Body, t)  #make sure rv and μ are in km and km³/s²
    return R2BPdynamics(rv, G*prim.m, t)
end

"""
    R2BPdynamics!(rvdot, rv, μ, t)

In-place version of `R2BPdynamics(rvdot, rv, μ, t)`.
"""
function R2BPdynamics!(rvdot, rv, μ, t)  #make sure rv and μ are in km and km³/s²
    rvdot[:] = R2BPdynamics(rv,μ,t)
    return nothing
end

"""
    R2BPdynamics!(rvdot, rv, prim::Body, t)

In-place version of `R2BPdynamics(rvdot, rv, prim::Body, t)`.
"""
function R2BPdynamics!(rvdot, rv, prim::Body, t)  #make sure rv and μ are in km and km³/s²
    rvdot[:] = R2BPdynamics(rv, prim, t)
    return nothing
end

"""
    CR3BPdynamics(rv, μ, t)

Compute time derivative of state vector `rv = [r; v]` {NON, NON} in the rotating frame of
the normalized CR3BP where `μ` is the CR3BP mass parameter μ₂/(μ₁+μ₂) {NON} and `t` is time
{NON}.
"""
function CR3BPdynamics(rv,μ,t) #Three body dynamics in Earth/Moon System
    x,y,z,ẋ,ẏ,ż = rv
    r₁³= ((x + μ)^2     + y^2 + z^2)^1.5 # distance to m1, LARGER MASS
    r₂³= ((x - 1 + μ)^2 + y^2 + z^2)^1.5 # distance to m2, smaller mass

    rvdot = zeros(6)
    rvdot[1:3] = [ẋ;ẏ;ż]
    rvdot[4] = -((1 - μ)*(x + μ)/r₁³) - (μ*(x - 1 + μ)/r₂³) + x + 2*ẏ
    rvdot[5] = -((1 - μ)*y      /r₁³) - (μ*y          /r₂³) + y - 2*ẋ
    rvdot[6] = -((1 - μ)*z      /r₁³) - (μ*z          /r₂³)
    return rvdot
end

"""
    CR3BPdynamics(rv, sys::System, t)

Compute time derivative of state vector `rv = [r; v]` {NON, NON} in the rotating frame of
the normalized CR3BP where `sys` is the CR3BP system and `t` is time {NON}.
"""
function CR3BPdynamics(rv,sys::System,t) #Three body dynamics in Earth/Moon System
    return CR3BPdynamics(rv,sys.μ,t)
end

"""
    CR3BPdynamics(rv, p::Array, t)

Compute time derivative of state vector `rv = [r; v]` {km, km/s} in the rotating frame of
the non-normalized CR3BP where `p = [μ₁;μ₂;d]` {km³/s²; km³/s²; km} contains the
gravitational parameters of the first and second primary bodies as well as the distance
between them. `t` is time {s}.
"""
function CR3BPdynamics(rv,p::Array,t) #Three body dynamics in Earth/Moon System
    x,y,z,ẋ,ẏ,ż = rv
    μ₁,μ₂,d = p # parameters
    d₁, d₂ = computed1d2(p) # distances from primaries to barycenter
    ωₛ = sqrt((μ₁ + μ₂)/d^3) #rotation rate of system
    r₁³= ((x+d₁)^2 + y^2 + z^2)^1.5 # distance to m1, LARGER MASS
    r₂³= ((x-d₂)^2 + y^2 + z^2)^1.5 # distance to m2, smaller mass

    rvdot = zeros(6)
    rvdot[1:3] = [ẋ;ẏ;ż]
    rvdot[4]   = -(μ₁*(x+d₁)/r₁³) - (μ₂*(x-d₂)/r₂³) + 2*ωₛ*ẏ + ωₛ^2*x
    rvdot[5]   = -(μ₁*y     /r₁³) - (μ₂*y     /r₂³) - 2*ωₛ*ẋ + ωₛ^2*y
    rvdot[6]   = -(μ₁*z     /r₁³) - (μ₂*z     /r₂³)
    return rvdot
end

"""
    CR3BPdynamics!(rvdot, rv, μ, t)

In-place version of `CR3BPdynamics(rv, μ, t)`.
"""
function CR3BPdynamics!(rvdot,rv,μ,t) #Three body dynamics in Earth/Moon System
    rvdot[:] = CR3BPdynamics(rv,μ,t)
    return nothing
end

"""
    CR3BPdynamics!(rvdot, rv, sys::System, t)

In-place version of `CR3BPdynamics(rv, sys::System, t)`.
"""
function CR3BPdynamics!(rvdot,rv,sys::System,t) #Three body dynamics in Earth/Moon System
    rvdot[:] = CR3BPdynamics(rv,sys,t)
    return nothing
end

"""
    CR3BPdynamics!(rvdot, rv, p::Array, t)

In-place version of `CR3BPdynamics(rv, p::Array, t)`.
"""
function CR3BPdynamics!(rvdot,rv,p::Array,t) #Three body dynamics in Earth/Moon System
    rvdot[:] = CR3BPdynamics(rv,p,t)
    return nothing
end

"""
    CR3BPjac(rv, μ)

Compute Jacobian of time derivative w.r.t. state vector [6x6]
"""
function CR3BPjac(rv, μ)
    x,y,z = rv[1:3]

    r₁³= ((x+μ)^2   + y^2 + z^2)^1.5 # distance to m1, LARGER MASS
    r₂³= ((x-1+μ)^2 + y^2 + z^2)^1.5 # distance to m2, smaller mass
    r₁⁵= ((x+μ)^2   + y^2 + z^2)^2.5 # distance to m1, LARGER MASS
    r₂⁵= ((x-1+μ)^2 + y^2 + z^2)^2.5 # distance to m2, smaller mass

    Ωxx = 1 - (1-μ)*(1/r₁³ - 3*(x+μ)^2/r₁⁵) - μ*(1/r₂³ - 3*(x-1+μ)^2/r₂⁵)
    Ωyy = 1 - (1-μ)*(1/r₁³ - 3*    y^2/r₁⁵) - μ*(1/r₂³ - 3*      y^2/r₂⁵)
    Ωzz =   - (1-μ)*(1/r₁³ - 3*    z^2/r₁⁵) - μ*(1/r₂³ - 3*      z^2/r₂⁵)
    
    Ωxy = 3*(1-μ)*(x+μ)*y/r₁⁵ + 3*μ*(x-1+μ)*y/r₂⁵
    Ωxz = 3*(1-μ)*(x+μ)*z/r₁⁵ + 3*μ*(x-1+μ)*z/r₂⁵
    Ωyz = 3*(1-μ)*    y*z/r₁⁵ + 3*μ*      y*z/r₂⁵


    F = [   0     0     0     1     0	 0 ;
            0     0     0     0 	1 	 0 ;
            0	  0     0     0     0    1 ;
          Ωxx   Ωxy   Ωxz     0     2 	 0 ;
          Ωxy   Ωyy   Ωyz    -2     0 	 0 ;
          Ωxz   Ωyz   Ωzz     0	    0	 0 ]
end

"""
    CR3BPjac(rv, sys)

Compute Jacobian of time derivative w.r.t. state vector [6x6]
"""
CR3BPjac(rv, sys::System) = CR3BPjac(rv, sys.μ)

"""
    CR3BPstm(w, μ, t)

Compute time derivative of state vector `w = [r; v; vec(Φ)]` {NON; NON; NON} in the rotating
frame of the normalized CR3BP. `vec(Φ)` is the vectorized state transition matrix while `μ`
is the CR3BP mass parameter μ₂/(μ₁+μ₂) {NON} and `t` is time {NON}.
"""
function CR3BPstm(w,μ,t) #Three body dynamics in Earth/Moon System
    rv = w[1:6]
    Φ = reshape(w[7:42],6,6)
    # x,y,z,ẋ,ẏ,ż = rv

    # r₁³= ((x + μ)^2     + y^2 + z^2)^1.5 # distance to m1, LARGER MASS
    # r₂³= ((x - 1 + μ)^2 + y^2 + z^2)^1.5 # distance to m2, smaller mass
    # r₁⁵= ((x + μ)^2     + y^2 + z^2)^2.5 # distance to m1, LARGER MASS
    # r₂⁵= ((x - 1 + μ)^2 + y^2 + z^2)^2.5 # distance to m2, smaller mass

    # Ωxx = 1 - (1-μ)*(1/r₁³ - 3*(x + μ)^2/r₁⁵) - μ*(1/r₂³ - 3*(x - 1 + μ)^2/r₂⁵)
    # Ωxy = 3*(1-μ)*(x + μ)*y/r₁⁵ + 3*μ*(x - 1 + μ)*y/r₂⁵
    # Ωxz = 3*(1-μ)*(x + μ)*z/r₁⁵ + 3*μ*(x - 1 + μ)*z/r₂⁵
    # Ωyy = 1 - (1-μ)*(1/r₁³ - 3*y^2/r₁⁵) - μ*(1/r₂³ - 3*y^2/r₂⁵)
    # Ωyz = 3*(1-μ)*y*z/r₁⁵ + 3*μ*y*z/r₂⁵
    # Ωzz = - (1-μ)*(1/r₁³ - 3*z^2/r₁⁵) - μ*(1/r₂³ - 3*z^2/r₂⁵)


    # F = [   0     0     0     1     0	 0 ;
    #         0     0     0     0 	1 	 0 ;
    #         0	  0     0     0     0    1 ;
    #       Ωxx   Ωxy   Ωxz     0     2 	 0 ;
    #       Ωxy   Ωyy   Ωyz    -2     0 	 0 ;
    #       Ωxz   Ωyz   Ωzz     0	    0	 0 ]

    F = CR3BPjac(rv, μ)

    Φdot = F*Φ
    wdot = zeros(42)
    wdot[1:6] = CR3BPdynamics(rv,μ,t)
    wdot[7:42] = reshape(Φdot, 36, 1)
    return wdot
end

"""
    CR3BPstm(w, sys, t)

Compute time derivative of state vector `w = [r; v; vec(Φ)]` {NON; NON; NON} in the rotating
frame of the normalized CR3BP. `vec(Φ)` is the vectorized state transition matrix while
`sys` is the CR3BP system and `t` is time {NON}.
"""
function CR3BPstm(w,sys::System,t)
    return CR3BPstm(w,sys.μ,t)
end

"""
    CR3BPstm!(wdot, w, μ, t)

In-place version of `CR3BPstm(w, μ, t)`.
"""
function CR3BPstm!(wdot,w,μ,t)
    wdot[:] = CR3BPstm(w,μ,t)
    return nothing
end

"""
    CR3BPstm!(wdot, w, sys::System, t)

In-place version of `CR3BPstm(w, sys::System, t)`.
"""
function CR3BPstm!(wdot,w,sys::System,t) #Three body dynamics in Earth/Moon System
    wdot[:] = CR3BPstm(w,sys,t)
    return nothing
end

"""
    CR3BPinert(rv,μ,t)

Compute time derivative of state vector `rv = [r; v]` {NON, NON} in the inertial frame of
the normalized CR3BP where `μ` is the CR3BP mass parameter μ₂/(μ₁+μ₂) {NON} and `t` is time
{NON}.
"""
function CR3BPinert(rv,μ,t)
    x,y,z,ẋ,ẏ,ż = rv
    r₁ = [x +       μ*cos(t); y +       μ*sin(t); 0]
    r₂ = [x - (1 - μ)*cos(t); y - (1 - μ)*sin(t); 0]
    r₁³ = norm(r₁)^3
    r₂³ = norm(r₂)^3
    rvdot = zeros(6)
    rvdot[1:3] = [ẋ;ẏ;ż]
    rvdot[4:6] = -(1 - μ)*r₁/r₁³ - μ*r₂/r₂³
    return rvdot
end

"""
    CR3BPinert(rv, sys::System, t)

Compute time derivative of state vector `rv = [r; v]` {NON, NON} in the inertial frame of
the normalized CR3BP where `sys` is the CR3BP system and `t` is time {NON}.
"""
function CR3BPinert(rv,sys::System,t)
    return CR3BPinert(rv,sys.μ,t)
end

"""
    CR3BPinert(rvdot, rv, p::Array, t)

Compute time derivative of state vector `rv = [r; v]` {NON, NON} in the inertial frame of
the non-normalized CR3BP where `p = [μ₁;μ₂;d]` {km³/s²; km³/s²; km} contains the
gravitational parameters of the first and second primary bodies as well as the distance
between them. `t` is time {s}.
"""
function CR3BPinert(rv,p::Array,t)
    x,y,z,ẋ,ẏ,ż = rv
    μ₁,μ₂,d = p # parameters
    d₁,d₂ = computed1d2(p)
    ωₛ = sqrt((μ₁ + μ₂)/d^3)
    r₁ = [x + d₁*cos(ωₛ*t); y + d₁*sin(ωₛ*t); 0] # distance to m1, LARGER MASS
    r₂ = [x - d₂*cos(ωₛ*t); y - d₂*sin(ωₛ*t); 0] # distance to m2, smaller mass
    r₁³= norm(r₁)^3
    r₂³= norm(r₂)^3
    rvdot = zeros(6)
    rvdot[1:3] = rv[4:6]
    rvdot[4:6] = -μ₁*r₁/r₁³ - μ₂*r₂/r₂³
    return rvdot
end

"""
    CR3BPinert!(rvdot, rv, μ, t)

In-place version of `CR3BPinert(rv, μ, t)`.
"""
function CR3BPinert!(rvdot,rv,μ,t)
    rvdot[:] = CR3BPinert(rv,μ,t)
    return nothing
end

"""
    CR3BPinert!(rvdot, rv, sys::System, t)

In-place version of `CR3BPinert(rv, sys::System, t)`.
"""
function CR3BPinert!(rvdot,rv,sys::System,t)
    rvdot[:] = CR3BPinert(rv,sys,t)
    return nothing
end

"""
    CR3BPinert!(rvdot, rv, p::Array, t)

In-place version of `CR3BPinert(rv, p::Array, t)`.
"""
function CR3BPinert!(rvdot,rv,p::Array,t)
    rvdot[:] = CR3BPinert(rv,p,t)
    return nothing
end

"""
    CWdynamics(rv, n, t)

Clohessy-Wiltshire equations. Compute time derivative of state vector `rv = [δr; δv]`
{km; km/s} where `n` {rad/s} is the mean motion of the chief and `t` is time {s}.
"""
function CWdynamics(rv,n,t)
    x,y,z,ẋ,ẏ,ż = rv
    rvdot = zeros(6)
    rvdot[1:3] = rv[4:6]
    rvdot[4] = 2*n*ẏ + 3*n^2*x
    rvdot[5] = -2*n*ẋ
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
See G. Gómez, C. Simó, J. Llibre, and R. Martínez, Dynamics and mission design near
libration points. Vol. II, vol. 3. 2001.

    BCPdynamics(rv, μ, m₃, n₃, t)

Compute time derivative of state vector `rv = [r; v]` {km; km/s} in the normalized
Bicircular Four-Body Problem (BCP). `μ` {NON} is the BCP mass parameter and `m₃` {NON} and
`n₃` {NON} are the normalized mass and mean motion of the tertiary body. `t` is time {NON}.
"""
function BCPdynamics(rv, μ, m₃, n₃, t)
    x,y,z,ẋ,ẏ,ż = rv

    # Xₑ = μ*cos(t);      Yₑ = μ*sin(t);      Zₑ = 0
    # Xₘ = (μ-1)*cos(t);  Yₘ = (μ-1)*sin(t);  Zₘ = 0
    # Xₛ = a₃*cos(n₃*t);  Yₛ = a₃*sin(n₃*t);  Zₛ = 0

    a₃ = (1+m₃)^(1/3)/n₃^(2/3)
    θ = (n₃-1)*t # This is changed from 1-n₃ so now θ will be negative for the Sun and positive for the Moon (BCP2)
    x₃ =  a₃*cos(θ)
    y₃ =  a₃*sin(θ) # I had this as negative before, as shown in the Gomez book

    r₁³ = (  (x+μ)^2 +      y^2 + z^2)^1.5 # distance to m1, LARGER MASS
    r₂³ = ((x-1+μ)^2 +      y^2 + z^2)^1.5 # distance to m2, smaller mass
    r₃³ = ( (x-x₃)^2 + (y-y₃)^2 + z^2)^1.5 # distance to m3, the tertiary body

    rvdot = zeros(6)
    rvdot[1:3] = [ẋ;ẏ;ż]
    rvdot[4] = -((1-μ)*(x+μ)/r₁³ + μ*(x-1+μ)/r₂³ + m₃*(x-x₃)/r₃³ + m₃*cos(θ)/a₃^2) + x + 2*ẏ
    rvdot[5] = -((1-μ)  *  y/r₁³ + μ   *   y/r₂³ + m₃*(y-y₃)/r₃³ + m₃*sin(θ)/a₃^2) + y - 2*ẋ
    rvdot[6] = -((1-μ)  *  z/r₁³ + μ   *   z/r₂³ + m₃  *   z/r₃³)
return rvdot
end

"""
    BCPdynamics(rv, sys::BicircularSystem, t)

Compute time derivative of state vector `rv = [r; v]` {km; km/s} in the normalized
Bicircular Four-Body Problem (BCP). `sys` is the BCP system and `t` is time {NON}.
"""
function BCPdynamics(rv, sys::BicircularSystem, t)
    return BCPdynamics(rv, sys.μ, sys.m₃, sys.n₃, t)
end

"""
    BCPdynamics!(rvdot, rv, μ, m₃, n₃, t)

In-place version of `BCPdynamics(rv, μ, m₃, n₃, t)`.
"""
function BCPdynamics!(rvdot, rv, μ, m₃, n₃, t)
    rvdot[:] = BCPdynamics(rv, μ, m₃, n₃, t)
    return nothing
end

"""
    BCPdynamics!(rvdot, rv, sys::BicircularSystem, t)

In-place version of `BCPdynamics(rv, sys::BicircularSystem, t)`.
"""
function BCPdynamics!(rvdot, rv, sys::BicircularSystem, t)
    rvdot[:] = BCPdynamics(rv, sys, t)
    return nothing
end



"""

    BCPdynamics2(rv, μ, m₃, n₃, t)

Compute time derivative of state vector `rv = [r; v]` {km; km/s} in the normalized
Bicircular Four-Body Problem (BCP). `μ` {NON} is the BCP mass parameter and `m₃` {NON} and
`n₃` {NON} are the normalized mass and mean motion of the tertiary body. `t` is time {NON}.
in BCPdynamics2, the tertiary is assumed to orbit the secondary rather than the barycenter.
For example, if working in the Sun/Earth barycenter with the Moon orbiting the Earth.
"""
function BCPdynamics2(rv, μ, m₃, n₃, t)
    x,y,z,ẋ,ẏ,ż = rv

    # Xₑ = μ*cos(t);      Yₑ = μ*sin(t);      Zₑ = 0
    # Xₘ = (μ-1)*cos(t);  Yₘ = (μ-1)*sin(t);  Zₘ = 0
    # Xₛ = a₃*cos(n₃*t);  Yₛ = a₃*sin(n₃*t);  Zₛ = 0
    
    a₃ = (1+m₃)^(1/3)/n₃^(2/3)
    θ = (n₃-1)*t # This is changed from 1-n₃ so now θ will be negative for the Sun and positive for the Moon (BCP2)
    x₃ = 1 - μ + a₃*cos(θ)
    y₃ = a₃*sin(θ) # I had this as negative before, as shown in the Gomez book

    r₁³ = (  (x+μ)^2 +      y^2 + z^2)^1.5 # distance to m1, LARGER MASS
    r₂³ = ((x-1+μ)^2 +      y^2 + z^2)^1.5 # distance to m2, smaller mass
    r₃³ = ( (x-x₃)^2 + (y-y₃)^2 + z^2)^1.5 # distance to m3, the tertiary body

    rvdot = zeros(6)
    rvdot[1:3] = [ẋ;ẏ;ż]
    rvdot[4] = -((1-μ)*(x+μ)/r₁³ + μ*(x-1+μ)/r₂³ + m₃*(x-x₃)/r₃³) + x + 2*ẏ
    rvdot[5] = -((1-μ)  *  y/r₁³ + μ   *   y/r₂³ + m₃*(y-y₃)/r₃³) + y - 2*ẋ
    rvdot[6] = -((1-μ)  *  z/r₁³ + μ   *   z/r₂³ + m₃  *   z/r₃³)
    return rvdot
end

"""
    BCPdynamics(rv, sys::BicircularSystem, t)

Compute time derivative of state vector `rv = [r; v]` {km; km/s} in the normalized
Bicircular Four-Body Problem (BCP). `sys` is the BCP system and `t` is time {NON}.
"""
function BCPdynamics2(rv, sys::BicircularSystem, t)
    return BCPdynamics2(rv, sys.μ, sys.m₃, sys.n₃, t)
end

"""
    BCPdynamics!(rvdot, rv, μ, m₃, n₃, t)

In-place version of `BCPdynamics(rv, μ, m₃, n₃, t)`.
"""
function BCPdynamics2!(rvdot, rv, μ, m₃, n₃, t)
    rvdot[:] = BCPdynamics2(rv, μ, m₃, n₃, t)
    return nothing
end

"""
    BCPdynamics!(rvdot, rv, sys::BicircularSystem, t)

In-place version of `BCPdynamics(rv, sys::BicircularSystem, t)`.
"""
function BCPdynamics2!(rvdot, rv, sys::BicircularSystem, t)
    rvdot[:] = BCPdynamics2(rv, sys, t)
    return nothing
end


"""
    BCPstm(w, μ, m₃, n₃, t)

Compute time derivative of state vector `w = [r; v; vec(Φ)]` {NON; NON; NON} in the
normalized Bicircular Four-Body Problem (BCP). `vec(Φ)` is the vectorized state transition
matrix. `μ` {NON} is the BCP mass parameter and `m₃` {NON} and `n₃` {NON} are the normalized
mass and mean motion of the tertiary body. `t` is time {NON}.
"""
function BCPstm(w, μ, m₃, n₃, t) #Three body dynamics in Earth/Moon System
    rv = w[1:6]
    Φ = reshape(w[7:42],6,6)
    x,y,z,ẋ,ẏ,ż = rv

    a₃ = (1+m₃)^(1/3)/n₃^(2/3)
    θ = (n₃-1)*t # This is changed from 1-n₃ so now θ will be negative for the Sun and positive for the Moon (BCP2)
    x₃ =  a₃*cos(θ)
    y₃ =  a₃*sin(θ) # I had this as negative before, as shown in the Gomez book
    
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
    wdot = zeros(42)
    wdot[1:6] = BCPdynamics(rv,μ,m₃,n₃,t)
    wdot[7:42] = reshape(Φdot, 36, 1)
    return wdot
end

"""
    BCPstm(w, μ, m₃, n₃, t)

Compute time derivative of state vector `w = [r; v; vec(Φ)]` {NON; NON; NON} in the
normalized Bicircular Four-Body Problem (BCP). `vec(Φ)` is the vectorized state transition
matrix, `sys` is the BCP system and `t` is time {NON}.
"""
function BCPstm(w,sys::BicircularSystem,t) #Three body dynamics in Earth/Moon System
    return BCPstm(w, sys.mu, sys.m3, sys.n3 ,t)
end

"""
    BCPstm!(wdot, w, μ, m₃, n₃, t)
    
In-place version of `BCPstm(w, μ, m₃, n₃, t)`.
"""
function BCPstm!(wdot,w,μ,m₃,n₃,t) #Three body dynamics in Earth/Moon System
    wdot[:] = BCPstm(w,μ,m₃,n₃,t)
    return nothing
end

"""
    BCPstm!(wdot, w, sys::BicircularSystem, t)
    
In-place version of `BCPstm(w, sys::BicircularSystem, t)`.
"""
function BCPstm!(wdot,w,sys::BicircularSystem,t) #Three body dynamics in Earth/Moon System
    wdot[:] = BCPstm(w,sys,t)
    return nothing
end

#TODO: Should I make t=0 the default for all these dynamics functions?