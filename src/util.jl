"""
     computeR1R2(μ)
"""
function computeR1R2(μ::Number)
    R₁ = -μ
    R₂ = 1-μ
    return R₁, R₂
end

"""
    computeR1R2(p::Array)
"""
function computeR1R2(p::Array)
    μ₁,μ₂,d = p
    R₁ = d*μ₂/(μ₁+μ₂)
    R₂ = d*μ₁/(μ₁+μ₂)
    return R₁, R₂
end

"""
    computer1r2(rv,μ)
"""
function computer1r2(rv,μ)
    x,y,z = rv[1:3]
    r₁ = (x + μ)^2      + y^2 + z^2
    r₂ = (x - 1 + μ)^2  + y^2 + z^2
    return r₁, r₂
end

"""
    computer1r2(rv,μ)
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

Compute 3D L1 in a normalized CR3BP given μ, the system mass ratio.
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

Compute 3D L2 in a normalized CR3BP given μ, the system mass ratio.
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

Compute 3D L3 in a normalized CR3BP given μ, the system mass ratio.
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

ComputeL4 3D L4 in a normalized CR3BP given μ, the system mass ratio.
"""
computeL4(μ;tol=1e-15) = [0.5-μ; √3/2; 0]

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

Compute 3D L5 in a normalized CR3BP given μ, the system mass ratio.
"""
computeL5(μ;tol=1e-15) = [0.5-μ; -√3/2; 0]

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

Compute 3D Lagrange points in a normalized CR3BP given μ, the system mass ratio.
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
    computeLpts(sys::System;tol=1e-15)

Compute 3D Lagrange points in a non-normalized CR3BP given system S.
"""
function computeLpts(sys::System;tol=1e-15)
    μ = sys.μ
    Lpts = computeLpts(μ,tol=tol)
    return Lpts
end

"""
    computeUeff(rv,μ)

Compute effective potential given normalized state rv {NON} and mass ratio {NON}
"""
function computeUeff(rv,μ)
    x,y,z = rv[1:3]
    r₁,r₂ = computer1r2(rv,μ)
    Ueff = -(x^2 + y^2)/2 - (1-μ)/r₁ - μ/r₂;
    return Ueff
end

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

Compute Jacobi constant given normalized state rv {NON} and mass ratio {NON}
"""
function computeC(rv,μ)
    v = norm(rv[4:6])
    Ueff = computeUeff(rv,μ)
    C = -2*Ueff - v^2
    return C
end

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
    computeC(rv,sys::System)

Compute Jacobi constant given normalized state rv {NON} and system
"""
function computeC(rv,sys::System)
    C = computeC(rv, sys.μ)
    return C
end

"""
    stability_index(Φ)

Compute the stability index for a trajectory given its state transition matrix Φ
"""
function stability_index(Φ)
    return Φ
end
