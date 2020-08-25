"""
#     findR1R2(μ)
# """
function findR1R2(μ::Number)
    R₁ = -μ
    R₂ = 1-μ
    return R₁, R₂
end

"""
    findR1R2(p::Array)
"""
function findR1R2(p::Array)
    μ₁,μ₂,d = p
    R₁ = d*μ₂/(μ₁+μ₂)
    R₂ = d*μ₁/(μ₁+μ₂)
    return R₁, R₂
end

"""
    findr1r2(rv,μ)
"""
function findr1r2(rv,μ)
    x,y,z = rv[1:3]
    r₁ = (x + μ)^2      + y^2 + z^2
    r₂ = (x - 1 + μ)^2  + y^2 + z^2
    return r₁, r₂
end

"""
    findr1r2(rv,μ)
"""
function findr1r2(rv,p::Array)
    x,y,z = rv[1:3]
    R₁,R₂ = findR1R2(p)
    r₁ = (x + R₁)^2 + y^2 + z^2
    r₂ = (x - R₂)^2 + y^2 + z^2
    return r₁, r₂
end


"""
    findL1(μ;tol=1e-15)

Find 3D L1 in a normalized CR3BP given μ, the system mass ratio.
"""
function findL1(μ;tol=1e-15)
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
    findL1(p::Array;tol=1e-15)

Find 3D L1 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function findL1(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L1 = findL1(μ;tol=tol)*d
    return L1
end

"""
    findL2(μ;tol=1e-15)

Find 3D L2 in a normalized CR3BP given μ, the system mass ratio.
"""
function findL2(μ;tol=1e-15)
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
    findL2(p::Array;tol=1e-15)

Find 3D L2 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function findL2(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L2 = findL2(μ;tol=tol)*d
    return L2
end

"""
    findL3(μ;tol=1e-15)

Find 3D L3 in a normalized CR3BP given μ, the system mass ratio.
"""
function findL3(μ;tol=1e-15)
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
    findL3(p::Array;tol=1e-15)

Find 3D L3 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function findL3(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L3 = findL3(μ;tol=tol)*d
    return L3
end

"""
    findL4(μ;tol=1e-15)

Find 3D L4 in a normalized CR3BP given μ, the system mass ratio.
"""
findL4(μ;tol=1e-15) = [0.5-μ; √3/2; 0]

"""
    findL4(p::Array;tol=1e-15)

Find 3D L4 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function findL4(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L4 = findL4(μ;tol=tol)*d
    return L4
end

"""
    findL5(μ;tol=1e-15)

Find 3D L5 in a normalized CR3BP given μ, the system mass ratio.
"""
findL5(μ;tol=1e-15) = [0.5-μ; -√3/2; 0]

"""
    findL5(p::Array;tol=1e-15)

Find 3D L5 in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function findL5(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L5 = findL5(μ;tol=tol)*d
    return L5
end

"""
    findLpts(μ;tol=1e-15)

Find 3D Lagrange points in a normalized CR3BP given μ, the system mass ratio.
"""
function findLpts(μ; tol=1e-15)
    Lpts = [
        findL1(μ,tol=tol),
        findL2(μ,tol=tol),
        findL3(μ,tol=tol),
        findL4(μ,tol=tol),
        findL5(μ,tol=tol),
        ]
    return Lpts
end

"""
    findLpts(p::Array;tol=1e-15)

Find 3D Lagrange points in a non-normalized CR3BP given p = [μ₁,μ₂,d], which are the
gravitational parameters of the first and second primary bodies [km³/s²]and the
distance between them [km].
"""
function findLpts(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    Lpts = findLpts(μ, tol=tol)*d
    return Lpts
end

"""
    findUeff(rv,μ)

Find effective potential given normalized state rv {NON} and mass ratio {NON}
"""
function findUeff(rv,μ)
    x,y,z = rv[1:3]
    r₁,r₂ = findr1r2(rv,μ)
    Ueff = -(x^2 + y^2)/2 - (1-μ)/r₁ - μ/r₂;
    return Ueff
end

"""
    findUeff(rv,p::Array)

Find effective potential given state rv = [r; v] {km; km/s} and p = [μ₁,μ₂,d],
which are the gravitational parameters of the first and second primary bodies
[km³/s²] and the distance between them [km].
"""
function findUeff(rv,p::Array)
    x,y,z = rv[1:3]
    R₁,R₂ = findR1R2(p)
    r₁,r₂ = findr1r2(rv,p)
    μ₁,μ₂,d = p
    ωₛ = sqrt((μ₁ + μ₂)/d^3)
    Ueff = -(x^2 + y^2)*ωₛ^2/2 - μ₁/r₁ - μ₂/r₂;
    return Ueff
end

"""
    findC(rv,μ)

Find Jacobi constant given normalized state rv {NON} and mass ratio {NON}
"""
function findC(rv,μ)
    v = norm(rv[4:6])
    Ueff = findUeff(rv,μ)
    C = -2*Ueff - v^2
    return C
end

"""
    findC(rv,p::Array)

Find Jacobi constant given state rv = [r; v] {km; km/s} and p = [μ₁,μ₂,d],
which are the gravitational parameters of the first and second primary bodies
[km³/s²] and the distance between them [km].
"""
function findC(rv,p::Array)
    μ₁,μ₂,d = p
    v = norm(rv[4:6])
    Ueff = findUeff(rv,p)
    C = -2*Ueff - v^2
    return C
end

"""
    stability_index(Φ)

Find the stability index for a trajectory given its state transition matrix Φ
"""
function stability_index(Φ)
    return Φ
end
