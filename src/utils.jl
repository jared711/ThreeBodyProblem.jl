"""
xdist(μ₁,μ₂,d)
"""
function xdist(μ₁,μ₂,d)
    R1 = d*mu₂/(mu₁+mu₂)
    R2 = d*mu₁/(mu₁+mu₂)
    return R₁, R₂
end

"""
    findL1(p;tol=1e-15)
"""
function findL1(p::Array;tol=1e-15)
    μ₁,μ₂,d = p
    μ = μ₂/(μ₁ + μ₂)
    L1 = findL1(μ;tol=tol)*d
end

"""
    different method
"""
function findL1(μ;tol=1e-15)
    α = (μ/3*(1 - μ))^(1/3)
    dα = 1
    while abs(dα) > tol
        α₀ = α
        α = (μ*(1 - α)^2 / (3 - 2*μ - α*(3 - μ - α)))^(1/3);
        dα = α - α₀
    end
    L1 = [1-μ-α;0;0]
end
