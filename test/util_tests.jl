using ThreeBodyProblem
using Test

μ₁ = 398600     # Earth gravitational parameter
μ₂ = 4902       # Moon gravitational paramter
d = 384400      # Avg distance between Earth and Moon
p = [μ₁; μ₂; d]
μ = μ₂/(μ₁ + μ₂)
@test findR1R2(μ) == (-μ, 1-μ)
@test findL1(μ) == [0.8369247061700067;0;0]
@test findL1(p) == [321713.8570517506;0;0]
@test findL2(μ) == [1.1556746768583537;0;0]
@test findL2(p) == [444241.3457843512;0;0]
@test findL3(μ) == [-1.005061834632099;0;0]
@test findL3(p) == [-386345.76923257887;0;0]
@test findL4(μ) == [0.48785136133154233;0.8660254037844386;0]
@test findL4(p) == [187530.06329584488;332900.16521473817;0]
@test findL5(μ) == [0.48785136133154233;-0.8660254037844386;0]
@test findL5(p) == [187530.06329584488;-332900.16521473817;0]

@test findLpts(μ,tol=1e-3) != findLpts(μ)
@test findLpts(p,tol=1e-3) != findLpts(p)

rv = [1,0,0,0,1,0]
@test findUeff(rv,p) < 0
@test findUeff([-μ;0;0],μ) == -Inf

@test_broken findC(rv,p) > 0
@test_broken findC([-μ;0;0],μ) == Inf

Φ = ones(6,6)
@test_broken stability_index(ϕ) > 1
