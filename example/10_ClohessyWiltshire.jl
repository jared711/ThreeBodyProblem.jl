using ThreeBodyProblem
using DifferentialEquations
using Plots

μ₁ = 398600 # {km³/s²} gravitational parameter of Earth
μ₂ = 4903   # {km³/s²} gravitational parameter of the Moon
d = 384400  # {km} average distance between Earth and the Moon
p = [μ₁, μ₂, d];

R₁,R₂ = computeR1R2(p)
ωₛ = sqrt((μ₁ + μ₂)/d^3);

Rₑ = 6378.0 # {km} radius of the Earth
Rₘ = 1738.0 # {km} radius of the Moon
h = 200.0   # {km} altitude of parking orbit

L1 = computeL1(p)

vᵢ₀ = 10.92367104   # {km/s} synodic frame injection velocity
ϕ₀ = 47.70061087    # {°}

vᵢ = vᵢ₀
ϕ = ϕ₀

r₀ = [-R₁ - (Rₑ + 200)*cosd(ϕ); -(Rₑ + 200)*sind(ϕ); 0]
v₀ = vᵢ*[sind(ϕ); -cosd(ϕ); 0];
rv₀ = [r₀;v₀]

tspan = [0.,86400*3] # 3 days

prob = ODEProblem(CR3BPdynamics!,rv₀,tspan,p)
sol = solve(prob,reltol=1e-6)

plot(sol,vars=(1,2),title="Initial Guess",label="")
# scatter!(L1[1],L1[2],marker="x")
plot!(circle(Rₑ,[-R₁;0]),color="blue",label="Earth")
plot!(circle(Rₘ,[R₂;0]),color="gray",label="Moon",aspect_ratio=1)
