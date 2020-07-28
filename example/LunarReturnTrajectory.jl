using ThreeBodyProblem
using DifferentialEquations
using Plots


μ₁ = 398600 # {km³/s²} gravitational parameter of Earth
μ₂ = 4902   # {km³/s²} gravitational parameter of the Moon
d = 384400  # {km} average distance between Earth and the Moon
p = [μ₁, μ₂, d];

Lpts = findLpts(p)
L4 = Lpts[4]
L1 = findL1(p)

R₁,R₂ = findR1R2(p)
ωₛ = sqrt((μ₁ + μ₂)/d^3)

N = 100
X = range(0,1.2*d,length=N)
Y = range(-d/2,d/2,length=N)

f(x,y) = begin
    -(x^2 + y^2)*ωₛ^2/2 - μ₁/sqrt((x+R₁)^2 + y^2) - μ₂/sqrt((x-R₂)^2 + y^2)
    # rv = [x; y; 0; 0; 0; 0]
    # findUeff(rv,p)
end

contour(X,Y,f,levels=200,fill=true)

Rₑ = 6378.0 # {km} radius of the Earth
Rₘ = 1738.0 # {km} radius of the Moon
h = 200.0   # {km} altitude of parking orbit

vᵢ = 10.92367104    # {km/s} injection velocity in rotating frame
ϕ = 47.70061087     # {°}

r₀ = [-R₁ - (Rₑ + 200)*cosd(ϕ); -(Rₑ + 200)*sind(ϕ); 0]
v₀ = vᵢ*[sind(ϕ); -cosd(ϕ); 0];
rv₀ = [r₀;v₀]

tspan = (0.,86400*6.) # {sec} 6 day span

prob = ODEProblem(CR3BPdynamics!,rv₀,tspan,p)
sol = solve(prob,reltol=1e-6)

plot(sol,vars=(1,2),title="Lunar return trajectory in rotating frame",label="")
plot_circle(Rₑ,[-R₁;0],color="blue",label="Earth")
plot_circle(Rₘ,[R₂;0],color="gray",label="Moon")

for i = 1:length(sol)
    S2I!(sol.u[i],sol.t[i],p)
end

plot(sol,vars=(1,2),title="Lunar return trajectory in inertial frame",label="")
# plot_circle(Rₑ,[-R₁*cos(ωₛ*tₚ);],color="blue",label="Earth")
# plot_circle(Rₘ,[R₂;0],color="gray",label="Moon")

S2I!(rv₀,0,p)
prob = ODEProblem(CR3BPinert!,rv₀,tspan,p)
sol = solve(prob,reltol=1e-6)
plot(sol,vars=(1,2),title="Lunar return trajectory in inertial frame",label="")
# plot_circle(Rₑ,[-R₁*cos(ωₛ*tₚ);],color="blue",label="Earth")
# plot_circle(Rₘ,[R₂;0],color="gray",label="Moon")
