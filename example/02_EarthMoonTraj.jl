@time using ThreeBodyProblem
using OrdinaryDiffEq
using Plots
pyplot() # This command sets pyplot as our plotting backend

μ₁ = 398600 # {km³/s²} gravitational parameter of Earth
μ₂ = 4902   # {km³/s²} gravitational parameter of the Moon
d = 384400  # {km} average distance between Earth and the Moon

p = [μ₁, μ₂, d] # parameters that define a CR3BP system

L1 = computeL1(p) # calculate the Lagrange point L1 for the given system parameters p

Lpts = computeLpts(p)
L4 = Lpts[4]

R₁,R₂ = computed1d2(p) # {km} distances of Primary and Secondary bodies from origin
ωₛ = sqrt((μ₁ + μ₂)/d^3) # {rad/s} rotation rate of system


N = 100
X = range(-1.5*d,1.5*d,length=N)
Y = range(-d,d,length=N)

f(x,y) = begin
    -(x^2 + y^2)*ωₛ^2/2 - μ₁/sqrt((x+R₁)^2 + y^2) - μ₂/sqrt((x-R₂)^2 + y^2)
    # rv = [x, y, 0, 0, 0, 0]
    # computeUeff(rv,p)
end
contour(X,Y,f,levels=200,fill=true)

Rₑ = 6378.0 # {km} radius of the Earth
Rₘ = 1738.0 # {km} radius of the Moon

h = 200.0           # {km} altitude of parking orbit
vᵢ = 10.92367104    # {km/s} injection velocity in rotating frame
ϕ = 47.70061087     # {°} injection angle, measured from +y


r₀ = [R₁ - (Rₑ + 200)*cosd(ϕ); -(Rₑ + 200)*sind(ϕ); 0]
v₀ = vᵢ*[sind(ϕ); -cosd(ϕ); 0];
rv₀ = [r₀;v₀] # {km; km/s} our initial state

tspan = (0.,86400*6.) # {sec} 6 day span


prob = ODEProblem(CR3BPdynamics!,rv₀,tspan,p) # CR3BPdynamics! is our in-place dynamics function
sol = solve(prob,Tsit5(),reltol=1e-6)

plot(sol,vars=(1,2),title="Lunar return trajectory in rotating frame")
plot!(circle(Rₑ,[R₁;0]),color="blue",label="Earth")
plot!(circle(Rₘ,[R₂;0]),color="gray",label="Moon",aspect_ratio=:equal)

rv_inert = [rot2inert(sol.u[i],sol.t[i],p) for i = 1:length(sol)]
x = [rv_inert[i][1] for i = 1:length(rv_inert)]
y = [rv_inert[i][2] for i = 1:length(rv_inert)]
plot(x,y,legend=:topleft)
