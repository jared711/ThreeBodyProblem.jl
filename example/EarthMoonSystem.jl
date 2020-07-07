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

N = 100
X = range(-1.5*d,1.5*d,length=N)
Y = range(-d,d,length=N)

Ueff = zeros(N,N)
for i = 1:N
    for j = 1:N
        rv = [X[i]; Y[j]; 0; 0; 0; 0]
        Ueff[i,j] = findUeff(rv, p)
    end
end

#Plot Ueff

Rₑ = 6378.0 # {km} radius of the Earth
Rₘ = 1738.0 # {km} radius of the Moon
h = 200.0   # {km} altitude of parking orbit

vᵢ = 10.92367104    # {km/s} injection velocity in rotating frame
ϕ = 47.70061087     # {°}

R₁,R₂ = findR1R2(p)
r₀ = [-R₁ - (Rₑ + 200)*cosd(ϕ); -(Rₑ + 200)*sind(ϕ); 0]
v₀ = vᵢ*[sind(ϕ); -cosd(ϕ); 0];
rv₀ = [r₀;v₀]

tspan = (0.,86400*6.) # {sec} 6 day span

prob = ODEProblem(CR3BPdynamics!,rv₀,tspan,p)
sol = solve(prob,reltol=1e-6)

plot(sol,vars=(1,2),title="Lunar return trajectory in rotating frame",)
