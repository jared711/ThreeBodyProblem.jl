using ThreeBodyProblem
using DifferentialEquations
using Plots
gr()

sys = jupiter_europa()

rv₀ = [1.000045941880842, -0.000083429669677, 0.002323145484474, -0.093200789960739, 0.099227116703640, -0.033175567182715]

tspan = (0.,-22.)
prob = ODEProblem(CR3BPdynamics!,rv₀,tspan,sys)
sol = solve(prob,reltol=1e-6)

c = [1-sys.μ, 0, 0]
xyz_hemisphere = deserno_hemisphere(100, c)

plot(sol,vars=(1,2),title="Lunar return trajectory in rotating frame")
plot!(sys)
