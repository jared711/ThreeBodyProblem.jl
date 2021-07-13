using ThreeBodyProblem
using DifferentialEquations
using Plots
using LinearAlgebra
pyplot()

sys = jupiter_europa()
rv0 = [1.000045941880842, -0.000083429669677,0.002323145484474,-0.093200789960739,0.099227116703640,-0.033175567182715]
r0 = rv0[1:3]
v0 = rv0[4:6]

# event function
condition(u, t, integrator) = u[2] + max(u[1],0.)
affect!(integrator) = terminate!(integrator)
cb = DifferentialEquations.ContinuousCallback(condition, affect!)

tspan = (0.,-100.)
prob = ODEProblem(CR3BPdynamics!,rv0,tspan,sys)
sol = solve(prob,TsitPap8(),reltol=1e-12,callback=cb)
plot(sol,vars=(1,2),flip=false,aspect_ratio=:equal,legend=:outerright,label="nominal",line=3)

# Set up initial conditions
# R = 1 #[km]
# R /= sys.RUNIT
# rs,N = deserno_sphere(100)
# rvs = [[R*rs[1:3,i]+r0; v0] for i = 1:N]

# Set up initial conditions
c = [1-sys.μ,0,0]
rs = spherical_ring(c,r0-c,1)
rvs = [[rs[i]; v0] for i = 1:N]

funnel = []
for i = 1:N
    prob = ODEProblem(CR3BPdynamics!,rvs[i],tspan,sys)
    sol = solve(prob,TsitPap8(),reltol=1e-12,callback=cb)
    plot!(sol,vars=(1,2),label="",color=:red)
    push!(funnel,sol)
end
plot!(flip=false,aspect_ratio=:equal,legend=:outerright,xlim=[0.995,1.03],ylim=[-0.01,0.02])
