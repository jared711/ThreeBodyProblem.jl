using ThreeBodyProblem
using DifferentialEquations
using Plots
using LinearAlgebra

# Our goal today is to compute a halo orbit in the Sun Earth system

## The first step is to define the system. ThreeBodyProblem.jl has built-in functions for common systems like Sun/Jupiter.
sys = sun_earth()
sys.name

## If you want to create a custom system, then you would use the following syntax
sys = System(ThreeBodyProblem.SUN, ThreeBodyProblem.EARTH)

## Now, we are going to compute our first guess at a Halo Orbit using Richardson's expansion
### The details are complicated, but if you're curious check out this paper
#### D. L. Richardson, “Analytic Construction Of Periodic Orbits About The Collinear Points,” Celest. Mech., vol. 22, no. 3, pp. 241–253, 1980, doi: 10.1007/BF01229511.
Az = 0.001
Lpt = 2 # Which libration point will I center about
NS = 1
npts = 100
t, rvs, T, Ax = rich3(sys, Az, Lpt, NS, npts)

## Let's plot the orbit to see what it looks like
plot(rvs[:,1],rvs[:,2],rvs[:,3], label="Richardson")

## Remember, this is a third order approximation
## What will happen if we actually integrate this trajectory?
tspan = (0., T)
prob = ODEProblem(CR3BPdynamics!,rvs[1,:],tspan,sys)
sol = solve(prob, reltol=1e-12)
plot!(sol,vars=(1,2,3),label="Actual",linecolor=:red)

## We can see that this orbit diverges before completing a period
## We need to use a differential corrector to hone in on the true periodic orbit
rv₀, ttf = differential_corrector(sys, rvs[1,:], myconst=3, tf=T)


Φ₀ = I(6)

# event function
condition(u, t, integrator) = u[2]
affect!(integrator) = terminate!(integrator)
cb = DifferentialEquations.ContinuousCallback(condition, affect!)
# sol = solve(prob, Tsit5(), calback=cb)

w₀ = vcat(rv₀, reshape(Φ₀,36,1))

prob = ODEProblem(CR3BPstm!,w₀,(0.0,20.0),sys)
tol = 1e-12
sol = solve(prob, reltol=tol, calback=cb)
plot(sol,vars=(1,2),xlims=[0.98,1.02],ylims=[-0.02,0.02])

w = sol[end]
rv = w[1:6]
Φ = reshape(w[7:42],6,6)
global T = 2*sol.t[end]




## Declare state vectors for L1 and L2 (with zero velocity)
L2 = computeL2(sys.μ,tol=1e-15)
rv2 = [L2; zeros(3)] # state at Lagrange point L2

## Declare state vectors for L1 and L2 (with zero velocity)
L1 = computeL1(sys.μ,tol=1e-15)
rv1 = [L1; zeros(3)] # state at Lagrange point L1

## The Lagrange points are equilibrium points in our dynamics. This means that an object placed there perfectly will stay there forever. But, Lagrange points L1, L2, and L3 are unstable, meaning if you perturb the object slightly, it may fall away. This is how we will generate the invariant manifolds.
## We first need to linearize the dynamics about these Lagrange points to determine which directions are stable and unstable

Φ₀ = I(6)
w₀ = [reshape(Φ₀,36,1);rv1]
tspan = (0.,1.)
prob = ODEProblem(CR3BPstm!,w₀,tspan,sys)
sol = solve(prob, reltol=1e-6)
Φₜ = Matrix(reshape(sol.u[end][1:36],6,6))
rvₜ = sol.u[end][37:42] # last time step should match very closely with rv₁
D,V = eigen(Φₜ,sortby=isreal)
Yw = real(V[:,findall(isreal, D)])
D = D[findall(isreal, D)]
Yws = Yw[:,findmin(real(D))[2]]
Ywu = Yw[:,findmax(real(D))[2]]

## We perturb our states in the stable and unstable directions
α = 1e-6

# L1 points
rv1up = rv1 + α*Ywu # unstable manifold + side
rv1un = rv1 - α*Ywu # unstable manifold - side
rv1sp = rv1 + α*Yws # stable manifold + side
rv1sn = rv1 - α*Yws # stable manifold - side




## Next, we integrate our initial conditions forward in time for the +x and -x directions and backward in time for the +y and -y directions
# Integrate trajectories for 2 Jupiter periods
tf = 2*2π
tspan_forward = (0.,tf) # make sure to add decimal points so tspan contains floating point values, not integers
tspan_back = (0.,-tf)
myreltol = 1e-12

# Set up the ODE Problems
prob1up = ODEProblem(CR3BPdynamics!,rv1up,tspan_forward,sys) # Unstable positive
prob1un = ODEProblem(CR3BPdynamics!,rv1un,tspan_forward,sys) # Unstable negative
prob1sp = ODEProblem(CR3BPdynamics!,rv1sp,tspan_back,sys) # Stable positive
prob1sn = ODEProblem(CR3BPdynamics!,rv1sn,tspan_back,sys) # Stable negative

# Solutions to the ODEs
sol1up = solve(prob1up, reltol=myreltol)
sol1un = solve(prob1un, reltol=myreltol)
sol1sp = solve(prob1sp, reltol=myreltol)
sol1sn = solve(prob1sn, reltol=myreltol)

plot(sys,scaled=true)
plot!(sol1up,vars=(1,2),label="Wu+",linecolor=:red)
plot!(sol1un,vars=(1,2),label="Wu-",linecolor=:magenta)
plot!(sol1sp,vars=(1,2),label="Ws+",linecolor=:blue)
plot!(sol1sn,vars=(1,2),label="Ws-",linecolor=:cyan)
plot!(aspect_ratio=1,ylims=[-1,1],xlims=[-0.6,1.3],legend=:topleft,flip=false)

## We can see two strands moving in towards the sun, but the rest is a jumbled mess. Let's zoom in to see what's going on. Also, I need to make sure things are scaled correctly, or Jupiter will look way too big.
plot(sys)
plot!(sol1up,vars=(1,2),label="Wu+",linecolor=:red)
plot!(sol1un,vars=(1,2),label="Wu-",linecolor=:magenta)
plot!(sol1sp,vars=(1,2),label="Ws+",linecolor=:blue)
plot!(sol1sn,vars=(1,2),label="Ws-",linecolor=:cyan)
plot!(aspect_ratio=1,ylims=[-0.04,0.04],xlims=[0.93,1.06],legend=:outerright,flip=false)

# L2 points
rv2up = rv2 + α*Ywu # unstable manifold + side
rv2un = rv2 - α*Ywu # unstable manifold - side
rv2sp = rv2 + α*Yws # stable manifold + side
rv2sn = rv2 - α*Yws # stable manifold - side

## Now we look at the L2 point manifolds
prob2up = ODEProblem(CR3BPdynamics!,rv2up,tspan_forward,sys) # Unstable positive
prob2un = ODEProblem(CR3BPdynamics!,rv2un,tspan_forward,sys) # Unstable negative
prob2sp = ODEProblem(CR3BPdynamics!,rv2sp,tspan_back,sys) # Stable positive
prob2sn = ODEProblem(CR3BPdynamics!,rv2sn,tspan_back,sys) # Stable negative

sol2up = solve(prob2up, reltol=1e-6)
sol2un = solve(prob2un, reltol=1e-6)
sol2sp = solve(prob2sp, reltol=1e-6)
sol2sn = solve(prob2sn, reltol=1e-6)

plot(sys, scaled=true)
plot!(sol2up,vars=(1,2),label="Wu+",linecolor=:red)
plot!(sol2un,vars=(1,2),label="Wu-",linecolor=:magenta)
plot!(sol2sp,vars=(1,2),label="Ws+",linecolor=:blue)
plot!(sol2sn,vars=(1,2),label="Ws-",linecolor=:cyan)
plot!(aspect_ratio=1,ylims=[-1.0,1.0],xlims=[-0.3,1.5],legend=:topright,flip=false)

## This time there are two strands moving out away from the Sun, along with another jumble around Jupiter. Let's zoom in again to take a look.
plot(sys)
plot!(sol2up,vars=(1,2),label="Wu+",linecolor=:red)
plot!(sol2un,vars=(1,2),label="Wu-",linecolor=:magenta)
plot!(sol2sp,vars=(1,2),label="Ws+",linecolor=:blue)
plot!(sol2sn,vars=(1,2),label="Ws-",linecolor=:cyan)
plot!(aspect_ratio=1,ylims=[-0.01,0.01],xlims=[0.975,1.025],legend=:outerright,flip=false)
