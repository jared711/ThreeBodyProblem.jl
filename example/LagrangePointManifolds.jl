using ThreeBodyProblem
using DifferentialEquations
using Plots
using LinearAlgebra
# using ForwardDiff

# The first step is to define the system. ThreeBodyProblem.jl has built-in functions for common systems like Sun/Jupiter.
sys = sun_jupiter()

## There's even a recipe to plot the system for us!
plot(sys)

## We can see the Lagrange Points, but where are the Sun and Jupiter? Turns out that the distance between them is so much larger than the radius of either one that it makes it hard to see them on a truly scaled picture. Let's scale the Sun and Jupiter to make them visible.
plot(sys, scaled = true)

## That's better! Now let's compute the L1 and L2 Lagrange points of our system

L1, L2 = computeLpts(sys)
# Declare state vectors for L1 and L2 (with zero velocity)
rv1 = [L1; zeros(3)] # state at Lagrange point L1 (zero velocity added on)
rv2 = [L2; zeros(3)] # state at Lagrange point L2

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
plot!(aspect_ratio=1,ylims=[-0.04,0.04],xlims=[0.93,1.07],legend=:outerright,flip=false)
