using ThreeBodyProblem
using DifferentialEquations
using Plots
using LinearAlgebra

# Our goal today is to compute a halo orbit in the Sun Earth system

## The first step is to define the system. ThreeBodyProblem.jl has built-in functions for common systems like Sun/Jupiter.
sys = earth_moon()
sys.name

## Use the following syntax
sys = System(ThreeBodyProblem.SUN, ThreeBodyProblem.EARTH)

## Now, we are going to compute our first guess at a Halo Orbit using Richardson's expansion
### The details are complicated, but if you're curious check out this paper
#### D. L. Richardson, “Analytic Construction Of Periodic Orbits About The Collinear Points,” Celest. Mech., vol. 22, no. 3, pp. 241–253, 1980, doi: 10.1007/BF01229511.
Az = 0.001 # The amplitude of the orbit in the z direction
Lpt = 2 # Which libration point will I center about
NS = 1
npts = 100
t, rvs, T, Ax = rich3(sys, Az, Lpt, NS, npts)

## Let's plot the orbit to see what it looks like
plot(rvs[:,1],rvs[:,2],rvs[:,3],label="Richardson Approximation")

## Remember, this is a third order approximation
## What will happen if we actually integrate this trajectory?
tspan = (0., T)
prob = ODEProblem(CR3BPdynamics!,rvs[1,:],tspan,sys)
sol = solve(prob, reltol=1e-12)
plot!(sol,vars=(1,2,3),label="Richardson Integrated",linecolor=:red)

## We can see that this orbit diverges before completing a period
## We need to use a differential corrector to hone in on the true periodic orbit
rv₀, T = differential_corrector(sys, rvs[1,:], myconst=3, tf=T)

plot(sys)
plot!(sol,vars=(1,2),color=:black,label="Actual Halo")
plot!(aspect_ratio=:equal,ylims=[-0.01,0.01],xlims=[0.985,1.015],legend=:outerright,flip=false)

# This function requires the initial condition and period of a periodic orbit
Wsp, Wsn, Wup, Wun = invariant_manifolds(sys,rv₀,T,tf=10.,nPts=100)

for i = 1:length(Wsp)
    plot!(Wsp[i],vars=(1,2),label="",linecolor=:blue)
    plot!(Wsn[i],vars=(1,2),label="",linecolor=:cyan)
    plot!(Wup[i],vars=(1,2),label="",linecolor=:red)
    plot!(Wun[i],vars=(1,2),label="",linecolor=:magenta)
end
plot!(flip=false,aspect_ratio=:equal,legend=:outerright,xlim=[0.995,1.03],ylim=[-0.02,0.02])