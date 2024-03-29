using ThreeBodyProblem
using RecipesBase
using Test
using Plots

### circle, sphere, torus ###
# circle
r = 3
x,y = circle(r)
@test x[1]^2 + y[1]^2 == r^2
# sphere
x,y,z = sphere(r)
@test x[1]^2 + y[1]^2 + z[1]^2 == r^2
# torus
a, b, c = 1, 5, [1,1,1]
x,y,z = torus(a,b,c)
@test x[1] == a + b + c[1]

### plot recipes ###
# system
sys = earth_moon()
p = plot(sys);
@test !isempty(p)
p = plot(sys,scaled=true);
@test !isempty(p)
p = plot(sys,center=[1,0,0]);
@test !isempty(p)
@test_throws ErrorException plot(sys,:planar=>false); # ("Cannot convert System to series data for plotting")
# RecipesBase method
@test !isempty(RecipesBase.apply_recipe(Dict{Symbol, Any}(), sys))
@test_throws ErrorException RecipesBase.apply_recipe(Dict{Symbol, Any}(:planar=>false), sys)
# body
p = plot(sys.prim);
@test !isempty(p)
p = plot(sys.sec, planar=false);
@test !isempty(p)
@test !isempty(RecipesBase.apply_recipe(Dict{Symbol, Any}(), sys.prim))
@test !isempty(RecipesBase.apply_recipe(Dict{Symbol, Any}(:planar=>false), sys.prim))
# trajectory
Lpts = computeLpts(sys)
t = LinRange(0,2π,10)
traj = [rot2inert([Lpts[1];zeros(3)],t[i],sys) for i ∈ eachindex(t)]
traj4D= [traj[i][[1,2,4,5]] for i ∈ eachindex(t)]
p = plot(traj);
@test !isempty(p)
p = plot(traj4D);
@test !isempty(p)
p = plot(traj, vel=true);
@test !isempty(p)
p = plot(traj4D);
@test !isempty(p)
@test !isempty(RecipesBase.apply_recipe(Dict{Symbol, Any}(), traj))
@test !isempty(RecipesBase.apply_recipe(Dict{Symbol, Any}(), traj4D))
@test !isempty(RecipesBase.apply_recipe(Dict{Symbol, Any}(:vel=>true), traj))
@test !isempty(RecipesBase.apply_recipe(Dict{Symbol, Any}(:vel=>true), traj4D))

#Lpts = computeLpts(sys)
#tt = 0:π/10:2π
#xx = [[Lpts[1];zeros(3)] for t in tt]
#xx_inert = [rot2inert(xx[i],tt[i],sys) for i in eachindex(tt)]
#@test norm(xx - [inert2rot(xx_inert[i],tt[i],sys) for i = eachindex(tt)]) < 1e-12


# seczoom
xlims, ylims, zlims = seczoom(sys)
@test  diff(xlims) ≈ diff(ylims) ≈ diff(zlims)