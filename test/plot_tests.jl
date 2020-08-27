using ThreeBodyProblem
using RecipesBase
using Test

r = rand(Float64)
x,y = circle(r)
@test x[1]^2 + y[1]^2 == r^2

x,y,z = sphere(r)
@test x[1]^2 + y[1]^2 + z[1]^2 == r^2

PRIM, SEC, SYS = earth_moon()
rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), SYS)