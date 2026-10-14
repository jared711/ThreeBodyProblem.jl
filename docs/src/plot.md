# [Plotting](@id plot)
```@meta
CurrentModule = ThreeBodyProblem
```

## Overview
ThreeBodyProblem.jl defines plot recipes, so once Plots.jl is loaded a `System`, a `Body` or a trajectory
(a `Vector{Vector{<:Real}}` of states) can be passed straight to `plot`:

```julia
using ThreeBodyProblem, Plots
sys = earth_moon()
plot(sys)               # primaries and Lagrange points in the rotating frame
plot(sys, planar=false) # 3D
```

## Recipes
```@autodocs
Modules = [ThreeBodyProblem]
Pages = ["plot.jl"]
Filter = t -> t isa Function && nameof(t) === :apply_recipe
```

## Functions
```@docs
circle
sphere
torus
seczoom
```

## Internals
```@docs
trajectory_2D
trajectory_3D
velocity_2D
```
