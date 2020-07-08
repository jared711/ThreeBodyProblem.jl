include("util.jl")
"""
    S2I!(rv,θ)
Synodic (rotating) frame to inertial frame
"""
# function rot2inert!(rv, θ, μ; origin=0)
function S2I!(rv,θ)
    cθ, sθ = cos(θ), sin(θ)
    A = [cθ -sθ  0;
         sθ  cθ  0;
          0   0  1]
    B = [-sθ -cθ  0;
          cθ -sθ  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]
    rv[1:6] = C*rv
    return nothing
end

# function rot2inert!(rv, θ, p::Array; origin=0)
function S2I!(rv,t,p::Array)
    μ₁,μ₂,d = p # parameters
    ωₛ = sqrt((μ₁ + μ₂)/d^3);
    cθ, sθ = cos(ωₛ*t), sin(ωₛ*t)
    A = [cθ -sθ  0;
         sθ  cθ  0;
          0   0  1]
    B = [-sθ -cθ  0;
          cθ -sθ  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]
    rv[1:6] = C*rv
    return nothing
end
