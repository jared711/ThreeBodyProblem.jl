include("utils.jl")

function rot2inert!(rv, θ, μ; origin=0)
    n = size(rv)[1]
    cθ, sθ = cos(θ), sin(θ)
    A = [cθ -sθ  0;
         sθ  cθ  0;
          0   0  1]
    B = [-sθ -cθ  0;
          cθ -sθ  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]
    rv = C*rv
    if origin == 0      #Barycentric Origin
    elseif origin == 1  #Primary Body Origin
        r_offset = -μ*[cθ;sθ;0]
        v_offset = -μ*[-sθ;cθ;0]
        rv_offset = [r_offset;v_offset]
        rv = rv - rv_offset
    elseif origin == 2 #Secondary Body Origin
        r_offset = (1-μ)*[cθ;sθ;0]
        v_offset = (1-μ)*[-sθ;cθ;0]
        rv_offset = [r_offset;v_offset]
        rv = rv - rv_offset
    else
        error("'origin' must be 0, 1, or 2")
    end
end

function rot2inert!(rv, θ, p::Array; origin=0)
    n = size(rv)[1]
    R₁,R₂ = findR1R2(p)
    cθ, sθ = cos(θ), sin(θ)
    A = [cθ -sθ  0;
         sθ  cθ  0;
          0   0  1]
    B = [-sθ -cθ  0;
          cθ -sθ  0;
           0   0  0]
    C = [A zeros(3,3);
         B A]
    rv = C*rv
    if origin == 0      #Barycentric Origin
    elseif origin == 1  #Primary Body Origin
        r_offset = -R₁*[cθ;sθ;0]
        v_offset = -R₁*[-sθ;cθ;0]
        rv_offset = [r_offset;v_offset]
        rv = rv - rv_offset
    elseif origin == 2 #Secondary Body Origin
        r_offset = R₂*[cθ;sθ;0]
        v_offset = R₂*[-sθ;cθ;0]
        rv_offset = [r_offset;v_offset]
        rv = rv - rv_offset
    else
        error("'origin' must be 0, 1, or 2")
    end
end
