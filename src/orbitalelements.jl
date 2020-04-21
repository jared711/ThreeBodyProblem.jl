using LinearAlgebra

"""
    AP2a(rA, rP)

Compute the semimajor axis given apoapsis and periapsis distances rA and rP.

# Example
```jldoctest
julia> AP2a(1e6,1e5)
550000.0
```
"""
function AP2a(rA, rP)
    a = (rA + rP)/2
end

"""
    findrP(rv::Array,p::Array)
"""
function findrP(rv::Array,p::Array=[0,0,0])
    if size(rv[1]) == (6,) || size(rv[1]) == (3,)
        r = [rv[i][1:3] - p for i = 1:length(rv)]
        rP = minimum(norm,r)
    else
        rP = norm(rv[1:3] - p)
    end
end
