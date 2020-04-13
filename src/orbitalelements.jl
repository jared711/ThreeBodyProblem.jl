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
