# using DifferentialEquations
# using LinearAlgebra

"""
    invariantManifolds(SYS, rv0, T, tf, nPts)

Compute invariant manifolds of a periodic orbit in the CR3BP.
"""
function invariantManifolds(sys::System, rv₀, T, tf, nPts)
    Φ₀ = I(6)
    w₀ = [reshape(Φ₀,36,1);rv₀]
    tspan = (0,)
    prob = ODEProblem(CR3BPstm!,w₀,(),sys)
    sol = solve(prob,reltol=1e-6)
    return sol
end

function differentialCorrector()
end
