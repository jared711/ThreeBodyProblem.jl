using DifferentialEquations
using LinearAlgebra

"""
    invariantManifolds(SYS, rv0, T, tf, nPts)

Compute invariant manifolds of a periodic orbit in the CR3BP.
"""
function invariant_manifolds(sys::System, rv₀, T, tf, nPts)
    Φ₀ = I(6)
    w₀ = [reshape(Φ₀,36,1);rv₀]
    tspan = (0.,T)
    prob = ODEProblem(CR3BPstm!,w₀,tspan,sys)
    sol = solve(prob, reltol=1e-6)
    Φₜ = Matrix(reshape(sol.u[end][1:36],6,6))
    rvₜ = sol.u[end][37:42] # last time step
    Λ,V = eigen(Φₜ,sortby=isreal) # Λ is vector of eigenvalues, V is matrix of eigenvectors
    Yw = real(V[:,findall(isreal, Λ)]) # Eigenvectors corresponding to real eigenvalues
    Λ = Λ[findall(isreal, Λ)] # Purely real eigenvalues (have 0.0 imaginary component)
    Yws = Yw[:,findmin(real(Λ))[2]] # Eigenvector associated with stable eigenvalue λ < 1
    Ywu = Yw[:,findmax(real(Λ))[2]] # Eigenvector associated with unstable eigenvalue λ > 1

    α = 1e-6
    rv0sp = rvₜ + α*Φₜ*Yws/norm(Φₜ*Yws);
    rv0sn = rvₜ - α*Φₜ*Yws/norm(Φₜ*Yws);
    rv0up = rvₜ + α*Φₜ*Ywu/norm(Φₜ*Ywu);
    rv0un = rvₜ - α*Φₜ*Ywu/norm(Φₜ*Ywu);

    tspan_forward = (0., tf)
    tspan_backward = (0., -tf)

    myreltol = 1e-12
    prob_sp = ODEProblem(CR3BPdynamics!,rv0sp,tspan_backward,sys) # stable positive
    Wsp = solve(prob_sp, reltol=myreltol)
    prob_sn = ODEProblem(CR3BPdynamics!,rv0sn,tspan_backward,sys) # stable negative
    Wsn = solve(prob_sn, reltol=myreltol)
    prob_up = ODEProblem(CR3BPdynamics!,rv0up,tspan_backward,sys) # unstable positive
    Wup = solve(prob_up, reltol=myreltol)
    prob_un = ODEProblem(CR3BPdynamics!,rv0un,tspan_backward,sys) # unstable negative
    Wun = solve(prob_un, reltol=myreltol)

    return Wsp, Wsn, Wup, Wun, Φₜ
end

function differential_corrector()
end

"""
"""
function rich3()
end
