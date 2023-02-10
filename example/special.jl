using LinearAlgebra
using OrdinaryDiffEq

"""
    invariant_manifolds(sys::System, rv0, T; tf=1., nPts=20, α=1e-6, reltol=1e-12, integrator=TsitPap8())

Compute invariant manifolds of a periodic orbit in the CR3BP.
"""
function invariant_manifolds(sys::System, rv₀, T; tf=1., nPts=20, α=1e-6, reltol=1e-12, integrator=TsitPap8())

    Φ₀ = I(6)
    w₀ = vcat(rv₀,reshape(Φ₀,36,1))
    tspan = (0.,T)
    prob = ODEProblem(CR3BPstm!,w₀,tspan,sys)
    nomTraj = solve(prob, integrator, reltol=reltol)
    Φₜ = Matrix(reshape(nomTraj.u[end][7:42],6,6))
    rvₜ = nomTraj.u[end][1:6] # last time step
    Λ,V = eigen(Φₜ,sortby=isreal) # Λ is vector of eigenvalues, V is matrix of eigenvectors
    Yw = real(V[:,findall(isreal, Λ)]) # Eigenvectors corresponding to real eigenvalues
    Λ = Λ[findall(isreal, Λ)] # Purely real eigenvalues (have 0.0 imaginary component)
    Yws = Yw[:,findmin(real(Λ))[2]] # Eigenvector associated with stable eigenvalue λ < 1
    Ywu = Yw[:,findmax(real(Λ))[2]] # Eigenvector associated with unstable eigenvalue λ > 1

    tspan_forward = (0., tf)
    tspan_backward = (0., -tf)

    # Wsp = []
    # Wsn = []
    # Wup = []
    # Wun = []
    # for i = 1:length(nomTraj)
    #   Φᵢ = Matrix(reshape(nomTraj.u[i][7:42],6,6)) # do I need Matrix()?
    #   rv0sp = nomTraj[i][1:6] + α*Φᵢ*Yws/norm(Φᵢ*Yws);
    #   rv0sn = nomTraj[i][1:6] - α*Φᵢ*Yws/norm(Φᵢ*Yws);
    #   rv0up = nomTraj[i][1:6] + α*Φᵢ*Ywu/norm(Φᵢ*Ywu);
    #   rv0un = nomTraj[i][1:6] - α*Φᵢ*Ywu/norm(Φᵢ*Ywu);
    #
    #    prob_sp = ODEProblem(CR3BPdynamics!,rv0sp, tspan_backward, sys) # stable positive
    #    push!(Wsp, solve(prob_sp, integrator, reltol=reltol))
    #    prob_sn = ODEProblem(CR3BPdynamics!,rv0sn, tspan_backward, sys) # stable negative
    #    push!(Wsn, solve(prob_sn, integrator, reltol=reltol))
    #    prob_up = ODEProblem(CR3BPdynamics!,rv0up, tspan_forward, sys) # unstable positive
    #    push!(Wup, solve(prob_up, integrator, reltol=reltol))
    #    prob_un = ODEProblem(CR3BPdynamics!,rv0un, tspan_forward, sys) # unstable negative
    #    push!(Wun, solve(prob_un, integrator, reltol=reltol))
    # end

    Wsp = []
    Wsn = []
    Wup = []
    Wun = []
    for t = LinRange(0,T,nPts)
      Φᵢ = Matrix(reshape(nomTraj(t)[7:42],6,6)) # do I need Matrix()?
      rv0sp = nomTraj(t)[1:6] + α*Φᵢ*Yws/norm(Φᵢ*Yws);
      rv0sn = nomTraj(t)[1:6] - α*Φᵢ*Yws/norm(Φᵢ*Yws);
      rv0up = nomTraj(t)[1:6] + α*Φᵢ*Ywu/norm(Φᵢ*Ywu);
      rv0un = nomTraj(t)[1:6] - α*Φᵢ*Ywu/norm(Φᵢ*Ywu);

      prob_sp = ODEProblem(CR3BPdynamics!,rv0sp, tspan_backward, sys) # stable positive
      push!(Wsp, solve(prob_sp, integrator, reltol=reltol))
      prob_sn = ODEProblem(CR3BPdynamics!,rv0sn, tspan_backward, sys) # stable negative
      push!(Wsn, solve(prob_sn, integrator, reltol=reltol))
      prob_up = ODEProblem(CR3BPdynamics!,rv0up, tspan_forward, sys) # unstable positive
      push!(Wup, solve(prob_up, integrator, reltol=reltol))
      prob_un = ODEProblem(CR3BPdynamics!,rv0un, tspan_forward, sys) # unstable negative
      push!(Wun, solve(prob_un, integrator, reltol=reltol))
   end

    return Wsp, Wsn, Wup, Wun, Φₜ, nomTraj
end

""" condition(u, t, integrator) = u[2]
 affect!(integrator) = terminate!(integrator)
 cb = ContinousCallback(condition, affect!)
 sol = solve(prob, Tsit5(), calback=cb)

    rich3()
"""

"""
   differential_corrector(sys::System, rv₀; myconst=3, iter=100, plot=false, t0=0., tf=1., dir=1, tol=1e-12)

Given a sufficiently close guess at a periodic trajectory rv₀, returns the corrected initial condition and period
"""
function differential_corrector(sys::System, rv₀; myconst=3, iter=100, plot=false, t0=0., tf=1., dir=1, tol=1e-12)
   Φ₀ = I(6)

   tspan = (t0, tf)

   # event function
   condition(u, t, integrator) = u[2]
   affect!(integrator) = terminate!(integrator)
   cb = OrdinaryDiffEq.ContinuousCallback(condition, affect!)

   for k = 1:iter
      w₀ = vcat(rv₀, reshape(Φ₀,36,1))

      prob = ODEProblem(CR3BPstm!,w₀,tspan,sys)
      sol = solve(prob, TsitPap8(), reltol=tol, callback=cb)

      w = sol[end]
      rv = w[1:6]
      Φ = reshape(w[7:42],6,6)
      global T = 2*sol.t[end]

      ∂ẋ = 0 - rv[4]
      ∂ż = 0 - rv[6]

      if norm([∂ẋ, ∂ż]) < tol
         break
      end
      if k == iter
         @warn("Differential corrector did not converge")
      end

      ẋ,ẏ,ż,ẍ,ÿ,z̈ = CR3BPdynamics(rv, sys, 0)

      if myconst == 1
         ∂Φ = [Φ[4,3] Φ[4,5];
               Φ[6,3] Φ[6,5]]
         dyad = [ẍ;z̈]*[Φ[2,3] Φ[2,5]]

         ∂z, ∂ẏ = (∂Φ - dyad/ẏ)\[∂ẋ; ∂ż]
         rv₀[3] += ∂z
         rv₀[5] += ∂ẏ
      elseif myconst == 3
         ∂Φ = [Φ[4,1] Φ[4,5];
               Φ[6,1] Φ[6,5]]
         dyad = [ẍ;z̈]*[Φ[2,1] Φ[2,5]]

         ∂x, ∂ẏ = (∂Φ - dyad/ẏ)\[∂ẋ; ∂ż]
         rv₀[1] += ∂x
         rv₀[5] += ∂ẏ
      elseif myconst == 5
         ∂Φ = [Φ[4,1] Φ[4,3];
               Φ[6,1] Φ[6,3]]
         dyad = [ẍ;z̈]*[Φ[2,1] Φ[2,5]]

         ∂x, ∂z = (∂Φ - dyad/ẏ)\[∂ẋ; ∂ż]
         rv₀[1] += ∂x
         rv₀[3] += ∂ẏ
      else
         error("myconst should be 1, 3, or 5")
      end

   end

   return rv₀, T
end

""" condition(u, t, integrator) = u[2]
 affect!(integrator) = terminate!(integrator)
 cb = ContinousCallback(condition, affect!)
 sol = solve(prob, Tsit5(), calback=cb)
"""

"""
   fourier_operator(N)

Returns the fourier operator D, the vector of indices k, and the vector of angles θ
"""
function fourier_operator(N)
   if iseven(N);  error("N must be odd"); end

   θ = 2π*(0:N-1)/N # Angles for the invariant circle
   k = Int(-(N-1)/2):Int((N-1)/2) # Vector of integers from -(N-1)/2 to (N-1)/2. N is odd, so (N-1)/2 is an integer.
   D = 1/N*exp.(-im*k*θ') # D is a constant matrix, as k and θ won't change
   return D, k, θ
end

"""
   rotation_operator(ρ, N)

Returns the rotation operator R(ρ) and its derivative ∂R/∂ρ(ρ) as well as Q(ρ) and ∂Q/∂ρ(ρ)
"""
function rotation_operator(ρ, N)
   D, k = fourier_operator(N) # D is the fft matrix and k is a vector of indices

   Q(ρ) = Diagonal(exp.(-im*k*ρ)) # Q is the rotation operator in the fourier domain. It's a diagonal matrix made up of exponential terms
   # since multiplying by e^(ikρ) rotates a point by kρ radians, we use e^(-ikρ) to rotate backwards by kρ radians.
   R(ρ) = real(inv(D)*Q(ρ)*D) # R is the rotation operator in the real domain. We use the similarity transform to convert it to the real domain.
   # We need to use real() to make sure each component is real. Multiplying UT by R(ρ) rotates the invariant circle by ρ radians
   ∂Q_∂ρ(ρ) = Diagonal(-im*k.*exp.(-im*k*ρ))
   ∂R_∂ρ(ρ) = real(inv(D)*∂Q_∂ρ(ρ)*D)

   return R, Q, ∂Q_∂ρ, ∂R_∂ρ, D
end


""" strob_map(rv₀, u₀, sys::System)
   
Given a base point rv₀ on a periodic orbit and and invariant circle u₀, returns the stroboscopic map uTR and the Jacobi constant
"""
function strob_map(rv₀, u, T, ρ, sys::System)
   N = length(u) # N is the number of points used in the invariant circle
   C = sum([computeC(rv₀+u[i],sys) for i in 1:N])/N # Compute Jacobi Constant of the invariant circle

   # First we integrate the orbit
   Φ₀ = I(6) # Initialization of the STM, Φ₀ = I
   w₀ = [rv₀; reshape(Φ₀,36,1)] # Reshape the matrix into a vector and append it to the state vector
   tspan = (0.,T) # integrate from 0 to T₀
   prob_halo = ODEProblem(CR3BPstm!,w₀,tspan,sys) # CR3BPstm! is our in-place dynamics function for state and STM
   function prob_func(prob, i, repeat) # function that defines the initial condition for each trajectory
      remake(prob, u0=[rv₀+u[i]; reshape(Φ₀,36,1)]) # perturb rv₀ by the ith point of the invariant circle and use that as the initial condition
   end # NOTE: Do not get confused! ODEProblems have an field called "u0" (e.g. prob_halo.u0) and the solved problems have a field called "u" (e.g. halo.u). Don't confuse these with the u₀ and u variables that we are defining here.
   prob_qpo = EnsembleProblem(prob_halo, prob_func=prob_func) # ODE problem with an ensemble of trajectories
   qpo = solve(prob_qpo, TsitPap8(), trajectories=N, abstol=1e-12, reltol=1e-12) # solve the problem

   uT = [qpo[i].u[end][1:6]-rv₀ for i in 1:N] # Invariant circle after integrating (make sure to subtract the base point rv₀)
   UT = reduce(vcat,uT') # UT is an Nx6 matrix form of uT

   R,_,_,_ = rotation_operator(ρ, N) # make sure to use the comma so R doesn't become an array of functions

   UTR = R(ρ)*UT # Rotate the integrated invariant circle back by ρ radians
   uTR = [UTR[i,:] for i in 1:N] # Convert back to a vector of vectors

   return uTR, C, qpo, UT
end

"""
   phase_constraints!(J, dγ, u, rv₀, u₀, T₀, ρ₀, N, sys)

Appends the phase constraints for the invariant circle
"""
function add_phase_constraints(J, dγ, u, rv₀, u₀, T₀, ρ₀, N, sys)
   n = length(rv₀)
   # Phase constraints
   D, k, θ = fourier_operator(N) # Compute the fast fourier transform matrix
   C0_tilde = D*reduce(vcat,u₀') # convert to an Nx6 matrix
   ∂u₀_∂θ₁ = im.*exp.(im.*θ*k')*Diagonal(k)*C0_tilde # Derivative of the initial invariant circle with respect to θ₁
   if maximum(imag(∂u₀_∂θ₁)) > eps() # If the derivative is complex, then we have a problem
      error("The derivative of the initial invariant circle with respect to θ₁ is complex")
   end
   ∂u₀_∂θ₁ = real(∂u₀_∂θ₁) # Enforce the derivative is real
   
   U̇₀ = zeros(N,n)
   for i = 1:N
      U̇₀[i,:] = CR3BPdynamics(rv₀ + u₀[i],sys,0)
   end
   ∂u₀_∂θ₀ = T₀/2π .* (U̇₀ - ρ₀/T₀.*∂u₀_∂θ₁) # Derivative of the initial invariant circle with respect to θ₀
   
   ∂u₀_∂θ₀ = reshape(∂u₀_∂θ₀',1,n*N) # Convert to a row vector
   ∂u₀_∂θ₁ = reshape(∂u₀_∂θ₁',1,n*N) # Convert to a row vector

   J = [J;        # Add the phase constraints to the Jacobian
   ∂u₀_∂θ₀ 0 0;
   ∂u₀_∂θ₁ 0 0]

   append!(dγ,∂u₀_∂θ₀*reduce(vcat,u)) # Add the first phase constraint to the error vector
   append!(dγ,∂u₀_∂θ₁*reduce(vcat,u)) # Add the second phase constraint to the error vector

   return J, dγ
end

"""
   add_state_constraints!(J, dγ, u; constraints=[])

add the state constraints to the Jacobian and error vector
"""
function add_state_constraints(J, dγ, u; constraints=[])
   N = length(u)
   n = length(u[1])
   for j in constraints
      row = zeros(1,n*N+2) # +2 for the time period and the phase constraint
      row[j] = 1 # Put a 1 in the row of zeros
      J = [J;row] # Add the row to the Jacobian
      append!(dγ,0) # Add a zero to the error vector
   end
   return J, dγ
end


"""
   differential_corrector_QPO(sys::System, rv₀, u₀, T₀, ρ₀; max_iter=10, plot_on=false, ϵ=1e-6, constraint::Symbol=:C)

Differential corrector for quasi-periodic orbits QPO problem. 
Takes in a System sys, periodic orbit state rv₀, invariant circle.
Returns the corrected state vector and the time period T.
"""
function differential_corrector_QPO(sys::System, rv₀, u₀, T₀, ρ₀; max_iter=10, plot_on=false, ϵ=1e-6, constraint::Symbol=:C)
   u = u₀ # u₀ is the invariant circle for which the phase constraint will be enforced
   T = T₀ # Initial guess at the period of the QPO
   ρ = ρ₀ # Initial guess at the rotation number of the QPO
   
   N = length(u) # Number of points along invariant circle 
   if iseven(N);  @error "N must be an odd number";   end # Should be an odd number

   n = length(rv₀) # Dimension of state vector (normally n = 6)

   C₀ = computeC(rv₀,sys) # Jacobi constant of central orbit
   C = sum([computeC(rv₀+u[i],sys) for i in 1:N])/N # Compute the average Jacobi constant across the invariant circle

   uTRs = []
   us = [u]
   Cs = [C]
   # while err > ϵ
   R, _, _, ∂R_∂ρ = rotation_operator(ρ,N) # ∂R_∂ρ is the fourth output from the rotation_operator function
   
   for _ in 1:max_iter
      # First we perform the stroboscopic mapping
      uTR, C, qpo, UT = strob_map(rv₀, u, T, ρ, sys) # uTR is the invariant circle after the stroboscopic map
      push!(uTRs,uTR) # Save the invariant circle after the stroboscopic map

      u_err = uTR - u # Compute the error between the initial and integrated/rotated invariant circles
      err = norm(u_err)
      if err < ϵ; # If the error is small enough, we're done
         @info "CONVERGED Differential corrector error = $err"
         break;
      else
         @info "Differential corrector error = $err"
      end

      # The error vector serves as our constraint. We want to drive dγ to zero.
      dγ = [reduce(vcat,u_err); # reduce(vcat,u_err) turns u_err into one big long vector instead of a vector of vectors 
                        C - C₀]
                        
      # The next step is to compute the Jacobian of the stroboscopic map with respect to the initial invariant circle u, rotation number ρ, and period T
      # Let's start with J₁ = ∂(uTR-u₀)/∂u
      Φ_tilde = zeros(n*N,n*N)
      for i = 1:N
         idx = (i-1)*n + 1:i*n
         Φ_tilde[idx,idx] = reshape(qpo[i].u[end][7:end],6,6) # Φ_tilde is made up of the state transition matrices of each point of the invariant circle (unrotated points, as the rotation operator doesn't depend on u and gets multiplied later)
      end
      ∂uTR_∂u = kron(R(ρ),I(6))*Φ_tilde # Compute the derivative of the integrated/rotated invariant circle with respect to the initial invariant circle (u₀) 
      J₁ = ∂uTR_∂u - I(n*N) # We subtract the identity because we really want ∂(uTR-u)/∂u

      # Next  J₂ = ∂uTR/∂T
      J₂ = zeros(n*N,1) # column vector of size n*N
      for i = 1:N
         idx = (i-1)*n + 1:i*n
         ẋ, ẏ, ż, ẍ, ÿ, z̈ = CR3BPdynamics(rv₀ + uTR[i],sys,0) # don't forget to add rv₀ to uTR[i]
         J₂[idx] = [ẋ, ẏ, ż, ẍ, ÿ, z̈] # The derivative with respect to time comes right from the equations of motion
      end
      
      # Next  J₃ = ∂uTR/∂ρ
      J₃ = ∂R_∂ρ(ρ)*UT # Compute the derivative of the rotation operator in the real domain
      J₃ = reshape(J₃',n*N,1) # Convert to a column vector

      # Finally J₄ = ∂C/∂u
      J₄ = zeros(1,n*N) # row vector of size n*N
      for i = 1:N
         idx = (i-1)*n + 1:i*n
         ẋ, ẏ, ż, ẍ, ÿ, z̈ = CR3BPdynamics(rv₀ + u[i],sys,0) # use u[i] here instead of uTR[i]
         Ωx = ẍ - 2ẏ
         Ωy = ÿ + 2ẋ
         Ωz = z̈
         J₄[idx] = [2Ωx, 2Ωy, 2Ωz, -2ẋ, -2ẏ, -2ż]/N # We divide by N because we want the derivative of C_avg with respect to u
      end

      J = [J₁  J₂  J₃;
           J₄   0   0];

      # J, dγ = add_phase_constraints(J, dγ, u, rv₀, u₀, T₀, ρ₀, N, sys)
      J, dγ = add_state_constraints(J, dγ, u, constraints = [1])

      dξ = -J\dγ
      du = [dξ[(i-1)*n + 1:i*n] for i = 1:N]
      dT = dξ[n*N+1]
      dρ = dξ[n*N+2]

      u += du
      T += dT
      ρ += dρ

      C = sum([computeC(rv₀+u[i],sys) for i in 1:N])/N # Compute the Jacobi constant for each state along the invariant circle

      push!(us,u)
      push!(Cs,C)

   end

   return u, T, ρ, C, us, uTRs, Cs
end


"""
   invariant_circle(rv, T, N, sys::System; α=1e-5)

Computes an approximation of the invariant circle around a given initial state `rv` 
for a periodic orbit with given period `T` using `N` points. The algorithm uses 
α as the step size. The system is given by `sys`.
"""
function invariant_circle(rv, T, N, sys::System; α=1e-5)
   Φₜ = monodromy(rv, T, sys) # Compute the monodromy matrix
   λ, V = eigen(Φₜ) # λ is a vector of eigenvalues and V is a matrix of eigenvectors

   sort_idx = sortperm(λ,by=x->(abs(imag(x)),abs(real(x)))) # Sort the eigenvalues by their imaginary magnitude, then by their real magnitude
   λ = λ[sort_idx]
   V = V[:,sort_idx]

   eig_idx = 6 # We want the 6th eigenvalue, which is the one with the largest imaginary part, corresponding to periodic motion
   ρ = real(-im*log(λ[eig_idx])) # Initial guess for the rotation number of the invariant circle

   θ = 2π*(0:N-1)/N # Angles for the invariant circle
   u = [α*(cos(θ[i])*real(V[:,eig_idx]) - sin(θ[i])*imag(V[:,eig_idx])) for i in 1:N] # Initial guess for the invariant circle
   
   return u, ρ
end

"""
   monodromy(rv₀, T, sys::System)

computes the monodromy matrix for a given periodic orbit with initial conditions
rv₀ and period T. Returns Φ(T,0).
"""
function monodromy(rv₀, T, sys::System)
   Φ₀ = I(6) # Initialization of the STM, Φ₀ = I
   w₀ = [rv₀; reshape(Φ₀,36,1)] # Reshape the matrix into a vector and append it to the state vector
   tspan = (0.,T) # integrate from 0 to T
   prob = ODEProblem(CR3BPstm!,w₀,tspan,sys) # CR3BPstm! is our in-place dynamics function for state and STM
   sol = solve(prob,TsitPap8(),abstol=1e-12,reltol=1e-12) # solve the problem
   Φₜ = reshape(sol[end][7:end],6,6) # The final STM or monodromy matrix M = Φ(T,0)
   return Φₜ
end

"""
    rich3()
"""
function rich3(sys::System, Az, Lpt, NS, npts=10)
    μ = sys.μ
    γ = gammaL(sys, Lpt)

    if Lpt == 1;        won =  1;   primary = 1-μ;
    elseif Lpt == 2;    won = -1;   primary = 1-μ;
    elseif Lpt == 3;    won =  1;   primary = -μ;
    end

    c = zeros(4)

    if Lpt == 3
      for N = 2:4
        c[N]= (1/γ^3)*( 1-μ + (-primary*γ^(N+1))/((1+γ)^(N+1)) );
      end
    else
      for N = 2:4
        c[N]= (1/γ^3)*( (won^N)*μ + ((-1)^N)*((primary)*γ^(N+1))/((1+(-won)*γ)^(N+1)) );
      end
    end

    # polylambda = [1, 0, (c[2]-2), 0, -(c[2]-1)*(1+2*c[2])];
    polylambda = [-(c[2]-1)*(1+2*c[2]), 0, c[2]-2, 0, 1];
    lambdaroots = roots(polylambda); # lambda = frequency of orbit
    ## Stuck here 6/22/21 (need to select the right root)
    λ = real(sort(lambdaroots, by = x -> abs(imag(x))))[1]

    # if Lpt==3
    #    λ = lambdaroots[1];
    # else
    #    λ = lambdaroots[1];
    # end

    k   = 2*λ/(λ^2 + 1 - c[2]);


    del = λ^2 - c[2];

    d1  = ((3*λ^2)/k)*(k*( 6*λ^2 -1) - 2*λ);
    d2  = ((8*λ^2)/k)*(k*(11*λ^2 -1) - 2*λ);

    a21 = (3*c[3]*(k^2 - 2))/(4*(1 + 2*c[2]));
    a22 = 3*c[3]/(4*(1 + 2*c[2]));
    a23 = -(3*c[3]*λ/(4*k*d1))*( 3*(k^3)*λ - 6*k*(k-λ) + 4);
    a24 = -(3*c[3]*λ/(4*k*d1))*( 2 + 3*k*λ );

    b21 = -(3*c[3]*λ/(2*d1))*(3*k*λ - 4);
    b22 = 3*c[3]*λ/d1;
    d21 = -c[3]/(2*λ^2);

    a31 = -(9*λ/(4*d2))*(4*c[3]*(k*a23 - b21) + k*c[4]*(4 + k^2)) + ((9*λ^2 + 1 -c[2])/(2*d2))*(3*c[3]*(2*a23 - k*b21) + c[4]*(2 + 3*k^2));
    a32 = -(1/d2)*( (9*λ/4)*(4*c[3]*(k*a24 - b22) + k*c[4]) + 1.5*(9*λ^2 + 1 - c[2])*( c[3]*(k*b22 + d21 - 2*a24) - c[4]) );

    b31 = (.375/d2)*( 8*λ*(3*c[3]*(k*b21 - 2*a23) - c[4]*(2 + 3*k^2)) + (9*λ^2 + 1 + 2*c[2])*(4*c[3]*(k*a23 - b21) + k*c[4]*(4 + k^2)) );
    b32 = (1/d2)*( 9*λ*(c[3]*(k*b22 + d21 - 2*a24) - c[4]) + .375*(9*λ^2 + 1 + 2*c[2])*(4*c[3]*(k*a24 - b22) + k*c[4]) );

    d31 = (3/(64*λ^2))*(4*c[3]*a24 + c[4]);
    d32 = (3/(64*λ^2))*(4*c[3]*(a23- d21) + c[4]*(4 + k^2));

    s1  = (1/(2*λ*(λ*(1+k^2) - 2*k)))*( 1.5*c[3]*(2*a21*(k^2 - 2)-a23*(k^2 + 2) - 2*k*b21) - .375*c[4]*(3*k^4 - 8*k^2 + 8) );
    s2  = (1/(2*λ*(λ*(1+k^2) - 2*k)))*( 1.5*c[3]*(2*a22*(k^2 - 2)+a24*(k^2 + 2) + 2*k*b22 + 5*d21) + .375*c[4]*(12 - k^2) );

    a1  = -1.5*c[3]*(2*a21+ a23 + 5*d21) - .375*c[4]*(12-k^2);
    a2  =  1.5*c[3]*(a24-2*a22) + 1.125*c[4];

    l1 = a1 + 2*(λ^2)*s1;
    l2 = a2 + 2*(λ^2)*s2;

    # ADDITIONAL TERMS FROM GEOMETRY CENTER PAPER
    b33 = -k/(16*λ)*(12*c[3]*(b21-2*k*a21+k*a23)+3*c[4]*k*(3*k^2-4)+16*s1*λ*(λ*k-1));
    b34 = -k/(8*λ)*(-12*c[3]*k*a22+3*c[4]*k+8*s2*λ*(λ*k-1));
    b35 = -k/(16*λ)*(12*c[3]*(b22+k*a24)+3*c[4]*k);

    deltan = 2 - NS;

    Ax  = sqrt( (-del - l2*Az^2)/l1 );
    omg = 1+s1*Ax^2+s2*Az^2;
    freq=λ*omg;
    period=2*pi/freq;

    rvi  = zeros(npts,6);
    ss   = zeros(npts,1);
    if npts > 1
       dtau1= 2*pi/(npts-1);
    else
       dtau1= 2*pi;
    end
    tau1 = 0;
    for i=1:npts
       x = a21*Ax^2 + a22*Az^2 - Ax*cos(tau1) + (a23*Ax^2 - a24*Az^2)*cos(2*tau1) + (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau1);
       y = k*Ax*sin(tau1) + (b21*Ax^2 - b22*Az^2)*sin(2*tau1) + (b31*Ax^3 - b32*Ax*Az^2)*sin(3*tau1);
       z = deltan*Az*cos(tau1) + deltan*d21*Ax*Az*(cos(2*tau1) - 3) + deltan*(d32*Az*Ax^2 - d31*Az^3)*cos(3*tau1);
       y_plus = (b33*Ax^3 + b34*Ax*Az^2 - b35*Ax*Az^2)*sin(tau1);
       y = y + y_plus;     # ADD EXTRA TERMS FROM G.C. PAPER

       xdot = freq*Ax*sin(tau1) - 2*freq*(a23*Ax^2-a24*Az^2)*sin(2*tau1) - 3*freq*(a31*Ax^3 - a32*Ax*Az^2)*sin(3*tau1);
       ydot = freq*(k*Ax*cos(tau1) + 2*(b21*Ax^2 - b22*Az^2)*cos(2*tau1) + 3*(b31*Ax^3 - b32*Ax*Az^2)*cos(3*tau1));
       zdot = - freq*deltan*Az*sin(tau1) - 2*freq*deltan*d21*Ax*Az*sin(2*tau1) - 3*freq*deltan*(d32*Az*Ax^2 - d31*Az^3)*sin(3*tau1);
       ydot_plus = freq*(b33*Ax^3 + b34*Ax*Az^2 - b35*Ax*Az^2)*cos(tau1);
       ydot = ydot_plus + ydot; # ADD EXTRA TERMS FROM G.C. PAPER

       rvi[i,:]= γ*[ (primary+γ*(-won+x))/γ, y, z, xdot, ydot, zdot];
       ss[i]   = tau1/freq;
       tau1=tau1+dtau1;
    end

    return ss, rvi, period, Ax
end
