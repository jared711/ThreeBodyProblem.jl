# Learn how to compute Quasi-periodic orbits in the CR3BP

using ThreeBodyProblem
using OrdinaryDiffEq
using Plots
using LinearAlgebra

# We'll first define a halo orbit in the Saturn/Enceladus system

# Define the Saturn/Enceladus CR3BP system
sys = saturn_enceladus()

# Initial Conditions of the halo orbit (obtained from https://ssd.jpl.nasa.gov/tools/periodic_orbits.html)
rv₀ = [ 1.002850044069033
                        0
        0.004808592996246
       -0.000000000000001
       -0.005747591930694
       -0.000000000000006]
# rv₀ = [ 0.996791592030909;
#                         0;
#         0.004595185232449;
#                         0;
#         0.006025851721922;
#         0.000000000000001]


# Period of the halo orbit
T₀ = 2.445783783357601

# We'll use the differential corrector to refine the initial conditions (our μ may be slightly different)
rv₀, T₀ = differential_corrector(sys, rv₀, tf=T₀) # Differential corrector is a function of ThreeBodyProblem.jl

# Jacobi constant of the halo orbit
C₀ = computeC(rv₀,sys)

# Integrate the halo orbit and its state transition matrix (STM)
Φ₀ = I(6) # Initialization of the STM, Φ₀ = I
w₀ = [rv₀; reshape(Φ₀,36,1)] # Reshape the matrix into a vector and append it to the state vector
tspan = (0.,T₀) # integrate from 0 to T₀
prob_halo = ODEProblem(CR3BPstm!,w₀,tspan,sys) # CR3BPstm! is our in-place dynamics function for state and STM
halo = solve(prob_halo,abstol=1e-12,reltol=1e-12) # solve the problem

# Plot the halo orbit
pxy = plot(halo,vars=(1,2),legend=false,title="Halo Orbit",xaxis="x",yaxis="y"); # plot the halo orbit in the x-y plane
pyz = plot(halo,vars=(2,3),legend=false,title="Halo Orbit",xaxis="y",yaxis="z"); # plot the halo orbit in the y-z plane
pxz = plot(halo,vars=(1,3),legend=false,title="Halo Orbit",xaxis="x",yaxis="z"); # plot the halo orbit in the x-z plane
pall = plot(halo,vars=(1,2,3),legend=false,title="Halo Orbit"); # plot the halo orbit in 3D
plot_halo = plot(pxy,pyz,pxz,pall,layout=(2,2)); # plot all of the plots in a 2x2 grid
display(plot_halo)

# OK, now we have the halo orbit and its STM. We want to approximate a quasi-halo orbit around the halo orbit.
# A quasi-halo orbit lives on the surface of a torus, so it is a 2-dimensional object rather than a 1-dimensional curve like the halo orbit.
# It's tricky working with surfaces. For example, how do you define how closely two surfaces match each other?
# We'd much rather work with a curve, so let's take a slice of the torus at one moment in time, say t=0.
# This slice will be a 1D curve called the invariant circle.
# We can approximate it by using the monodromy matrix, which is the STM at the end of the integration Φ(T,0)

# Let's pull out the monodromy matrix from our halo object
wf = halo[end] # The final state vector appended with the final STM
M = reshape(wf[7:end],6,6) # The final STM or monodromy matrix M = Φ(T,0)

# We need to compute the eigenvalue decomposition of the monodromy matrix
λ, V = eigen(M) # λ is a vector of eigenvalues and V is a matrix of eigenvectors

sort_idx = sortperm(λ,by=x->(abs(imag(x)),abs(real(x)))) # Sort the eigenvalues by their imaginary magnitude, then by their real magnitude
V = V[:,sort_idx]
λ = λ[sort_idx]
# after sorting, the first two eigenvalues are on the real axis, the next two are nearly identity, and the last two are complex conjugate pairs on the unit circle



# Plot these eigenvalues on the complex plan
plot_eig = scatter(real(λ),imag(λ),legend=false,title="Eigenvalues of the Monodromy Matrix",xaxis="Real",yaxis="Imaginary",aspect_ratio=1,marker=:x);
plot!(plot_eig, circle(), seriestype = [:path,], lw = 0.5, linecolor = :black, label = "Unit Circle"); # circle() is a function in ThreeBodyProblem.jl that produces the coordinates of the unit circle
# You can see that there are two eigenvalues on the real axis. They are a reciprocal pair and correspond to stable and unstable motion.
annotate!(plot_eig, real(λ[1]), imag(λ[1]), text("Stable",:top, :right));
annotate!(plot_eig, real(λ[2]), imag(λ[2]), text("Unstable",:bottom, :right));
# There is another pair very near to 1. Every purely periodic orbit has two eigenvalues that are 1, which correspond to motion along the periodic orbit itself. Any deviation from 1 is caused by numerical imprecision
annotate!(plot_eig, real(λ[3]), imag(λ[3]), text("Identity",:top, :left));
# The remaining two eigenvalues are complex conjugate pairs on the unit circle. They correspond to periodic motion around the halo orbit. These are the ones we will use.
annotate!(plot_eig, real(λ[6]), imag(λ[6]), text("Periodic",:top, :left))
display(plot_eig)

# Because we sorted the eigenvalues by their imaginary magnitude, we should always be able to use the last one to find periodic motion
eig_idx = 6 # Index of the eigenvalue we want to use
# Note, sometimes, there are two pairs of eigenvalues on the unit circle. In that case, this method will use the one with the largest imaginary magnitude
# If there are two pairs of eigenvalues on the unit circle, then the orbit is stable and there exist two periodic manifolds, so we could compute a 3-torus. But let's not worry about that for now.

# Compute the approximate invariant circle
N = 19 # Number of points on the invariant circle (THIS SHOULD BE AN ODD NUMBER!!!)
# N is also the number of frequencies that we will break our u function into
n = 6 # Number of dimensions of the system

θ = 2π*(0:N-1)/N # Angles for the invariant circle
α = 1e-5 # parameter to control the size of the invariant circle
u₀ = [α*(cos(θ[i])*real(V[:,eig_idx]) - sin(θ[i])*imag(V[:,eig_idx])) for i in 1:N] # Initial guess for the invariant circle
plot_u = plot(u₀, legend=true,label="u₀",title="Approximate Invariant Circle",linecolor=:blue); # Plot the invariant circle
scatter!(plot_u, [u₀[1][1]],[u₀[1][2]],[u₀[1][3]],label="u₀[1]",shape=:o,markercolor=:blue) # Plot a marker on the first point of the invariant circle
display(plot_u)

# An invariant circle is defined as a set of points that are mapped to themselves under the dynamics
# So if this was the true invariant circle, then integrating this ring of points should return the exact same points
# As an added point of subtlety, the points will have rotated by some angle ρ during that integration
# ρ is referred to as the rotation number of the invariant circle

# Use the periodic eigenvalue to guess the rotation number of the invariant circle
ρ₀ = real(-im*log(λ[eig_idx])) # Initial guess for the rotation number of the invariant circle
# ρ₀ = atan(imag(λ[eig_idx]),real(λ[eig_idx])) # Another method for computing the initial guess for the rotation number of the invariant circle
# ρ₀ is simply the angle of the eigenvalue in the complex plane

# Now integrate the invariant circle and plot the result
# In Julia, we can use the EnsembleProblem to solve the same problem N times with different initial conditions
function prob_func(prob, i, repeat)
    remake(prob, u0=[rv₀+u₀[i]; reshape(Φ₀,36,1)]) # perturb rv₀ by the ith point of the invariant circle and use that as the initial condition
end # NOTE: Do not get confused! ODEProblems have an field called "u0" (e.g. prob_halo.u0) and the solved problems have a field called "u" (e.g. halo.u). Don't confuse these with the u₀ and u variables that we are defining here.
prob_qpo = EnsembleProblem(prob_halo, prob_func=prob_func) # ODE problem with an ensemble of trajectories
qpo = solve(prob_qpo, trajectories=N, abstol=1e-12, reltol=1e-12) # solve the problem

uT = [qpo[i].u[end][1:6]-rv₀ for i in 1:N] # Invariant circle after integrating (make sure to subtract the base point rv₀)
plot!(plot_u, uT,legend=true,label="uT",linecolor=:red); # Plot the invariant circle after integrating
scatter!(plot_u, [uT[1][1]],[uT[1][2]],[uT[1][3]],label="uT[1]",shape=:o,markercolor=:red); # Plot the first point of the integrated invariant circle
display(plot_u)
# You can see that the invariant circle doesn't quite match the integrated invariant circle
# You can also see the rotation number at work, as the first point of the integrated invariant circle is rotated by about ρ₀ radians

# This is called a stroboscopic map. It takes points x(t) and maps them to the points x(t+T)
# The stroboscopic map of a periodic orbit is a single point, because x(t) = x(t+nT)
# For a quasi-periodic orbit, the stroboscopic map is a circle. Therefore, our differential corrector will have to compare circles
# If the initial circle and the integrated circle are the same, then the circle is "invariant" and we have found the quasi-periodic orbit

# To truly compare the two circles, we need to rotate the integrated invariant circle back by ρ radians
# We'll create a rotation operator R to do this. 

# We can't just use an ordinary rotation matrix, as the invariant circle may be elliptical and lives in 6D phase space rather than 3D position space.
# Instead, we'll use the discrete Fourier transform the invariant circle in the Fourier domain, apply the rotation there, then transform it back to the real domain.

# D is the discrete Fourier transform matrix (make sure to use exp.() rather than exp(), 
# as we are applying the exponential operator to each component rather than performing the matrix exponential)
k = Int(-(N-1)/2):Int((N-1)/2) # Vector of integers from -(N-1)/2 to (N-1)/2. N is odd, so (N-1)/2 is an integer.
D = 1/N*exp.(-im*k*θ') # D is a constant matrix, as k and θ won't change
# The D matrix is size NxN, so we would like to choose as small of an N as possible

# Let's look at what the D matrix does to the invariant circle
UT = reduce(vcat,uT') # UT is an Nx6 matrix form of uT
pθ = plot(θ, UT, xticks = ([0:π/2:2*π;], ["0","\\pi/2","\\pi","3 \\pi/2","2\\pi"]),xlabel="θ [rad]",ylabel="",legend=true,title="Angular Domain",label=["x" "y" "z" "dx/dt" "dy/dt" "dz/dt"]); # Plot the result
uT_fourier = D*UT # uT_fourier is the invariant circle in the Fourier domain
pf = plot(k*ρ₀, abs.(uT_fourier),xlabel="frequency [rad/s]",ylabel="magnitude",legend=true,title="Fourier Domain",label=["x" "y" "z" "ẋ" "ẏ" "ż"]); # Plot the result
plot_fourier = plot(pθ,pf,layout=(1,2),size=(1000,400)); # Plot the two plots side by side
display(plot_fourier)


# Now that we have our invariant circle converted to the Fourier domain, let's rotate it backwards by ρ radians.
# Note, we need to rotate it backwards becasuse ρ is the amount it has rotated forward through integration. We need to reverse this rotation.
Q(ρ) = Diagonal(exp.(-im*k*ρ)) # Q is the rotation operator in the fourier domain. It's a diagonal matrix made up of exponential terms
# since multiplying by e^(ikρ) rotates a point by kρ radians, we use e^(-ikρ) to rotate backwards by kρ radians.
R(ρ) = real(inv(D)*Q(ρ)*D) # R is the rotation operator in the real domain. We use the similarity transform to convert it to the real domain.
# We need to use real() to make sure each component is real. Multiplying UT by R(ρ) rotates the invariant circle by ρ radians

UTR = R(ρ₀)*UT # Rotate the integrated invariant circle back by ρ₀ radians
uTR = [UTR[i,:] for i in 1:N] # Convert back to a vector of vectors
# uTR is the integrated invariant circle rotated back by ρ radians
# So we should be able to compare u₀ with uTR. If they are the same, then we have found the quasi-periodic orbit
plot!(plot_u, uTR,legend=true,label="uTR",linecolor=:magenta); # Plot the invariant circle after integrating and rotating
scatter!(plot_u, [uTR[1][1]],[uTR[1][2]],[uTR[1][3]],label="uTR[1]",shape=:o,markercolor=:magenta) # Plot the first point of the rotated, integrated invariant circle
display(plot_u)

# We can see that the rotation operator has done its job, as the first point of the integrated invariant circle is now rotated back to the same angle as the first point of the initial invariant circle
# We define the error term u_err as the difference between the initial and integrated/rotated invariant circles
u_err = uTR - u₀ # Compute the error between the initial and integrated/rotated invariant circles

# Let's see how big the error is. If it's above a certain threshold, we'll have to do some differential correction
err = norm(u_err)

# That error is small but not small enough. We can see that the invariant circles don't line up perfectly yet.
# We started with a small perturbation α. If we increase α, the error will increase.
ϵ = α/10 # Set the threshold for the error to be 1/10 of α

# Now we're going to do some differential correction to refine the invariant circle and reduce the error
# while err > ϵ
    
# First we construct the constraint vector
Cs = [computeC(rv₀+uTR[i],sys) for i in 1:N] # Compute the Jacobi constant for each state along the invariant circle
Cavg = sum(Cs)/N # Compute the average Jacobi constant across the invariant circle
dγ = [reduce(vcat,u_err); # reduce(vcat,u_err) turns u_err into one big long vector instead of a vector of vectors 
               Cavg - C₀]
            
# Let's construct the free variable vector ξ (That's greek letter "xi", type "\xi<tab>" in the repl to get it)
# these are the variables that we can change to reduce the error
ξ = [reduce(vcat,u₀); # reduce(vcat,u₀) turns u₀ into one big long vector instead of a vector of vectors 
                  ρ₀;
                  T₀]
# TODO we actually don't need to create ξ, we just need to find dξ, but it's useful to visualize what's in the free variable vector

# Now we construct the Jacobian matrix J = ∂γ/∂ξ
# γ(ξ + dξ) = γ(ξ) + J*dξ -> dγ = J*dξ -> dξ = J\dγ
# J tells us how a small change in the free variables will affect the constraint vector
# We can use this to compute the change in the free variables that will reduce the error the most
# J will be of the form
# J = [∂(uTR-u₀)/∂u  ∂(uTR-u₀)/∂ρ  ∂(uTR-u₀)/∂T;
#             ∂C/∂u         ∂C/∂ρ         ∂C/∂T]
# Luckily, several of the terms cancel out to zero
# ∂u₀/∂ρ = 0, ∂u₀/∂T = 0, ∂C_∂ρ = 0, and ∂C_∂T = 0
# Leaving us with
# J = [∂(uTR-u₀)/∂u  ∂uTR/∂ρ  ∂uTR/∂T;
#             ∂C/∂u        0        0]
# Let's call each of these blocks J₁ -  J₄
# J = [J₁  J₂  J₃;
#       J₄  0  0]

# Let's start with J₁ = ∂(uTR-u₀)/∂u
Φ_tilde = zeros(6N,6N)
for i = 1:N
    idx = (i-1)*n + 1:i*n
    Φ_tilde[idx,idx] = reshape(qpo[i].u[end][7:end],6,6) # Φ_tilde is made up of the state transition matrices of each point of the invariant circle (unrotated points, as the rotation operator doesn't depend on u and gets multiplied later)
end
∂uTR_∂u = kron(R(ρ₀),I(6))*Φ_tilde # Compute the derivative of the integrated/rotated invariant circle with respect to the initial invariant circle (u₀) 
J₁ = ∂uTR_∂u - I(6N) # We subtract the identity because we really want ∂(uTR-u)/∂u

# Next  J₂ = ∂uTR/∂ρ
# Let's note that uTR = R(ρ)*F(u₀). The only thing dependent on ρ is the rotation operator R(ρ)
# ∂uTR/∂ρ = ∂R/∂ρ*F(u₀)
# Note also that R(ρ) = inv(D)*Q(ρ)*D, and D is not dependent on ρ either
# ∂R/∂ρ = inv(D)*∂Q/∂ρ*D
# Recall that Q(ρ) = Diagonal(exp.(-im*k*ρ)). Luckily taking derivatives of exponents is very easy
# ∂Q/∂ρ = Diagonal(-im*k.*exp.(-im*k*ρ))
∂Q_∂ρ(ρ) = Diagonal(-im*k.*exp.(-im*k*ρ)) # Compute the derivative of the rotation operator in the Fourier domain
 J₂ = real(inv(D)*∂Q_∂ρ(ρ₀)*D)*UT # Compute the derivative of the rotation operator in the real domain
 J₂ = reshape( J₂',6N,1) # Convert to a column vector

# Next  J₃ = ∂uTR/∂T and  J₄ = ∂C/∂u
 J₃ = zeros(6N,1)
 J₄ = zeros(1,6N)
for i = 1:N
    idx = (i-1)*n + 1:i*n
    # TODO Do I compute the time derivative of uTR or uT?
    # Compute the time derivative of each state along the integrated/rotated invariant circle
    ẋ, ẏ, ż, ẍ, ÿ, z̈ = CR3BPdynamics(rv₀ + uTR[i],sys,0) # don't forget to add rv₀ to uTR[i]
    J₃[idx] = [ẋ, ẏ, ż, ẍ, ÿ, z̈] # The derivative with respect to time comes right from the equations of motion
    
    # Recall that C = 2Ω(r) - vᵀv
    # Ω is dependent only on position r, and vᵀv is of course only dependent on velocity v
    # The CR3BP equations of motion are as follows
    # ̈x = Ωx + 2̇y 
    # ̈y = Ωy - 2̇x 
    # ̈z = Ωz 
    # So we can compute Ωx, Ωy, and Ωz from our dot and ddot terms
    Ωx = ẍ - 2ẏ
    Ωy = ÿ + 2ẋ
    Ωz = z̈
    J₄[idx] = [2Ωx 2Ωy 2Ωz 2ẋ -2ẏ -2ż]
end
J₄ = J₄/N # We divide by N because we want the derivative of C_avg with respect to u

J = [J₁  J₂  J₃;
      J₄  0  0]

dξ = J\dγ
du = [dξ[(i-1)*n + 1:i*n] for i = 1:N]
dρ = dξ[6N+1]
dT = dξ[6N+2]

u = u₀ + du
ρ = ρ₀ + dρ
T = T₀ + dT

Cs = [computeC(rv₀+u[i],sys) for i in 1:N]
C = sum(Cs)/N

plot!(plot_u, u,legend=true,label="u1",linecolor=:black); # Plot the invariant circle after integrating
scatter!(plot_u, [u[1][1]],[u[1][2]],[u[1][3]],label="u1[1]",shape=:o,markercolor=:black); # Plot the first point of the integrated invariant circle
display(plot_u)

# Let's see if that new guess for the invariant circle is any better
plot_check = plot(u,legend=true,label="u1",linecolor=:blue); # Plot the original guess for the invariant circle
scatter!(plot_check, [u[1][1]],[u[1][2]],[u[1][3]],label="u[1]",shape=:o,markercolor=:blue); # Plot the first point of the rotated, integrated invariant circle

# Integrate the invariant circle
function prob_func(prob, i, repeat)
    remake(prob, u0=[rv₀+u[i]; reshape(Φ₀,36,1)]) # perturb rv₀ by the ith point of the invariant circle and use that as the initial condition
end # NOTE: Do not get confused! ODEProblems have an field called "u0" (e.g. prob_halo.u0) and the solved problems have a field called "u" (e.g. halo.u). Don't confuse these with the u₀ and u variables that we are defining here.
prob_qpo = EnsembleProblem(prob_halo, prob_func=prob_func) # ODE problem with an ensemble of trajectories
qpo = solve(prob_qpo, trajectories=N, abstol=1e-12, reltol=1e-12) # solve the problem
uT = [qpo[i].u[end][1:6]-rv₀ for i in 1:N] # Invariant circle after integrating
UT = reduce(vcat,uT') # convert to a Nx6 matrix
UTR = R(ρ)*UT # Rotate the integrated invariant circle back by ρ radians
uTR = [UTR[i,:] for i in 1:N] # Convert back to a vector of N vectors
plot!(plot_check, uTR,legend=true,label="uTR",linecolor=:magenta); # Plot the invariant circle after integrating and rotating
scatter!(plot_check, [uTR[1][1]],[uTR[1][2]],[uTR[1][3]],label="uTR[1]",shape=:o,markercolor=:magenta); # Plot the first point of the rotated, integrated invariant circle
display(plot_check)

# what's the error for this plot?
u_err = uTR - u₀ # Compute the error between the initial and integrated/rotated invariant circles
err = norm(u_err)

# Hmm, the error still isn't small enough. Let's try iterating

# Plot each iteration in a different color
plot_iter = plot(u₀,legend=true,label="u₀",linecolor=:blue, title="Iterations"); # Plot the original guess for the invariant circle
plot!(plot_iter, u, legend=true,label="u1",linecolor=:black); # Plot the updated invariant circle

# while err > ϵ
for iter = 2:3
    # plot_iter = plot(u,legend=true,label="u₀",linecolor=:blue) # Plot the invariant circle after integrating
    # scatter!(plot_iter, [u1[1][1]],[u1[1][2]],[u1[1][3]],label="uT[1]",shape=:o,markercolor=:red) # Plot the first point of the integrated invariant circle


    function prob_func(prob, i, repeat)
        remake(prob, u0=[rv₀+u[i]; reshape(Φ₀,36,1)]) # perturb rv₀ by the ith point of the invariant circle and use that as the initial condition
    end # NOTE: Do not get confused! ODEProblems have an field called "u0" (e.g. prob_halo.u0) and the solved problems have a field called "u" (e.g. halo.u). Don't confuse these with the u₀ and u variables that we are defining here.
    prob_qpo = EnsembleProblem(prob_halo, prob_func=prob_func) # ODE problem with an ensemble of trajectories
    qpo = solve(prob_qpo, trajectories=N, abstol=1e-12, reltol=1e-12) # solve the problem    
    uT = [qpo[i].u[end][1:6]-rv₀ for i in 1:N] # Invariant circle after integrating
    UT = reduce(vcat,uT') # convert to a Nx6 matrix
    UTR = R(ρ)*UT # Rotate the integrated invariant circle back by ρ radians
    uTR = [UTR[i,:] for i in 1:N] # Convert back to a vector of N vectors
    
    u_err = uTR - u # Compute the error between the initial and integrated/rotated invariant circles
    err = norm(u_err)

    Cs = [computeC(rv₀+uTR[i],sys) for i in 1:N] # Compute the Jacobi constant for each state along the invariant circle
    Cavg = sum(Cs)/N # Compute the average Jacobi constant across the invariant circle
    dγ = [reduce(vcat,u_err); # reduce(vcat,u_err) turns u_err into one big long vector instead of a vector of vectors 
                    Cavg - C]
        
    # Let's start with J₁ = ∂(uTR-u₀)/∂u
    Φ_tilde = zeros(6N,6N)
    for i = 1:N
        idx = (i-1)*n + 1:i*n
        Φ_tilde[idx,idx] = reshape(qpo[i].u[end][7:end],6,6) # Φ_tilde is made up of the state transition matrices of each point of the invariant circle (unrotated points, as the rotation operator doesn't depend on u and gets multiplied later)
    end
    ∂uTR_∂u = kron(R(ρ),I(6))*Φ_tilde # Compute the derivative of the integrated/rotated invariant circle with respect to the initial invariant circle (u₀) 
    J₁ = ∂uTR_∂u - I(6N) # We subtract the identity because we really want ∂(uTR-u)/∂u

    # Next  J₂ = ∂uTR/∂ρ
    ∂Q_∂ρ(ρ) = Diagonal(-im*k.*exp.(-im*k*ρ)) # Compute the derivative of the rotation operator in the Fourier domain
     J₂ = real(inv(D)*∂Q_∂ρ(ρ)*D)*UT # Compute the derivative of the rotation operator in the real domain
     J₂ = reshape( J₂',6N,1) # Convert to a column vector

    # Next  J₃ = ∂uTR/∂T and  J₄ = ∂C/∂u
     J₃ = zeros(6N,1)
     J₄ = zeros(1,6N)
    for i = 1:N
        idx = (i-1)*n + 1:i*n
        ẋ, ẏ, ż, ẍ, ÿ, z̈ = CR3BPdynamics(rv₀ + uTR[i],sys,0) # don't forget to add rv₀ to uTR[i]
        J₃[idx] = [ẋ, ẏ, ż, ẍ, ÿ, z̈] # The derivative with respect to time comes right from the equations of motion
        Ωx = ẍ - 2ẏ
        Ωy = ÿ + 2ẋ
        Ωz = z̈
        J₄[idx] = [2Ωx 2Ωy 2Ωz 2ẋ -2ẏ -2ż]
    end

    J = [J₁  J₂  J₃;
         J₄   0   0]

    dξ = J\dγ
    du = [dξ[(i-1)*n + 1:i*n] for i = 1:N]
    dρ = dξ[6N+1]
    dT = dξ[6N+2]

    u += du
    ρ += dρ
    T += dT

    Cs = [computeC(rv₀+u[i],sys) for i in 1:N]
    C = sum(Cs)/N

    plot!(plot_iter,u,legend=true,label=string("u",iter)) # Plot the invariant circle after integrating

end
display(plot_iter)

# If we plot the trajectories from the invariant circle, it will look just like the periodic orbit. That's because the invariant circle is very small.
plot(qpo,vars=(1,2,3),legend=false,title="Approximate QPO") # Plot the approximate quasi-periodic orbit

u⁽ⁱ⁺¹⁾ = 1

# TODO add roman numerals to big headers!
