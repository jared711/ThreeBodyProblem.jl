

"""
    invariantManifolds(SYS, rv0, T; tf=1., nPts=20, α=1e-6, reltol=1e-12, integrator=TsitPap8())

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
