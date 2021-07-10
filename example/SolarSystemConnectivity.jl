using ThreeBodyProblem
using DifferentialEquations
using Plots

T = 1

sys = sun_earth()
Lpts = computeLpts(sys)
L1 = Lpts[1]
tf = 1*JY/sys.TUNIT
rv₀ = vcat(L1, zeros(3))
@time Wsp, Wsn, Wup, Wun, Φₜ = invariant_manifolds(sys, rv₀, T, tf)
# Juno.@enter invariant_manifolds(sys, rv₀, T, tf)
plot(sys)
plot!(Wup,vars=(1,2),label="Wu+",linecolor=:red)
plot!(Wun,vars=(1,2),label="Wu-",linecolor=:magenta)
plot!(Wsp,vars=(1,2),label="Ws+",linecolor=:blue)
plot!(Wsn,vars=(1,2),label="Ws-",linecolor=:cyan)
plot!(aspect_ratio=1,ylims=[-0.04,0.04],xlims=[0.93,1.06],legend=:outerright,flip=false)

W78 = []
C = []
a = []
e = []
t = []
systems = (sun_mercury(), sun_venus(), sun_earth(), sun_mars(), sun_jupiter(), sun_saturn(), sun_uranus(), sun_neptune())
for sys in systems
    L1, L2 = computeLpts(sys)
    C1 = []
    a1 = []
    e1 = []
    t1 = []
    for Lpt = (L1, L2)
        tf = 100*JY/sys.TUNIT # integrate for 100 years
        rv₀ = vcat(Lpt, zeros(3))
        Wsp, Wsn, Wup, Wun, Φₜ = invariant_manifolds(sys, rv₀, T, tf)
        push!(W78,Wsp,Wsn,Wup,Wun)

        C2 = []
        a2 = []
        e2 = []
        t2 = []
        for W0 = (Wsp, Wsn, Wup, Wun)
            push!(C2, [computeC(W0.u[j], sys) for j = 1:length(W0.u)])
            rv_inert = [rot2inert(W0.u[j], W0.t[j], sys.μ, origin=:prim) for j = 1:length(W0.u)]
            rv_dim = [dimensionalize(rv_inert[i], sys) for i = 1:length(rv_inert)]
            push!(a2, [cart2oe(rv_dim[i], sys.μ₁)[1] for i = 1:length(rv_dim)])
            push!(e2, [cart2oe(rv_dim[i], sys.μ₁)[2] for i = 1:length(rv_dim)])
            push!(t2, [W0.t[j]*sys.TUNIT/JY for j = 1:length(W0.u)])
        end
        push!(C1, C2)
        push!(a1, a2)
        push!(e1, e2)
        push!(t1, t2)
    end
    push!(C,C1)
    push!(a,a1)
    push!(e,e1)
    push!(t,t1)

        # plot(sys)
        # plot!(Wsp,vars=(1,2),label="Ws+",linecolor=:blue)
        # plot!(Wsn,vars=(1,2),label="Ws-",linecolor=:cyan)
        # plot!(Wup,vars=(1,2),label="Wu+",linecolor=:red)
        # plot!(Wun,vars=(1,2),label="Wu-",linecolor=:magenta)
        # plot!(aspect_ratio=1,legend=:outerright,flip=false)
end

# Mercury L1 SP
plot(t[1][1][1], C[1][1][1],label="Wₛ⁺",xlabel="time [years]", ylabel="Jacobi Constant, C", title="Mercury L1 RK78",color=:blue)
# Mercury L2 SP
plot(t[1][2][1], C[1][2][1],label="Wₛ⁺",xlabel="time [years]", ylabel="Jacobi Constant, C", title="Mercury L2 RK78",color=:blue)
# Mercury L1 SN
plot(t[1][1][2], C[1][1][2],label="Wₛ⁻",xlabel="time [years]", ylabel="Jacobi Constant, C", title="Mercury L1 RK78",color=:cyan)
# Mercury L2 SN
plot(t[1][2][2], C[1][2][2],label="Wₛ⁻",xlabel="time [years]", ylabel="Jacobi Constant, C", title="Mercury L2 RK78",color=:cyan)
plot(sys)
plot!(W78[1],vars=(1,2),label="Ws+",linecolor=:blue)
plot!(W78[2],vars=(1,2),label="Ws-",linecolor=:cyan)
plot!(W78[5],vars=(1,2),label="Ws+",linecolor=:blue)
plot!(W78[6],vars=(1,2),label="Ws-",linecolor=:cyan)

scatter(a[1][1][1]/AU,e[1][1][1], color=:blue, label="Mercury", markerstrokewidth=0)
scatter!(a[1][1][2]/AU,e[1][1][2], color=:cyan, markerstrokewidth=0)
scatter!(a[1][1][3]/AU,e[1][1][3], color=:red, markerstrokewidth=0)
scatter!(a[1][1][4]/AU,e[1][1][4], color=:magenta, markerstrokewidth=0)

scatter(a[1][2][1]/AU,e[1][2][1], color=:blue, label="Mercury", markerstrokewidth=0)
scatter!(a[1][2][2]/AU,e[1][2][2], color=:cyan, markerstrokewidth=0)
scatter!(a[1][2][3]/AU,e[1][2][3], color=:red, markerstrokewidth=0)
scatter!(a[1][2][4]/AU,e[1][2][4], color=:magenta, markerstrokewidth=0)
#L1 SN UP
#L2 SP UN

scatter(a[1][1][1]/AU,e[1][1][1], color=:grey, label="Mercury", markerstrokewidth=0)
scatter!(a[2][1][1]/AU,e[2][1][1], color=:orange, label="Venus", markerstrokewidth=0)
scatter!(a[3][1][1]/AU,e[3][1][1], color=:blue, label="Earth", markerstrokewidth=0)
scatter!(a[4][1][1]/AU,e[4][1][1], color=:red, label="Mars", markerstrokewidth=0)
xlabel!("a [AU]")
ylabel!("e")
scatter!(a[5][1][1]/AU,e[5][1][1], color=:blue, label="Mercury", markerstrokewidth=0)
scatter!(a[6][1][1]/AU,e[6][1][1], color=:blue, label="Mercury", markerstrokewidth=0)
scatter!(a[7][1][1]/AU,e[7][1][1], color=:blue, label="Mercury", markerstrokewidth=0)
scatter!(a[8][1][1]/AU,e[8][1][1], color=:blue, label="Mercury", markerstrokewidth=0)



for (i, sys) in enumerate(systems)
    display(i)
    display(sys)
    plot(sys)
    plot!(W[(i-1)*8+1],vars=(1,2),label="Wu+",linecolor=:red)
    plot!(W[(i-1)*8+2],vars=(1,2),label="Wu-",linecolor=:magenta)
    plot!(W[(i-1)*8+3],vars=(1,2),label="Ws+",linecolor=:blue)
    plot!(W[(i-1)*8+4],vars=(1,2),label="Ws-",linecolor=:cyan)
    plot!(W[(i-1)*8+5],vars=(1,2),label="Wu+",linecolor=:red)
    plot!(W[(i-1)*8+6],vars=(1,2),label="Wu-",linecolor=:magenta)
    plot!(W[(i-1)*8+7],vars=(1,2),label="Ws+",linecolor=:blue)
    plot!(W[(i-1)*8+8],vars=(1,2),label="Ws-",linecolor=:cyan)
    plot!(aspect_ratio=1,legend=:outerright,flip=false)
end


plot(sys)
plot!(W[17],vars=(1,2),label="Wu+",linecolor=:red)
plot!(W[18],vars=(1,2),label="Wu-",linecolor=:magenta)
plot!(W[19],vars=(1,2),label="Ws+",linecolor=:blue)
plot!(W[20],vars=(1,2),label="Ws-",linecolor=:cyan)
plot!(aspect_ratio=1,legend=:outerright,flip=false)
