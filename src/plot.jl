using RecipesBase
using DifferentialEquations

function circle(r=1,c=[0,0,0];color="blue",label="",npts=100)
    u = range(0,stop=2π,length=npts)
    x = r*cos.(u) .+ c[1]
    y = r*sin.(u) .+ c[2]
    return x,y
end

function sphere(r=1,c=[0,0,0],col='b',n=100)
    θ = range(0,stop=2π,length=n)
    ψ = range(0,stop=π,length=n)
    x = r*cos.(θ) * sin.(ψ)';
    y = r*sin.(θ) * sin.(ψ)';
    z = r*ones(n) * cos.(ψ)';
    # const_color = cgrad( [ RGB{Float64}(0.,0.,1.) for _ in 1:2 ] )
    # s = surface!(r*x,r*y,r*z,colorbar=false,color=const_color,width=1)
    return x,y,z
end

function trajectory(rv)
    if length(rv[1]) == 2 || length(rv[1]) == 4
        x = [rv[i][1] for i = 1:length(rv)]
        y = [rv[i][2] for i = 1:length(rv)]
        return x,y
    elseif length(rv[1]) == 3 || length(rv[1]) == 6
        x = [rv[i][1] for i = 1:length(rv)]
        y = [rv[i][2] for i = 1:length(rv)]
        z = [rv[i][3] for i = 1:length(rv)]
        return x,y,z
    else
        error("rv should be 2xN, 3xN, 4xN, or 6xN")
    end
end


@recipe function f(sys::System; prim=true, sec=true, Lpts=true, scaled=false)
    name1, name2 = split(sys.name, '/')
    color1 = sys.prim.color
    color2 = sys.sec.color
    L1,L2,L3,L4,L5 = computeLpts(sys.μ)

    legend := true
    legend := :topleft
    aspect_ratio --> :equal

    if Lpts
        @series begin
            markershape --> :x
            seriestype --> :scatter
            markercolor --> :black
            label := "Lagrange Points"
            x = [L1[1], L2[1], L3[1], L4[1], L5[1]]
            y = [L1[2], L2[2], L3[2], L4[2], L5[2]]
            x,y
        end
    end

    if prim
        @series begin
            label := name1
            seriescolor := color1
            seriestype --> :shape
            fillalpha --> 0.5
            if scaled
                x,y = circle(0.1, [-sys.μ,0,0])
            else
                x,y = circle(sys.R₁/sys.RUNIT, [-sys.μ,0,0])
            end
        end
    end

    if sec
        @series begin
            label := name2
            seriescolor := color2
            seriestype --> :shape
            fillalpha --> 0.5
            if scaled
                x,y = circle(0.025, [1-sys.μ,0,0])
            else
                x,y = circle(sys.R₂/sys.RUNIT, [1-sys.μ,0,0])
            end
        end
    end
end

# function plot_Lpts(sys::System)
#     plot(sys,prim=false,sec=false)
# end
#
# function plot_Lpts!(sys::System)
#     plot!(sys,prim=false,sec=false)
# end
#
# function plot_prim(sys::System)
#     plot(sys,Lpts=false,sec=false)
# end
#
# function plot_prim!(sys::System)
#     plot!(sys,Lpts=false,sec=false)
# end
#
# function plot_sec(sys::System)
#     plot(sys,Lpts=false,prim=false)
# end
#
# function plot_sec!(sys::System)
#     plot!(sys,Lpts=false,prim=false)
# end


@recipe function f(body::Body; legend=true, center=[0,0,0], scalar=1)
    legend := legend
    legend := :topleft
    aspect_ratio --> :equal

    @series begin
        label := body.name
        seriescolor := body.color
        seriestype --> :shape
        fillalpha --> 0.5
        x,y = circle(body.R*scalar, center)
    end
end

@recipe function f(rv::Vector{Vector{T}} where T<:Real; label="trajectory", color=:black, planar=false)

    @series begin
        label := label
        seriescolor := color
        if length(rv[1]) == 2 || length(rv[1]) == 4 || planar
            x,y = trajectory(rv)
        elseif length(rv[1]) == 3 || length(rv[1]) == 6
            x,y,z = trajectory(rv)
        end
    end
end

# function seczoom()
#     plot!(xlim=[0.995,1.005],ylim=[-0.005,0.005])
#     return nothing
# end

# @recipe function f(odevec::Vector{ODESolution}; label="trajectory", color=:black)
#
#     for i = 1:length(odevec)
#         @series begin
#             label := label
#             seriescolor := color
#             if length(odevec[i][1]) == 2 || length(odevec[i][1]) == 4
#                 x,y = trajectory(odevec[i])
#             elseif length(odevec[i][1]) == 3 || length(odevec[i][1]) == 6
#                 x,y,z = trajectory(odevec[i])
#             end
#         end
#     end
# end

# @recipe function f(ode::ODESolution; label="trajectory", color=:black)
#
#     @series begin
#         label := label
#         seriescolor := color
#         if length(ode[1]) == 2 || length(ode[1]) == 4
#             x,y = trajectory(ode)
#         elseif length(ode[1]) == 3 || length(ode[1]) == 6
#             x,y,z = trajectory(ode)
#         end
#     end
# end

#
# function plot_earth(;CR3BP=true,col='b',n=100)
#     CR3BP ? r = earth.r/moon.a      : r = earth.r
#     CR3BP ? c = [-earth_moon.μ,0,0] : c = [0,0,0]
#     plot_sphere(r,c,col,n)
# end
#
#
# function plot_rv(rv;dim=3)
#     r,v = rv[1:3], rv[4:6]
#     if dim == 2
#         p = plot!(r[1,:],r[2,:])
#     elseif dim == 3
#         p = plot!(r[1,:],r[2,:],r[3,:])
#     else
#         error("dim should be 2 or 3")
#     end
#     return p
# end
