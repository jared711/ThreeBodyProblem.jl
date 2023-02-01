using Plots: cgrad

"""
    circle(r=1, c=[0,0,0]; n=100)

Outputs the x and y coordinates of a sphere with radius 'r', centered at 'c', and with 'n' points in each direction.
To plot, just use 'plot(circle(), seriescolor=:blue, seriestype=:shape, fillalpha=0.5, lims=:auto)'
"""
function circle(r=1,c=[0,0,0]; n=100)
    u = range(0,stop=2π,length=n)
    x = r*cos.(u) .+ c[1]
    y = r*sin.(u) .+ c[2]
    return x,y
end

# """
#     plot_circle(r=1, c=[0,0,0]; col=:blue, n=100)

# Outputs a plot of the circle produced by circle(). Default color is blue.
# """
# function plot_circle(r=1,c=[0,0,0]; col=:blue,n=100,label="")
#     x,y = circle(r,c,n=n)
#     plot_c = plot(x,y,seriescolor=col, seriestype=:shape, fillalpha=0.5, label=label, lims=:auto)
#     return plot_c
# end


# """
#     plot_circle!(r=1, c=[0,0,0]; col=:blue, n=100)

# Overlays a plot of the circle produced by circle(). Default color is blue.
# """
# function plot_circle!(r=1,c=[0,0,0]; col=:blue,n=100,label="")
#     x,y= circle(r,c,n=n)
#     plot!(x,y,seriescolor=col, seriestype=:shape, fillalpha=0.5, label=label,lims=:auto)
#     return nothing
# end

"""
    sphere(r=1, c=[0,0,0]; n=100)

Outputs the x, y, and z coordinates of a sphere with radius 'r', centered at 'c', and with 'n' points in each direction.
To plot, just use 'surface(sphere(),colorbar=false,color=cgrad([:blue, :blue]),width=1,alpha=1).'

"""
function sphere(r=1, c=[0,0,0]; n=100)
    θ = range(0,stop=2π,length=n)
    ψ = range(0,stop=π,length=n)
    x = r*cos.(θ) * sin.(ψ)' .+ c[1];
    y = r*sin.(θ) * sin.(ψ)' .+ c[2];
    z = r*ones(n) * cos.(ψ)' .+ c[3];
    return x, y, z
end

# """
#     plot_sphere(r=1, c=[0,0,0]; col=:blue, n=100)

# Outputs a plot of the sphere produced by sphere(). Default color is blue.
# """
# function plot_sphere(r=1,c=[0,0,0]; col=:blue,n=100)
#     x,y,z = sphere(r,c,n=n)
#     const_color = cgrad([col, col])
#     plot_s = surface(x,y,z,colorbar=false,color=const_color,width=1,alpha=1)
#     return plot_s
# end

# """
#     plot_sphere!(r=1, c=[0,0,0]; col=:blue, n=100)

# Overlays a plot of the sphere produced by sphere(). Default color is blue.
# """
# function plot_sphere!(r=1, c=[0,0,0]; col=:blue,n=100)
#     x,y,z = sphere(r,c,n=n)
#     const_color = cgrad([col, col])
#     surface!(x,y,z,colorbar=false,color=const_color,width=1,alpha=1)
#     return nothing
# end


"""
    torus(r=1, c=[0,0,0]; n=100)

Outputs the x, y, and z coordinates of a torus with primary radius a (radius of the circle being rotated around the z-axis)
and secondary radius b (radius of the circle perpendicular to the z axis), centered at 'c', and with 'n' points in each direction.
To plot the torus, just use 'surface(torus(),colorbar=false,color=cgrad([:blue, :blue]),lims=(-10,10)).'
"""
function torus(a=5, b=10, c=[0,0,0]; n=50) 
    ψ = ϕ = range(0, 2π, length=n+1) # ψ is longitude ϕ is latitude
    x = @. (b + a * cos(ψ')) * cos(ϕ) + c[1] 
    y = @. (b + a * cos(ψ')) * sin(ϕ) + c[2]   
    z = @.      a * sin(ψ')  * one(ϕ) + c[3]
    return x,y,z
    # surface(x, y, z; lims, colorbar=false)
end


# """
#     plot_torus(r=1, c=[0,0,0]; col=:blue, n=100)

# Outputs a plot of the torus produced by torus(). Default color is blue.
# """
# function plot_torus(a=5, b=10; col=:blue, n=50)
#     x,y,z = torus(a,b,n=n)
#     const_color = cgrad([col, col])
#     plot_t = surface(x,y,z,colorbar=false,color=const_color,lims=(-10,10))
#     return plot_t
# end

# """
#     plot_torus(r=1, c=[0,0,0]; col=:blue, n=100)

# Overlays a plot of the torus produced by torus(). Default color is blue.
# """
# function plot_torus!(a=5, b=10; col=:blue, n=50)
#     x,y,z = torus(a,b,n=n)
#     const_color = cgrad([col, col])
#     surface!(x,y,z,colorbar=false,color=const_color,lims=(-10,10))
#     return nothing
# end

"""
    trajectory_2D(rv)

Outputs the x and y coordinates of the trajectory 'rv', which is a vector of vectors.
"""
function trajectory_2D(rv::Vector{Vector{Float64}})
    N = length(rv)
    x = [rv[i][1] for i = 1:N]
    y = [rv[i][2] for i = 1:N]
    return x,y
end

"""
    trajectory_2D(rv)

Outputs the x, y, and z coordinates of the trajectory 'rv', which is a vector of vectors.
"""
function trajectory_3D(rv::Vector{Vector{Float64}})
    N = length(rv)
    x = [rv[i][1] for i = 1:N]
    y = [rv[i][2] for i = 1:N]
    z = [rv[i][3] for i = 1:N]
    return x,y,z
end

"""
    trajectory_2D(rv)

Outputs the ẋ and ẏ components of the velocity of the trajectory 'rv', which is a vector of vectors.
"""
function velocity_2D(rv)
    N = length(rv)
    if length(rv[1]) == 4
        ẋ = [rv[i][3] for i = 1:N]
        ẏ = [rv[i][4] for i = 1:N]
    elseif length(rv[1]) == 6
        ẋ = [rv[i][4] for i = 1:N] # If the trajectory is 3D, we still only want the x and y components of the velocity
        ẏ = [rv[i][5] for i = 1:N]
    end
    return ẋ,ẏ
end

function velocity_3D(rv)
    N = length(rv)
    ẋ = [rv[i][4] for i = 1:N]
    ẏ = [rv[i][5] for i = 1:N]
    ż = [rv[i][6] for i = 1:N]
    return ẋ,ẏ,ż
end

"""
    Recipe for plotting systems
"""
@recipe function f(sys::System; planar=true, prim=true, sec=true, Lpts=true, scaled=false, center=[])
    name1, name2 = split(sys.name, '/')
    color1 = sys.prim.color
    color2 = sys.sec.color
    L1,L2,L3,L4,L5 = computeLpts(sys.μ)

    legend --> true # colon equals := forces the operator even if it exists
    legend --> :topleft
    xguide --> "x [NON]"
    yguide --> "y [NON]"
    # box_aspect := [1,1,1]

    if planar
        aspect_ratio --> :equal
        if Lpts
            @series begin
                markershape --> :x # arrow --> sets the operator only when it doesn't exist
                seriestype --> :scatter
                markercolor --> :black
                label --> "Lagrange Points"
                x = [L1[1], L2[1], L3[1], L4[1], L5[1]]
                y = [L1[2], L2[2], L3[2], L4[2], L5[2]]
                x,y
            end
        end

        if prim
            @series begin
                label --> name1
                seriescolor := color1
                seriestype --> :shape
                fillalpha --> 0.5
                if isempty(center)
                    center = [-sys.μ,0,0]
                end
                if scaled
                    x,y = circle(0.1, center)
                else
                    x,y = circle(sys.R₁/sys.RUNIT, center)
                end
            end
        end

        if sec
            @series begin
                label --> name2
                seriescolor := color2
                seriestype --> :shape
                fillalpha --> 0.5
                if isempty(center)
                    center = [1-sys.μ,0,0]
                end
                if scaled
                    x,y = circle(0.025, center)
                else
                    x,y = circle(sys.R₂/sys.RUNIT, center)
                end
            end
        end
    else
        zguide := "Z [NON]"
        # error("Non-planar not implemented yet.")
        xlims --> (-1.,1.)
        ylims --> (-1.,1.)
        zlims --> (-1.,1.)
        if Lpts
            @series begin
                markershape --> :x # arrow --> sets the operator only when it doesn't exist
                seriestype --> :scatter
                markercolor --> :black
                label --> "Lagrange Points"
                x = [L1[1], L2[1], L3[1], L4[1], L5[1]]
                y = [L1[2], L2[2], L3[2], L4[2], L5[2]]
                z = zeros(5)
                x,y,z
            end
        end

        if prim
            @series begin
                label --> name1
                seriescolor := cgrad( [color1, color1] )
                colorbar --> false
                seriestype --> :surface
                if scaled
                    x,y,z = sphere(0.1, [-sys.μ,0,0])
                else
                    x,y,z = sphere(sys.R₁/sys.RUNIT, [-sys.μ,0,0])
                end
            end
        end

        if sec
            @series begin
                label --> name2
                seriescolor := cgrad( [color2, color2] )
                colorbar --> false
                seriestype --> :surface
                if scaled
                    x,y,z = sphere(0.025, [1-sys.μ,0,0])
                else
                    x,y,z = sphere(sys.R₂/sys.RUNIT, [1-sys.μ,0,0])
                end
            end
        end
    end
    
end

# I'm leaning away from these functions, but I'll keep them for now. I prefer just using
# the native plot function and inputting the system as an argument along with parameters
# that determine whether to plot the desired components.

# function plot_Lpts(sys::System; planar=true, scaled=false)
#     plot(sys,prim=false,sec=false,planar=planar,scaled=scaled,lims=:auto)
# end

# function plot_Lpts!(sys::System; planar=true, scaled=false)
#     plot!(sys,prim=false,sec=false,planar=planar,scaled=scaled,lims=:auto)
# end

# function plot_prim(sys::System; planar=true, scaled=false)
#     plot(sys,Lpts=false,sec=false,planar=planar,scaled=scaled,lims=:auto)
# end

# function plot_prim!(sys::System; planar=true, scaled=false)
#     plot!(sys,Lpts=false,sec=false,planar=planar,scaled=scaled,lims=:auto)
# end

# function plot_sec(sys::System; planar=true, scaled=false)
#     plot(sys,Lpts=false,prim=false,planar=planar,scaled=scaled,lims=:auto)
# end

# function plot_sec!(sys::System; planar=true, scaled=false)
#     plot!(sys,Lpts=false,prim=false,planar=planar,scaled=scaled,lims=:auto)
# end

"""
    Recipe for plotting bodies

This is a recipe for plotting bodies. It is called by the plot function. It always the body at the origin.
If planar is true, it plots the body as a circle. If planar is false, it plots the body as a sphere.
The radius of the body is the true dimensional radius in km. Pass in 1/sys.RUNIT as the scalar parameter
to scale the body to the CR3BP.
"""
@recipe function f(body::Body; planar=true, legend=true, center=[0,0,0], scalar=1)
    legend := legend
    legend := :topleft
    
    if planar
        aspect_ratio --> :equal
        @series begin
            label := body.name
            seriescolor := body.color
            seriestype --> :shape
            fillalpha --> 0.5
            x,y = circle(body.R*scalar, center)
        end
    else
        @series begin
            label := body.name
            seriescolor := cgrad( [ RGB{Float64}(0.,0.,1.) for _ in 1:2 ] )
            colorbar --> false
            seriestype --> :surface
            x,y,z = sphere(body.R*scalar, center)
        end
    end
end

"""
    Recipe for plotting trajectories
"""
@recipe function f(rv::Vector{Vector{T}} where T<:Real; label="trajectory", color=:black, planar=false, vel = false)

    @series begin
        label := label
        seriescolor := color
        if vel
            if length(rv[1]) == 4 || planar
                x,y = velocity_2D(rv)
            elseif length(rv[1]) == 6
                x,y,z = velocity_3D(rv)
            else
                error("Velocity must be 2D or 3D")
            end
        else
            if length(rv[1]) == 2 || length(rv[1]) == 4 || planar
                x,y = trajectory_2D(rv)
            elseif length(rv[1]) == 3 || length(rv[1]) == 6
                x,y,z = trajectory_3D(rv)
            else
                error("Trajectory must be 2D or 3D")
            end
        end
    end
end

"""
    Produces the appropriate limits for a plot to show the secondary body
"""
function seczoom(sys::System)
    c = [1-sys.μ, 0, 0] # center of secondary body, center of plot
    xlims = c[1] .+ []
    return nothing
end

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
