using RecipesBase

function plot_circle(r=1,c=[0,0,0];color="blue",label="",npts=100)
    u = range(0,stop=2π,length=npts)
    x = r*cos.(u) .+ c[1]
    y = r*sin.(u) .+ c[2]
    return x,y
end

@recipe function f(sys::System)
    name1, name2 = split(sys.name, '/')
    L1,L2,L3,L4,L5 = findLpts(sys.μ)

    legend := true

    @series begin
        markershape --> :x
        seriestype --> :scatter
        label := "Lagrange Points"
        x = [L1[1], L2[1], L3[1], L4[1], L5[1]]
        y = [L1[2], L2[2], L3[2], L4[2], L5[2]]
        x,y
    end

    @series begin
        label := name1
        x,y = plot_circle(sys.R₁/sys.RUNIT, [-sys.μ,0,0])
    end

    @series begin
        label := name2
        x,y = plot_circle(sys.R₂/sys.RUNIT, [1-sys.μ,0,0], label="Secondary")
    end
end

# function plot_sphere(r=1,c=[0,0,0],col='b',n=100)
#     u = range(0,stop=2π,length=n)
#     v = range(0,stop=π,length=n)
#     x = cos.(u) * sin.(v)';
#     y = sin.(u) * sin.(v)';
#     z = ones(n) * cos.(v)';
#     const_color = cgrad( [ RGB{Float64}(0.,0.,1.) for _ in 1:2 ] )
#     s = surface!(r*x,r*y,r*z,colorbar=false,color=const_color,width=1)
#     return s
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
