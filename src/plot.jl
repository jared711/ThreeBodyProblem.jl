using Plots
plotly()

function plot_sphere(r=1,c=[0,0,0],col='b',n=100)
    u = range(0,stop=2π,length=n)
    v = range(0,stop=π,length=n)
    x = cos.(u) * sin.(v)';
    y = sin.(u) * sin.(v)';
    z = ones(n) * cos.(v)';
    const_color = cgrad( [ RGB{Float64}(0.,0.,1.) for _ in 1:2 ] )
    s = surface!(r*x,r*y,r*z,colorbar=false,color=const_color,width=1)
    return s
end

function plot_earth(;CR3BP=true,col='b',n=100)
    CR3BP ? r = earth.r/moon.a      : r = earth.r
    CR3BP ? c = [-earth_moon.μ,0,0] : c = [0,0,0]
    plot_sphere(r,c,col,n)
end

function plot_circle(r=1,c=[0,0,0];color="blue",label="",npts=100)
    u = range(0,stop=2π,length=npts)
    x = r*cos.(u)
    y = r*sin.(u)
    c = plot!(x .+ c[1], y .+ c[2] ,label=label, color=color)
    return c
end

function plot_rv(rv;dim=3)
    r,v = rv[1:3], rv[4:6]
    if dim == 2
        p = plot!(r[1,:],r[2,:])
    elseif dim == 3
        p = plot!(r[1,:],r[2,:],r[3,:])
    else
        error("dim should be 2 or 3")
    end
    return p
end
