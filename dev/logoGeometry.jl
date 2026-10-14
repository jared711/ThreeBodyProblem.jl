ϕ = (1+sqrt(5))/2   # golden ratio

d = [1, ϕ, ϕ^2, ϕ^3]
r = d/2
A = π*r.^2

a = (d[3] + d[2])/2
b = (d[3] + d[1])/2
c = (d[2] + d[1])/2

C = acosd((a^2 + b^2 - c^2)/(2*a*b))

x = [b*sind(C), 0, 0, 0]
y = [(d[3] - d[4])/2 + b*cosd(C), (d[4] - d[2])/2,(d[3] - d[4])/2, 0]

x_balance = zeros(360)
y_balance = zeros(360)
xy_balance = zeros(360)
x'*A
y'*A
x'*A + y'*A
for θ = 1:360 # {deg}
    R = [cosd(θ) -sind(θ); sind(θ) cosd(θ)]

    xynew = R*[x';y']
    xnew = xynew[1,:]
    ynew = xynew[2,:]

    x_balance[θ] = xnew'*A
    y_balance[θ] = ynew'*A
    xy_balance[θ] = xnew'*A + ynew'*A
end

using PyPlot
plot(x_balance,label="x_balance")
plot(y_balance,label="y_balance")
plot(xy_balance,label="xy_balance")

grid("on")
legend()

θˣ = findall(x->x==minimum(abs.(x_balance),2),x_balance)

θᶜᶜ = 323.5 # {deg} counterclockwise
θᶜ = 360 - θᶜᶜ # {deg} clockwise

x_final = zeros(4)
y_final = zeros(4)
for θ = θᶜᶜ # {deg}
    R = [cosd(θ) -sind(θ); sind(θ) cosd(θ)]
    global x_final, y_final
    xynew = R*[x';y']
    x_final = xynew[1,:]
    y_final = xynew[2,:]
end

x_final
y_final
