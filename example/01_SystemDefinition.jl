@time using ThreeBodyProblem
@time using Plots

### Defining New Bodies ###

## We use the Body() constructor to create a new planet and a new moon

# Create the ice planet Hoth
m = 1e24 # {kg} mass
R = 1e6 # {km} mean radius
a = 1.5*AU # {km} mean semimajor axis
T = 0.6*JY # {sec} mean orbital period
name = "Hoth"
color = :blue
newplanet = Body(m, R, a, T, name, color)

# Create the forest moon Endor
m = 1e22 # {kg} mass
R = 1e6 # {km} mean radiusasdf
a = 0.2*AU # {km} mean semimajor axis
T = 2π*sqrt(a^3/(G*newplanet.m)) # {sec} mean orbital period
name = "Endor"
color = :green
newmoon = Body(m, R, a, T, name, color)

### Note that we can use the constants AU, JY, JD, and G to define our bodies

println("AU = ",AU) # km in an astronomical unit
println("JY = ",JY) # seconds in a Julian Year
println("JD = ",JD) # seconds in a Julian day
println("G = ",G)  # universal gravitational constant

## Defining a New System

### We use the System() constructor with newplanet and newmoon as our primary and secondary bodies.

sys = System(newplanet, newmoon)

### We can access system parameters with the dot operator, including the body structs in the system

println("System period is ",sys.T, " seconds")
println("The Primary Body is ", sys.prim.name)

### In the CR3BP, the dynamics of a system are completely defined by its mass parameter, μ.

sys.μ

### While it looks pretty to use greek letters μ and other unicode operators like superscripts r¹ and subscripts r₁, 
# they may cause trouble with plotting and in other situations. So there are always non-unicode alternatives, like sys.mu instead of sys.μ

sys.mu

## Plotting a System

### We can plot a System object using the plot() command

plot(sys)

# Notice that the plot is made in the rotating frame, with non-dimensionalized coordinates

###  We can also plot individual bodies

plot(newmoon)

# Note that bodies plotted alone will be centered at the origin and be plotted in dimensional coordinates

### We can plot just the bodies in the CR3BP by setting Lpts to false, or vice versa

plot(sys, Lpts=false)

plot(sys, prim=false, sec=false)

### There are pre-programmed variables and functions for many common systems.

println("EARTH = ",EARTH)
println("SUN = ",SUN)
println("MOON = ",MOON)
println("JUPITER = ",JUPITER)
println("EUROPA = ",EUROPA)
println("STYX = ",STYX)
sys = System(SUN,MERCURY)

plot(jupiter_europa())

# When you plot these real systems, it's interesting to see the relative size between the primary and secondary

plot(saturn_enceladus())

### You can set the 'scaled' argument to 'true' to make the bodies easier to see, but note that they are not their true size

plot(saturn_enceladus(),scaled=true)