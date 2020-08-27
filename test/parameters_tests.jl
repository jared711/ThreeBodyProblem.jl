using ThreeBodyProblem
using Test

PRIM, SEC, SYS = set_system(ThreeBodyProblem.JUPITER,ThreeBodyProblem.CALLISTO)
@test typeof(PRIM) == Body
@test typeof(SEC) == Body
@test typeof(SYS) == System

# Create the ice planet Hoth
m = 1e24 # {kg} mass
R = 1e6 # {km} mean radius
a = 1.5*AU # {km} mean semimajor axis
T = 0.6*JY # {sec} mean orbital period
name = "Hoth"
newplanet = Body(m, R, a, T, name)
@test typeof(newplanet) == Body
@test newplanet.m == m
@test newplanet.R == R
@test newplanet.a == a
@test newplanet.T == T
@test newplanet.name == name

# Create the forest moon Endor
m = 1e22 # {kg} mass
R = 1e5 # {km} mean radius
a = 0.2*AU # {km} mean semimajor axis
T = 2π*sqrt(a^3/(G*newplanet.m)) # {sec} mean orbital period
name = "Endor"
newmoon = Body(m, R, a, T, name)
@test typeof(newmoon) == Body
@test newmoon.m == m
@test newmoon.R == R
@test newmoon.a == a
@test newmoon.T == T
@test newmoon.name == name

PRIM, SEC, SYS = set_system(newplanet, newmoon)
@test typeof(PRIM) == Body
@test typeof(SEC) == Body
@test typeof(SYS) == System

PRIM, SEC, SYS = sun_mercury()
@test PRIM.name == "Sun"
@test SEC.name == "Mercury"
@test SYS.name == "Sun/Mercury"

PRIM, SEC, SYS = sun_venus()
PRIM, SEC, SYS = sun_earth()
PRIM, SEC, SYS = earth_moon()
PRIM, SEC, SYS = sun_mars()
PRIM, SEC, SYS = sun_jupiter()
PRIM, SEC, SYS = jupiter_europa()
PRIM, SEC, SYS = sun_saturn()
PRIM, SEC, SYS = saturn_enceladus()
PRIM, SEC, SYS = sun_uranus()
PRIM, SEC, SYS = sun_neptune()
