using ThreeBodyProblem
using Test

sys = System(ThreeBodyProblem.JUPITER,ThreeBodyProblem.CALLISTO)
@test typeof(sys.prim) == Body
@test typeof(sys.sec) == Body
@test typeof(sys) == System

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
newmoon = Body(m, R, a, T)
@test typeof(newmoon) == Body
@test newmoon.m == m
@test newmoon.R == R
@test newmoon.a == a
@test newmoon.T == T
@test newmoon.name == "NewPlanet"
@test newmoon.color == :blue

sys = System(newplanet, newmoon)
@test typeof(sys.prim) == Body
@test typeof(sys.sec) == Body
@test typeof(sys) == System

sys = sun_mercury()
@test sys.prim.name == "Sun"
@test sys.sec.name == "Mercury"
@test sys.name == "Sun/Mercury"

sys = sun_venus()
@test sys.name == "Sun/Venus"
sys = sun_earth()
@test sys.name == "Sun/Earth"
sys = earth_moon()
@test sys.name == "Earth/Moon"
sys = sun_mars()
@test sys.name == "Sun/Mars"
sys = sun_jupiter()
@test sys.name == "Sun/Jupiter"
sys = jupiter_europa()
@test sys.name == "Jupiter/Europa"
sys = sun_saturn()
@test sys.name == "Sun/Saturn"
sys = saturn_enceladus()
@test sys.name == "Saturn/Enceladus"
sys = sun_uranus()
@test sys.name == "Sun/Uranus"
sys = sun_neptune()
@test sys.name == "Sun/Neptune"
sys = sun_earth_moon()
@test sys.name == "Sun/Earth/Moon"
sys = jupiter_europa_sun()
@test sys.name == "Jupiter/Europa/Sun"
sys = saturn_enceladus_sun()
@test sys.name == "Saturn/Enceladus/Sun"
