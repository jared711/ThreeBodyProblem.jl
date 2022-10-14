# import ThreeBodyProblem.Constants
# include("constants.jl")
"""
    Body(m, R, a, T, name, color)
    Body(m, R, a, T, name)
    Body(m, R, a, T)

A planet, moon, or other gravitationally significant object
m::Float64      # mass {kg}
R::Float64      # mean radius {km}
a::Float64      # mean semimajor axis about parent body {km}
T::Float64      # sidereal orbital period {s}
name::String    # name of body (e.g. "Earth")
color::Symbol   # color of the body for plotting
"""
struct Body
    m::Float64      # mass {kg}
    R::Float64      # mean radius {km}
    a::Float64      # mean semimajor axis about parent body {km}
    T::Float64      # sidereal orbital period {s}
    name::String    # name of body (e.g. "Earth")
    color::Symbol   # color of the body for plotting
    # parent::Body    # parent body (e.g. SUN for EARTH)
end
Body(m::Float64, R::Float64, a::Float64, T::Float64, name::String) = Body(m, R, a, T, name, :blue)
Body(m::Float64, R::Float64, a::Float64, T::Float64) = Body(m, R, a, T, "NewPlanet", :blue)

"""
    System(prim, sec, μ₁, μ₂, μ, d, R₁, R₂, T, RUNIT, VUNIT, TUNIT, name)
    System(prim::Body, sec::Body)

A Circular Restricted Three-Body Problem System defined by primary and secondary bodies (prim and sec)
prim::Body      # Primary body
sec::Body       # Secondary body
μ₁::Float64     # {km^3/s^2} gravitational parameter of primary body
μ₂::Float64     # {km^3/s^2} gravitational parameter of secondary body
μ::Float64      # {} mass parameter
d::Float64      # {km} average distance between two primaries
R₁::Float64     # {km} Radius of primary body
R₂::Float64     # {km} Radius of secondary Body
T::Float64      # {s} sidereal orbital period
RUNIT::Float64  # {km} distance normalizing parameter
VUNIT::Float64  # {km/s} velocity normalizing parameter
TUNIT::Float64  # {s} time normalizing parameter
name::String    # name of system (e.g. "Earth/Moon")
"""
struct System
    prim::Body      # Primary body
    sec::Body       # Secondary body
    d::Float64      # {km} average distance between two primaries
    T::Float64      # {s} sidereal orbital period

    # non-unicode
    # Added these since it is unwise to require unicode operators from the user
    mu1::Float64     # {km^3/s^2} gravitational parameter of primary body
    mu2::Float64     # {km^3/s^2} gravitational parameter of secondary body
    mu::Float64      # {} mass parameter
    R1::Float64     # {km} Radius of primary Body
    R2::Float64     # {km} Radius of Secondary Body

    #unicode
    μ₁::Float64     # {km^3/s^2} gravitational parameter of primary body
    μ₂::Float64     # {km^3/s^2} gravitational parameter of secondary body
    μ::Float64      # {} mass parameter
    R₁::Float64     # {km} Radius of primary Body
    R₂::Float64     # {km} Radius of Secondary Body

    MUNIT::Float64  # {kg} mass normalizing parameter
    RUNIT::Float64  # {km} distance normalizing parameter
    TUNIT::Float64  # {s} time normalizing parameter
    VUNIT::Float64  # {km/s} velocity normalizing parameter
    AUNIT::Float64  # {km^2/s} acceleration normalizing parameter
    name::String    # name of system (e.g. "Earth/Moon")

end
System(prim::Body, sec::Body) = System(prim, sec, sec.a, sec.T,
    prim.m*G, sec.m*G, sec.m/(prim.m+sec.m), prim.R, sec.R,
    prim.m*G, sec.m*G, sec.m/(prim.m+sec.m), prim.R, sec.R,
    prim.m+sec.m, sec.a, sec.T/2π, sec.a/(sec.T/2π), sec.a/(sec.T/2π)^2, string(prim.name,"/",sec.name))

"""
    BicircularSystem(prim, sec, μ₁, μ₂, μ, d, R₁, R₂, T, RUNIT, VUNIT, TUNIT, name)
    BicircularSystem(prim::Body, sec::Body)

A BicircularSystem defined by primary, secondary, and tertiary bodies (prim, sec, and ter)
prim::Body      # Primary body
sec::Body       # Secondary body
ter::Body       # Tertiary body
μ₁::Float64     # {km^3/s^2} gravitational parameter of primary body
μ₂::Float64     # {km^3/s^2} gravitational parameter of secondary body
μ₃::Float64     # {km^3/s^2} gravitational parameter of tertiary body
μ::Float64      # {} mass parameter, μ₂/(μ₁+μ₂)
d::Float64      # {km} average distance between two primaries
R₁::Float64     # {km} Radius of primary body
R₂::Float64     # {km} Radius of secondary body
R₃::Float64     # {km} Radius of tertiary body
T::Float64      # {s} sidereal orbital period
RUNIT::Float64  # {km} distance normalizing parameter
VUNIT::Float64  # {km/s} velocity normalizing parameter
TUNIT::Float64  # {s} time normalizing parameter
name::String    # name of system (e.g. "Earth/Moon/Sun")
"""
struct BicircularSystem
    prim::Body      # Primary body
    sec::Body       # Secondary body
    ter::Body       # Tertiary body

    # non-unicode
    mu1::Float64    # {km^3/s^2} gravitational parameter of primary body
    mu2::Float64    # {km^3/s^2} gravitational parameter of secondary body
    mu3::Float64    # {km^3/s^2} gravitational parameter of tertiary body
    mu::Float64     # {} mass parameter, μ₂/(μ₁+μ₂)
    m3::Float64     # {} normalized mass of tertiary body
    n3::Float64     # {} normalized mean motion of tertiary body
    R1::Float64     # {km} Radius of primary Body
    R2::Float64     # {km} Radius of Secondary Body
    R3::Float64     # {km} Radius of tertiary body
    T3::Float64      # {s} sidereal orbital period of ter about prim/sec barycenter (or just about sec)

    # unicode
    μ₁::Float64     # {km^3/s^2} gravitational parameter of primary body
    μ₂::Float64     # {km^3/s^2} gravitational parameter of secondary body
    μ₃::Float64     # {km^3/s^2} gravitational parameter of tertiary body
    μ::Float64      # {} mass parameter, μ₂/(μ₁+μ₂)
    m₃::Float64     # {} normalized mass of tertiary body
    n₃::Float64     # {} normalized mean motion of tertiary body
    R₁::Float64     # {km} Radius of primary Body
    R₂::Float64     # {km} Radius of Secondary Body
    R₃::Float64     # {km} Radius of tertiary body
    T₃::Float64      # {s} sidereal orbital period of ter about prim/sec barycenter (or just about sec)

    T::Float64      # {s} sidereal orbital period of prim and sec about each other
    d::Float64      # {km} average distance between two primaries
    RUNIT::Float64  # {km} distance normalizing parameter, equivalent to d
    TUNIT::Float64  # {s} time normalizing parameter
    VUNIT::Float64  # {km/s} velocity normalizing parameter
    AUNIT::Float64  # {km^2/s} acceleration normalizing parameter
    name::String    # name of system (e.g. "Earth/Moon/Sun")
end
BicircularSystem(prim::Body, sec::Body, ter::Body) = BicircularSystem(prim, sec, ter,
    prim.m*G, sec.m*G, ter.m*G, sec.m/(prim.m+sec.m), ter.m/(prim.m+sec.m), sec.T/prim.T, prim.R, sec.R, ter.R, ter.T,
    prim.m*G, sec.m*G, ter.m*G, sec.m/(prim.m+sec.m), ter.m/(prim.m+sec.m), sec.T/prim.T, prim.R, sec.R, ter.R, ter.T,
    sec.T, sec.a, sec.a, sec.T/2π, sec.a/(sec.T/2π), sec.a/(sec.T/2π)^2, string(prim.name,"/",sec.name,"/",ter.name))


### Setting Parameters for Various Bodies and Systems
# Values from wikipedia
SUN         = Body(1.98847e30,  695700,     0.0,            0.0,            "Sun",  :yellow) # The sun does not orbit a central body

## Planets
# Values from:
# m, R, T - https://ssd.jpl.nasa.gov/?planet_phys_par,
# a - https://ssd.jpl.nasa.gov/txt/p_elem_t2.txt
MERCURY     = Body(0.330114e24,             2439.4,     0.38709843*AU,  0.2408467*JY,   "Mercury",  :grey)
VENUS       = Body(4.86747e24,              6051.8,     0.72332102*AU,  0.61519726*JY,  "Venus",    :orange)
EARTH       = Body(5.97237e24,              6371.0084,  1.00000018*AU,  1.0000174*JY,   "Earth",    :blue)
MARS        = Body(0.641712e24,             3389.50,    1.52371243*AU,  1.8808476*JY,   "Mars",     :red)
JUPITER     = Body(1898.187e24,             69911,      5.20248019*AU,  11.862615*JY,   "Jupiter",  :red)
SATURN      = Body(568.3174e24,             58232,      9.54149883*AU,  29.447498*JY,   "Saturn",   :yellow4)
URANUS      = Body(86.8127e24,              25362,      19.18797948*AU, 84.016846*JY,   "Uranus",   :cyan)
NEPTUNE     = Body(102.4126e24,             24622,      30.06952752*AU, 164.79132*JY,   "Neptune",  :blue4)
PLUTO       = Body(0.013030e24,             1188.3,     39.48686035*AU, 247.92065*JY,   "Pluto",    :grey4)

## Satellites
# Values for R, a, and T for all satellites from https://ssd.jpl.nasa.gov/?sat_phys_par, https://ssd.jpl.nasa.gov/?sat_elem
MOON        = Body(4902.801/G,              1737.5,     384400,         27.322*JD,      "Moon",     :grey)
# Values for m from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/mar097.2100-2600.txt,
PHOBOS      = Body(7.087546066894452E-04/G, 11.1,       9376.,          0.3189*JD,      "Phobos",   :grey)
DEIMOS      = Body(9.615569648120313E-05/G, 6.2,        23458.,         1.2624*JD,      "Deimos",   :grey)
# Values for m from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/jup310xl.txt
IO          = Body(5.959924010272514E+03/G, 1821.6,     421800.,        1.769*JD,       "Io",       :grey)
EUROPA      = Body(3.202739815114734E+03/G, 1560.8,     671100.,        3.551*JD,       "Europa",   :grey)
GANYMEDE    = Body(9.887819980080976E+03/G, 2631.2,     1070400.,       7.155*JD,       "Ganymede", :grey)
CALLISTO    = Body(7.179304867611079E+03/G, 2410.3,     1882700.,       16.69*JD,       "Callisto", :grey)
AMALTHEA    = Body(1.487604677404272E-01/G, 83.45,      181400.,        0.498*JD,       "Amalthea", :grey)
# Values for m from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/sat427l.txt
MIMAS       = Body(2.503617062809250E+00/G, 198.20,     185539.,        0.942*JD,       "Mimas",    :grey)
ENCELADUS   = Body(7.210497553340731E+00/G, 252.10,     238042.,        1.370*JD,       "Enceladus",:grey)
TETHYS      = Body(4.121405263872402E+01/G, 533.00,     294672.,        1.888*JD,       "Tethys",   :grey)
DIONE       = Body(7.311617801921636E+01/G, 561.70,     377415.,        2.737*JD,       "Dione",    :grey)
RHEA        = Body(1.539409077211430E+02/G, 764.30,     527068.,        4.518*JD,       "Rhea",     :grey)
TITAN       = Body(8.978137369591670E+03/G, 2574.73,    1221865.,       15.95*JD,       "Titan",    :grey)
HYPERION    = Body(3.704182596063880E-01/G, 135.00,     1500933.,       21.28*JD,       "Hyperion", :grey)
IAPETUS     = Body(1.205081845217891E+02/G, 735.60,     3560854.,       79.33*JD,       "Iapetus",  :grey)
PHOEBE      = Body(5.581081743011904E-01/G, 106.50,     12947918.,      548.02*JD,      "Phoebe",   :grey)
# Values for m from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/ura111.xl.txt
ARIEL       = Body(8.346344431770477E+01/G, 578.9,      190900.,        2.520*JD,       "Ariel",    :grey)
UMBRIEL     = Body(8.509338094489388E+01/G, 584.7,      266000.,        4.144*JD,       "Umbriel",  :grey)
TITANIA     = Body(2.269437003741248E+02/G, 788.9,      436300.,        8.706*JD,       "Titania",  :grey)
OBERON      = Body(2.053234302535623E+02/G, 761.4,      583500.,        13.46*JD,       "Oberon",   :grey)
MIRANDA     = Body(4.319516899232100E+00/G, 235.8,      129900.,        1.413*JD,       "Miranda",  :grey)
# Values for m from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/nep081xl.txt
TRITON      = Body(1.427598140725034E+03/G, 1353.4,     354759.,        5.877*JD,       "Triton",   :grey)
NEREID      = Body(2.060000000000000E+03/G, 170,        5513818.,       360.13*JD,      "Nereid",   :grey)
PROTEUS     = Body(2.588564729289845E+00/G, 210,        117646.,        1.122*JD,       "Proteus",  :grey)
# Values for m from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/plu055l.txt
CHARON      = Body(1.062509269522026E+02/G, 603.6,      19591.,         6.387*JD,       "Charon",   :grey)
NIX         = Body(2.150552267969335E-03/G, 23.0,       48671.,         24.85*JD,       "Nix",      :grey)
HYDRA       = Body(3.663917107480563E-03/G, 30.5,       48671.,         38.20*JD,       "Hydra",    :grey)
KERBEROS    = Body(4.540734312735987E-04/G, 14.0,       64698.,         32.17*JD,       "Kerberos", :grey)
STYX        = Body(2.000000000000000E-20/G, 10.0,       57729.,         20.16*JD,       "Styx",     :grey)


"""
    sun_mercury()
"""
sun_mercury() = System(SUN, MERCURY)
"""
    sun_venus()
"""
sun_venus() = System(SUN, VENUS)
"""
    sun_earth()
"""
sun_earth() = System(SUN, EARTH)
"""
    earth_moon()
"""
earth_moon() = System(EARTH, MOON)
"""
    sun_mars()
"""
sun_mars() = System(SUN, MARS)
"""
    sun_jupiter()
"""
sun_jupiter() = System(SUN, JUPITER)
"""
    sun_saturn()
"""
sun_saturn() = System(SUN, SATURN)
"""
    sun_uranus()
"""
sun_uranus() = System(SUN, URANUS)
"""
    sun_neptune()
"""
sun_neptune() = System(SUN, NEPTUNE)
"""
    jupiter_europa()
"""
jupiter_europa() = System(JUPITER, EUROPA)
"""
    saturn_enceladus()
"""
saturn_enceladus() = System(SATURN, ENCELADUS)

"""
    earth_moon_sun()
"""
earth_moon_sun() = BicircularSystem(EARTH, MOON, SUN)

"""
    sun_earth_moon()
"""
sun_earth_moon() = BicircularSystem(SUN, EARTH, MOON)

"""
    jupiter_europa_sun()
"""
jupiter_europa_sun() = BicircularSystem(JUPITER, EUROPA, SUN)

"""
    saturn_enceladus_sun()
"""
saturn_enceladus_sun() = BicircularSystem(SATURN, ENCELADUS, SUN)


## Extra values from ephemeris files
# MARS        = Body(4.282837362069909E+04/G,
# SYSTEM      = Body(4.282837442560939E+04/G,
#
# JUPITER     = Body(1.266865341960128E+08/G,
# SYSTEM      = Body(1.267127641334463E+08/G,
#
# SATURN      = Body(3.793120615901047E+07/G,
# SYSTEM      = Body(3.794058484179918E+07/G,
# HELENE      = Body(4.551624250415933E-04/G,
# TELESTO     = Body(0.000000000000000E+00/G,
# CALYPSO     = Body(0.000000000000000E+00/G,
# POLYDEUCES  = Body(0.000000000000000E+00/G,
# METHONE     = Body(0.000000000000000E+00/G,
#
# URANUS      = Body(5.793951322279009E+06/G,
# SYSTEM      = Body(5.794556465751799E+06/G,
#
# NEPTUNE     = Body(6.835099502439672E+06/G,
# SYSTEM      = Body(6.836527100580397E+06/G,
# NAIAD       = Body(8.028519741022843E-03/G,
# THALASSA    = Body(2.358873197992170E-02/G,
# DESPINA     = Body(1.179477866876596E-01/G,
# GALATEA     = Body(1.905015637145090E-01/G,
# LARISSA     = Body(2.551116091803392E-01/G,
#
# PLUTO       = Body(8.693390780381926E+02/G,
# SYSTEM      = Body(9.755962735332020E+02/G,
