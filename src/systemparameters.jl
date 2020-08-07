"""
    BODY
"""
struct BODY
    Rₑ # equatorial radius
    Rₘ # mean radius
    M # mass
    μ # gravitational parameter
    D # sidereal rotation period (day)
    T # sidereal orbital period (year)
end


"""
    CR3BP

A Circular Restricted Three-Body Problem System
"""
struct CR3BP
    μ₁
    μ₂
    μ
    d # {km} average distance between two primaries
    R1 # {km} Radius of primary BODY
    R2 # {km} Radius of Secondary Body
end

# Values from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/mar097.2100-2600.txt
# Phobos
μ₂ = 7.087546066894452E-04
# Deimos
μ₂ = 9.615569648120313E-05
# Mars
μ₂ = 4.282837362069909E+04
# System
μ₂ = 4.282837442560939E+04

# Values from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/sat427l.txt
# Mimas
μ₂ = 2.503617062809250E+00
# Enceladus
μ₂ = 7.210497553340731E+00
# Tethys
μ₂ = 4.121405263872402E+01
# Dione
μ₂ = 7.311617801921636E+01
# Rhea
μ₂ = 1.539409077211430E+02
# Titan
μ₂ = 8.978137369591670E+03
# Hyperion
μ₂ = 3.704182596063880E-01
# Iapetus
μ₂ = 1.205081845217891E+02
# Phoebe
μ₂ = 5.581081743011904E-01
# Helene
μ₂ = 4.551624250415933E-04
# Telesto
μ₂ = 0.000000000000000E+00
# Calypso
μ₂ = 0.000000000000000E+00
# Polydeuces
μ₂ = 0.000000000000000E+00
# Methone
μ₂ = 0.000000000000000E+00
# Saturn
μ₂ = 3.793120615901047E+07
# System
μ₂ = 3.794058484179918E+07

# Values from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/jup310xl.txt
# Io
μ₂ = 5.959924010272514E+03
# Europa
μ₂ = 3.202739815114734E+03
# Ganymede
μ₂ = 9.887819980080976E+03
# Callisto
μ₂ = 7.179304867611079E+03
# Amalthea
μ₂ = 1.487604677404272E-01
# Jupiter
μ₂ = 1.266865341960128E+08
# System
μ₂ = 1.267127641334463E+08

# Values from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/ura111.xl.txt
# Ariel
μ₂ = 8.346344431770477E+01
# Umbriel
μ₂ = 8.509338094489388E+01
# Titania
μ₂ = 2.269437003741248E+02
# Oberon
μ₂ = 2.053234302535623E+02
# Miranda
μ₂ = 4.319516899232100E+00
# Uranus
μ₂ = 5.793951322279009E+06
# System
μ₂ = 5.794556465751799E+06

# Values from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/nep081xl.txt
# Triton
μ₂ = 1.427598140725034E+03
# Naiad
μ₂ = 8.028519741022843E-03
# Thalassa
μ₂ = 2.358873197992170E-02
# Despina
μ₂ = 1.179477866876596E-01
# Galatea
μ₂ = 1.905015637145090E-01
# Larissa
μ₂ = 2.551116091803392E-01
# Proteus
μ₂ = 2.588564729289845E+00
# Neptune
μ₂ = 6.835099502439672E+06
# System
μ₂ = 6.836527100580397E+06

# Values from ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/plu055l.txt
# Charon
μ₂ = 1.062509269522026E+02
# Nix
μ₂ = 2.150552267969335E-03
# Hydra
μ₂ = 3.663917107480563E-03
# Kerberos
μ₂ = 4.540734312735987E-04
# Styx
μ₂ = 2.000000000000000E-20
# Pluto
μ₂ = 8.693390780381926E+02
# System
μ₂ = 9.755962735332020E+02
