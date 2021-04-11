module ThreeBodyProblem

include("constants.jl")
export JD, JY, AU, G

include("parameters.jl")
export  System,
        Body,
        set_system,
        sun_mercury,
        sun_venus,
        sun_earth,
        earth_moon,
        sun_mars,
        sun_jupiter,
        sun_saturn,
        sun_uranus,
        sun_neptune,
        jupiter_europa,
        saturn_enceladus

include("frames.jl")
export  rot2inert!,
        r2latlon,
        enu2ecef,
        ecef2enu,
        ecef2eci,
        eci2ecef,
        eci2sci,
        sci2eci

include("dynamics.jl")
export  R2BPdynamics!,
        R2BPdynamics,
        CR3BPdynamics!,
        CR3BPdynamics,
        CR3BPinert!,
        CR3BPinert,
        CWdynamics!,
        CWdynamics

include("orbitalelements.jl")
export  AP2a,
        cart2oe,
        oe2cart,
        nu2E,
        E2nu,
        E2M,
        M2E,
        azel2cart,
        cart2azel



include("util.jl")
export  computeR1R2,
        computer1r2,
        computeL1,
        computeL2,
        computeL3,
        computeL4,
        computeL5,
        computeLpts,
        computeUeff,
        computeC,
        rotx,
        rotxd,
        roty,
        rotyd,
        rotz,
        rotzd,
        rotlatlon,
        date2mjd,
        mjd2gmst,
        wrapto360,
        wrapto180,
        wrapto2pi,
        wraptopi,
        date2str


include("plot.jl")
export  circle,
        sphere

include("special.jl")
export  invariant_manifolds,
        differential_corrector

end # module
