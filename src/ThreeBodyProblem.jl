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

include("conversions.jl")
export  rot2inert!,
        r2latlon

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
        findrP



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
        computeC

include("plot.jl")
export  circle,
        sphere

include("special.jl")
export  invariantManifolds

end # module
