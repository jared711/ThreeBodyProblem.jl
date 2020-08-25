module ThreeBodyProblem

include("constants.jl")
export JD, JY, AU, G

include("conversions.jl")
export  rot2inert

include("dynamics.jl")
export  R2BPdynamics!,
        CR3BPdynamics!,
        CR3BPinert!,
        CWdynamics!

include("orbitalelements.jl")
export  AP2a,
        findrP


include("parameters.jl")
export  System,
        Body,
        PRIM,
        SEC,
        SYS,
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

include("util.jl")
export  findR1R2,
        findr1r2,
        findL1,
        findL2,
        findL3,
        findL4,
        findL5,
        findLpts,
        findUeff,
        findC

include("plot.jl")
# export plot_sphere,
#         plot_earth,
#         plot_circle,
#         plot_rv

end # module
