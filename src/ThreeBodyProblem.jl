module ThreeBodyProblem

# export modules # right now I don't have any modules within my module
#     Dynamics,
#     Problems,
#     Controllers

include("conversions.jl")
export S2I!

include("dynamics.jl")
export R2BPdynamics!,
        CR3BPdynamics!,
        CR3BPinert!,
        CWdynamics!

include("orbitalelements.jl")
export AP2a,
        findrP

include("plot.jl")
export plot_sphere,
        plot_earth,
        plot_circle,
        plot_rv


include("utils.jl")
export  findL1,
        findL2,
        findL3,
        findL4,
        findL5,
        findLpts,
        findUeff,
        findC,
        findR1R2,
        findr1r2

end # module
