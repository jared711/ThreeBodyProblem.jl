# [Frames](@id frames)
```@meta
CurrentModule = ThreeBodyProblem
```

```@contents
Pages = ["frames.md"]
```

## Overview
The most important frame conversion when dealing with the Circular Restricted Three-Body Problem (CR3BP) is between the rotating and inertial frames.

## Functions

```@docs
rot2inert
rot2inert!
inert2rot
inert2rot!
enu2ecef
ecef2enu
ecef2eci
eci2ecef
eci2sci
sci2eci
dimensionalize,
dimensionalize!,
nondimensionalize,
nondimensionalize!
```
