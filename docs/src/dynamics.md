# [Dynamics](@id dynamics)
```@meta
CurrentModule = ThreeBodyProblem
```
## Overview
The dynamics functions give the time derivative of a state in specific systems and are of the form ẋ = f(x). Corresponding in-place methods have names that end with !.

## Functions


```@docs
R2BPdynamics
R2BPdynamics!
CR3BPdynamics
CR3BPdynamics!
CR3BPinert
CR3BPinert!
CWdynamics
CWdynamics!
BCPdynamics
BCPdynamics!
```
