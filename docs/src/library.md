```@contents
Pages = ["library.md"]
Depth = 2
```

# Public API

## Airfoil geometries
```@docs
circle
ellipse
gaw1
naca00
```

## Airfoil motion
```@docs
uniform
heavepitch
circularmotion
```

## Running the simulation
```@docs
Profile(;id, profile::Function, N, position, dt, T, δ, ϵ, lumpargs)
profilerun
```

## Post processing
```@docs
History
getimpulse(p::Profile)
getboundpotential
plotaero
plotcps
plotϕs
```

# Private API

## Vortex Elements
```@docs
VSFlow.VortexPoint
VSFlow.ConstantPanel
VSFlow.LinearPanel
```

## Airfoil geometries
```@docs
VSFlow.initiator
VSFlow.carlingfish
VSFlow.splinefish
VSFlow.swimmerthickness
VSFlow.getbodymass
VSFlow.getinertia
VSFlow.getcenterofmass
VSFlow.getcenterofmassvelocity
VSFlow.correctvelocity
```

## Maths
```@docs
VSFlow.rotate
VSFlow.simps
VSFlow.trapz
VSFlow.backwarddifference
VSFlow.centereddifference
VSFlow.checkangle
VSFlow.gausslegendretriangle
VSFlow.get∫fdV
VSFlow.cubicspline
```
