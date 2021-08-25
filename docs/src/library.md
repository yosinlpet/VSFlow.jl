# Public API

## Airfoil geometries
```@docs
naca00
gaw1
circle
ellipse
```

## Airfoil motion
```@docs
heavepitch
circularmotion
```

## Running the simulation
```@docs
Profile(id, profile::Function, N, position, dt, T, δ, ϵ, rng, is_comb, lump, args...)
setinitvelocity
profilerun
```

## Quantities of interest
```@docs
getimpulse(p::Profile)
getboundpotential
getcoeffimpulse!
getcoeffpotential!
getnoca!
```

# Private API

## Vortex Elements
```@docs
VortexPoint
ConstantPanel
LinearPanel
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
