```@meta
CurrentModule = VSFlow
```

# VSFlow

*An inviscid planar flow model with vortex shedding for external aerodynamics*

## Installation

PotentialFlow can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run

```julia
pkg> add VSFlow
```
The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).

## Basic Usage

Let's simulate an impulsively started NACA0012 with a 10º pitch angle.

```@example startingnaca0012
using VSFlow

# Simulation parameters
N = 100          # number of vortex panels on the body surface
δ = 1e-2         # blob kernel cut-off width
T = 5            # horizon time
dt = 5e-2        # time step
animate = true   # generate animation
filewrite = true # save results to file

# Body geometry
shape = naca00(12)
initial_position = [0, 0, deg2rad(10)]
airfoil_ID = 01

# Body motion
function steady(t)
	return [-t, 0, 0], [-1, 0, 0], [0, 0, 0]
end
motion_args = ()

airfoil = Profile(id = airfoil_ID, profileshape = shape, N = N, position = initial_position, dt = dt, T = T, δ = δ, ϵ = ϵ)
setinitvelocity(airfoil, steady(-dt/2, motion_args...)[2]...)
@time profilerun(airfoil, steady, motion_args, filewrite, animate)
```

Documentation for [VSFlow](https://github.com/yosinlpet/VSFlow.jl/dev).
