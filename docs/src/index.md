```@meta
CurrentModule = VSFlow
```

# VSFlow

*An inviscid planar flow model with vortex shedding for external aerodynamics*
Documentation for [VSFlow](https://github.com/yosinlpet/VSFlow.jl).

## Installation

PotentialFlow can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run

```julia
pkg> add VSFlow
```
The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).

## Basic Usage

Let's simulate an impulsively started NACA0012 with a 10º pitch angle.

```@example startingnaca0012
ENV["GKSwstype"] = "100"
using VSFlow

# Simulation parameters
N = 100          # number of vortex panels on the body surface
T = 5            # horizon time
dt = 5e-2        # time step
animate = true   # generate animation

# Body geometry and motion
airfoil_ID = "My-naca0012"
shape = naca00(12.)
motion = uniform(1., 0., 0.)
initial_position = [0., 0., deg2rad(10.)]
initial_velocity = motion(0.)[2]

# Building the `Profile` object
airfoil = Profile(id = airfoil_ID,
                    profileshape = shape,
                    x0 = initial_position,
                    ẋ0 = initial_velocity,
                    N = N,
                    dt = dt,
                    T = T)

# Run the simulation
simulate(airfoil, motion, animate)

# Show aerodynamic coeffs
plotaero(airfoil.history)
```

A history of most useful data in the flow is saved at each time step in a data
structure called `history`. These data can be saved to a file or rendered with
[Plots.jl](http://docs.juliaplots.org/latest/).

The animation is saved as a `.gif` file. Its name contains all simulation
paramaters.

![Impulsively started NACA0012](assets/My-naca0012_np100_dt005_T5_dv001_eps001.gif)

