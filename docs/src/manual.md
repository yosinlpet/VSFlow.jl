# User Guide

# Example

A heaving and pitching NACA0013 airfoil.

```@example
ENV["GKSwstype"] = "100"
using VSFlow
N = 100
dt = 5e-2
T = 13.5
animate = false

motion_args = (1, 25, 90, 0, .3)
motion = heavepitch(motion_args...)
shape = naca00(13)

airfoil = Profile(id = "hp-naca0013",
                  profileshape = shape,
                  x0 = motion(0)[1],
                  ẋ0 = motion(0)[2],
                  dt = dt,
                  T = T,
                  N = N)
profilerun(airfoil, motion, animate)

plotaero(airfoil.history)
```
