# User Guide

# Example

A heaving and pitching GAW1 airfoil.

```@example
ENV["GKSwstype"] = "100"
using VSFlow
N = 100
dt = 2e-2
T = 13.34
animate = false

motion_args = (1, 25, 90, 0, .3)
motion = heavepitch(motion_args...)
shape = gaw1

airfoil = Profile(id = "hp-gaw1",
                  profileshape = shape,
                  x0 = motion(0)[1],
                  ẋ0 = motion(0)[2],
                  dt = dt,
                  T = T,
                  N = N)
profilerun(airfoil, motion, animate)

plotaero(airfoil.history)
```
