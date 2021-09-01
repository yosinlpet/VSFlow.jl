# User Guide

## Time marching
The default numerical scheme used for time marching is the standard 4th order
[Runge-Kutta](https://en.wikipedia.org/wiki/Runge–Kutta_methods) (RK4) method.
The order can be reduced to `1` when running the simulation using
`profilerun((args)...; is4thorder=false, kwargs...)`.
Note that the RK1 is the same as the Forward Euler scheme.

## Time derivatives
The default numerical scheme used for time derivatives is the standard
[backward difference](https://en.wikipedia.org/wiki/Finite_difference_coefficient#Backward_finite_difference)
method up to the 4th order.
The order can be adapted when running the simulation using
`profilerun((args)...; bdorder=4, kwargs...)`.

## Vortex Kernel
By default, the flow model uses a [Low-Algebraic]() vortex regularization kernel
with cut-off width of `1e-2`.
Nonetheless, the user can use a [Gaussian]() kernel if they so desire.
These two parameters can be specified in the `Profile` declaration.
`Profile(isgaussian=false, δ=1e-2, kwargs...)`.


## Custom shape
Any airfoil with a wedge angle can be used.
The template is the following:
```
my_shape(varargs...) = x -> begin
    top = my_f(x, varargs...)
    bot = my_g(x, varargs...)
    return [top, bot]
end
```
By default, `my_shape(varargs...)(0:dx:1)` is called at the airfoil
initialization. `dx` is computed to match the number of panels `N`.
The shape functions `my_f` and `my_g` must belong to the following class of
functions:

```math
f: [0, 1] \to \mathbb{R},\quad f(0) = f(1) = 0.
```

The trailing edge is placed at `x=1`.

## Airfoil discretization
The panel distribution is based on the flat projection of uniformly
spaced points on the unit circle centered at `x=0.5`.
Therefore, panels closest to both edges of the airfoil are smaller than panels
in the middle of the airfoil where its surface has a larger curvature radius.

# Examples

## A heaving and pitching airfoil
A heaving and pitching NACA0013 airfoil.
The heaving motion has an amplitude of 1, the pitching motion has an amplitued
of 25º, the Strouhal number is .3 and the
phase shift between heaving and pitching is 90º.

```@example
ENV["GKSwstype"] = "100"
using VSFlow
N = 100
dt = 5e-2
T = 13.5

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
profilerun(airfoil, motion)

## We could plot Cp and ϕ distribution at time t=12.
#plotcps(airfoil.history, 12.)
#plotϕs(airfoil.history, 12.)

plotaero(airfoil.history)
```
