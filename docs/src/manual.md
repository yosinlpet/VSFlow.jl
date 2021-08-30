# User Guide

## Time marching
The default numerical scheme used for time marching is the standard
[Runge-Kutta](https://en.wikipedia.org/wiki/Runge–Kutta_methods) (RK) method
up to the 4th order.

The order can be adapted when running the simulation using
`profilerun((args)...; simorder=4, kwargs...)`.

## Time derivatives
The default numerical scheme used for time derivatives is the standard
[backward difference](https://en.wikipedia.org/wiki/Finite_difference_coefficient#Backward_finite_difference)
method up to the 4th order.

The order can be adapted when running the simulation using
`profilerun((args)...; bdorder=4, kwargs...)`.

## Custom shape
Any airfoil with a wedge angle can be used.
```
my_shape(varargs...) = x -> begin
    top = my_f(x, varargs...)
    bot = my_g(x, varargs...)
    return [top, bot]
end
```
By default, `my_shape(varargs...)(0:dx:1)` is called at the airfoil
initialization. `dx` is computed to match the number of panels `N`.
The shape function has to be zero at `x=0` and `x=1`.
The trailing edge is placed at `x=1`.

## Custom body motion
The body motion has to match the structure of the following template:
```
my_motion(varargs...) = t -> begin
    x = my_position_vector(t, varargs...) #x = [-X, Y, \alpha]
    v = my_velocity_vector(t, varargs...)
    a = my_acceleration_vector(t, varargs...)
return [x, v, a]
end
```

# Examples

## A heaving and pitching airfoil
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
