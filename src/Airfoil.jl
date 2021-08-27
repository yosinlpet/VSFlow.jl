#!/usr/bin/env julia
# File              : Airfoil.jl
# Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
# Date              : 18.03.2019
# Last Modified Date: 03.03.2021

"""
    gaw1(x)

Return the `y` coordinates of a GAW1 airfoil.
"""
function gaw1(x)
	ye = @. .278536sqrt(x) - .148567x + .006397x^2 - .220980x^3 + .081084x^4
	yi = @. -.190361sqrt(x)  + .161628x - .341176x^2 + .897631x^3 - .531252x^4
	return [ye, yi]
end

"""
    naca00(ZZ)(x)

Return the `y` coordinates of a NACA00ZZ airfoil.
"""
naca00(ZZ) = x -> begin
    ye = @. .05ZZ*(.2969sqrt(x) - .1260x - .3516x^2 + .2843x^3 - .1036x^4)
    yi = -ye
    return [ye, yi]
end

"""
    circle(R)(x)

Return the `y` coordinates of a circle of radius `R`.
"""
circle(R) = x -> begin
	ye = @. sqrt(R^2 - (x-R)^2)
	yi = -ye
	return [ye, yi]
end

"""
    ellipse(B)(x)

Return the `y` coordinates of an ellipse with half height `B`.
"""
ellipse(B) = x -> begin
	ye = @. .5b*sqrt(1 - 4(x-.5)^2)
	yi = -ye
	return [ye, yi]
end

"""
    uniform(vx, vy, α̇)(t)

Return the vector `[x, ẋ, ẍ]` corresponding to a uniform straight motion in
the `[-X, Y, -θ]` direction  with velocity `[vx, vy, α̇]`.
"""
uniform(vx, vy, α̇) = t -> ([-vx*t, vy*t, -α̇*t], [-vx, vy, -α̇], [0, 0, 0])

"""
    heavepitch(h0, αmax, ψ1, ψ2, strouhal)(t)

Return the vector `[x, ẋ, ẍ]` corresponding to a heaving & pitching motion where
`x = [X, Y, α]`.
"""
heavepitch(h0, αmax, ψ1, ψ2, strouhal) = t -> begin
	ω = strouhal*pi/h0
	Y = h0*sin(ω*t + deg2rad(ψ1))
	Ẏ = ω*h0*cos(ω*t + deg2rad(ψ1))
	Ÿ = -(ω^2)*h0*sin(ω*t + deg2rad(ψ1))
	Y⃛ = -(ω^3)*h0*cos(ω*t + deg2rad(ψ1))

	αh = atan(-Ẏ)
	α̇h = -Ÿ / (1+Ẏ^2)
	α̈h = 2Ẏ*Ÿ^2/(1+Ẏ^2)^2 - Y⃛/(1 + Ẏ^2)

	α = deg2rad(αmax)*sin(ω*t + deg2rad(ψ2)) - αh
	α̇ = ω*deg2rad(αmax)*cos(ω*t + deg2rad(ψ2)) - α̇h
	α̈ = -ω^2*deg2rad(αmax)*sin(ω*t + deg2rad(ψ2)) - α̈h
	return [[-t, Y, α], [-1, Ẏ, α̇], [0, Ÿ, α̈]]
end

"""
    circularmotion(R, strouhal)(t)

Return the vector `[x, ẋ, ẍ]` corresponding to a circular motion of radius `R` where
`x = [X, Y, α]`.
"""
circularmotion(R, strouhal) = t -> begin
	ω = .5strouhal*pi/R
	ωt = ω*t

	X = R*cos(ωt)
	Ẋ = -ω*R*sin(ωt)
	Ẍ = -ω^2*R*cos(ωt)

	Y = R*sin(ωt)
	Ẏ = ω*R*cos(ωt)
	Ÿ = -ω^2*R*sin(ωt)

	return [[X, Y, -ωt+π/2], [Ẋ, Ẏ, -ω], [Ẍ, Ÿ, 0]]
end

#Constant values for the fish
const wh = .04
const sb = .04
const st = .95
const wt = .01

"""
    initiator(ξ, T)

Smooth increase from 0 to 1 in half a period `T` and its time derivative.
"""
function initiator(ξ, T)
	#= ξ >= T/2 && return [1, 0] =#
	#= return [.5 + .5sin(π*(2ξ/T - .5)), π/T*cos(π*(2ξ/T - .5))] =#
	ξ >= T && return [1, 0]
	return [.5 + .5sin(π*(ξ/T - .5)), .5π/T*cos(π*(ξ/T - .5))]
end

"""
    carlingfish(s, t, T)

Return the curvature of a curve of unit length associated with the motion of a
carling fish of period `T`. Its time derivative is computed too.
"""
function carlingfish(s, t, T)
	@assert 0 ≤ s ≤ 1 "Wrong coordinate 's'"
	a = .125
	b = .03125
	c = 1+b
	y = 2π*(s-t/T)
	κ1 = 4a/c*(π*cos(y) - π^2*(b+s)*sin(y))/(1+(a/c*(sin(y) + 2π*(b+s)*cos(y)))^2)^1.5
	κ2 = (75.3982*((1.0472a*c^2*(π*(b+s)*cos(y) + sin(y))) /((a^2*(2π*(b + s)*cos(y) + sin(y))^2)/c^2 + 1)^1.5
				   + (π*a^3*(2π*(b + s)*cos(y) + sin(y))*(cos(y) - 2π*(b + s)*sin(y))*(cos(y) - π*(b + s)*sin(y)))/((a^2*(2π*(b + s)*cos(y) + sin(y))^2)/c^2 + 1)^2.5))/(c^3*T)

	i1, i2 = initiator(t, T)
	return [κ1*i1; (κ1*i2 + κ2*i1)]
end

function anticarlingfish(s, t, T)
    -carlingfish(s, t, T)
end

function deadfish(s, t, T)
	zeros(2)
end

"""
    splinefish(s, t, T, X, K)

Return the curvature of a curve of unit length associated with the motion of a
fish following a cubic spline.
"""
function splinefish(s, t, T, X, K)
	@assert 0 ≤ s ≤ 1 "Wrong coordinate 's'"
	Ks = cubicspline(s, X, K)
	y = 2π*(t/T - 1.72s)
	κ1 = Ks*sin(y)
	κ2 = 2π/T*Ks*cos(y)

	i1, i2 = initiator(t, T)
	return [κ1*i1; κ1*i2 + κ2*i1]
end

function efficientfish(s, t, T)
	X = [0, 1/3, 2/3, 1]
	K = [3.34, 1.67, 2π, 2π]
	return splinefish(s, t, T, X, K)
end

function fastfish(s, t, T)
	X = [0, 1/3, 2/3, 1]
	K = [1.51, .48, 5.74, 2.73]
	return splinefish(s, t, T, X, K)
end

"""
    swimmerthickness(s)

Return the thickness of the swimming body along its coordinate `s`.
"""
function swimmerthickness(s)
	s < sb && return [sqrt(s*(2wh - s)), -sqrt(s*(2wh - s))]
	#= sb ≤ s < st && return [wh - (wh-wt)*((s-sb)/(st-sb))^2, -wh + (wh-wt)*((s-sb)/(st-sb))^2] =#
	sb ≤ s < st && return [wh - (wh-wt)*(s-sb)/(st-sb), -wh + (wh-wt)*(s-sb)/(st-sb)]
	return [wt*(1-s)/(1-st), -wt*(1-s)/(1-st)]
end

