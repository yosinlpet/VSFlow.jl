#!/usr/bin/env julia
# File              : Utils.jl
# Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
# Date              : 14.03.2019
# Last Modified Date: 08.12.2020

#Constant values for Gauss quadrature in triangles
const wgauss = [4.40470137440027E-04 2.73540341982354E-03 6.06288503309780E-03 7.32940618020361E-03 4.35981888065413E-03 2.54894897110758E-03 1.58294584351310E-02 3.50852038616632E-02 4.24144130415724E-02 2.52297599892837E-02 4.35986399676418E-03 2.70755855459964E-02 6.00116828031788E-02 7.25479695591883E-02 4.31543838150997E-02 2.80482254237839E-03 1.74184820315194E-02 3.86071952834637E-02 4.66721394461879E-02 2.77624230060559E-02 3.69367248237134E-04 2.29384093975323E-03 5.08418385425357E-03 6.14625683304730E-03 3.65603514489964E-03]

const xgauss = [5.65222820508010E-03 5.65222820508010E-03 5.65222820508010E-03 5.65222820508010E-03 5.65222820508010E-03 7.34303717426523E-02 7.34303717426523E-02 7.34303717426523E-02 7.34303717426523E-02 7.34303717426523E-02 2.84957404462558E-01 2.84957404462558E-01 2.84957404462558E-01 2.84957404462558E-01 2.84957404462558E-01 6.19482264084778E-01 6.19482264084778E-01 6.19482264084778E-01 6.19482264084778E-01 6.19482264084778E-01 9.15758083004698E-01 9.15758083004698E-01 9.15758083004698E-01 9.15758083004698E-01 9.15758083004698E-01]

const ygauss = [5.62028052139780E-03 7.30153265243790E-02 2.83346760183808E-01 6.15980808959171E-01 9.10582009338909E-01 5.23718298680677E-03 6.80383522483882E-02 2.64032876322051E-01 5.73993451145053E-01 8.48513626543325E-01 4.04158392633041E-03 5.25058436021453E-02 2.03756682104520E-01 4.42956206000591E-01 6.54806036556072E-01 2.15077307947324E-03 2.79415588029271E-02 1.08431346378371E-01 2.35723988569175E-01 3.48462192391012E-01 4.76154539290863E-04 6.18591528127868E-03 2.40053580139315E-02 5.21863734710916E-02 7.71452164162586E-02]

"""
    rotate(X, Y, θ)

Rotate frame by angle `θ` in the clockwise direction.
"""
function rotate(X, Y, θ)
	x = X*cos(θ) - Y*sin(θ)
	y = Y*cos(θ) + X*sin(θ)
	return [x, y]
end

"""
    norm2(x)

Squared norm of `x`.
"""
norm2(x...) = sum(x.^2)

"""
    simps(f, x)

Integrate function `f` using Simpson's rule.
"""

function simps(f, x)
	N = length(x) - 1
	N < 2 && return trapz(f, x)
    h = diff(x)

    result = 0.
    for i in 2:2:N
        hph = h[i] + h[i - 1]
        result += f[i]*(h[i]^3 + h[i-1]^3 + 3h[i]*h[i - 1]*hph)/(6h[i]*h[i-1])
        result += f[i-1]*(2h[i-1]^3 - h[i]^3 + 3h[i]*h[i-1]^2)/(6h[i-1]*hph)
        result += f[i+1]*(2h[i]^3 - h[i-1]^3 + 3h[i-1]*h[i]^2)/(6h[i]*hph)
	end

    if (N + 1) % 2 == 0
        result += f[N]*(2h[N-1]^2 + 3h[N-2]*h[N-1])/(6(h[N-2] + h[N-1]))
        result += f[N-1]*(h[N-1]^2 + 3h[N-1]*h[N-2])/(6h[N-2])
        result -= f[N-2]*h[N-1]^3/(6*h[N-2]*(h[N-2] + h[N-1]))
	end
	return result
end

"""
    trapz(f, x)

Integrate function `f` using trapeze rule.
"""
function trapz(y, x)
	n = length(y)
	length(y) == length(x) || error("Dimension mismatch between x and f(x).")

	r = zero(zero(eltype(x))*zero(eltype(y)))
	for i in 2:n
		r += (x[i] - x[i-1])*(y[i] + y[i-1])
	end
	.5r
end

"""
    backwarddifference(x, old_x, dt)

Compute the backward difference of nth order according to the length of `old_x`.
"""
function backwarddifference(x, old_x, dt)
	n = length(old_x)
	0 < n < 5 || error("5+th order Backward Difference undefined.")
	A = [-1 1 0 0 0;
		 -1.5 2 -.5 0 0;
		 -11/6 3 -1.5 1/3 0;
		 -25/12 4 -3 4/3 -.25]
	-sum(A[n, 1:n+1].*[x, old_x...])/dt
end

"""
    customfourthorder(x, old_x, dt)

Compute the 4th order time derivative at t using
	t-1/2
	t-3/2
	t-5/2
	t-7/2.
"""
function customfourthorder(x, old_x, dt)
	A = [.5 1.5 2.5 3.5;
		 1 9 25 49;
		 1 27 125 343;
		 1 81 625 2401]
	a = A\[1 0 0 0]'
	return sum([sum(a), -a...].*[x, old_x...])/dt
end

"""
    centereddifference(neighbours, dx)

Centered differences. Order depends on length of neighbours.
"""
function centereddifference(neighbours, dx)
	n = length(neighbours)
	if n == 2
		return [-.5 .5]*neighbours /dx
	elseif n == 4
		return [1/12 -2/3 2/3 -1/12]*neighbours /dx
	else
		error("6+th order Backward Difference undefined.")
	end
end

"""
    checkangle(θ)

Assert the continuity of panel angles by subtracting 2π if necessary.
"""
function checkangle(θ)
	θ > .5π && return θ-2π
	return θ
end

"""
    heavepitch(t, h0, αmax, ψ1, ψ2, strouhal)

Return the vector `[x, ẋ, ẍ]` corresponding to a heaving & pitching motion where
`x = [X, Y, α]`.
"""
function heavepitch(t, h0, αmax, ψ1, ψ2, strouhal)
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
    circularmotion(t, R, strouhal)

Return the vector `[x, ẋ, ẍ]` corresponding to a circular motion of radius `R` where
`x = [X, Y, α]`.
"""
function circularmotion(t, R, strouhal)
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

"""
    gausslegendretriangle(f, fx, fy)

Computes the integral of function `f(x,y)`
on the standard triangle abc:
a (0; 0)
b (1; 0)
c (0; 1)
using 25 points per dimension.
fx, fy are the mapping from the standard triangle to the real one.
Coefficients are extracted from Rathod et al (2004)
"Gauss Legendre quadrature over a triangle"
"""
function gausslegendretriangle(f, fx, fy)
	x = fx.(xgauss, ygauss)
	y = fy.(xgauss, ygauss)
	return (f.(x, y)*wgauss')[1]
end

"""
    get∫fdV(XY, f)

Computes the integral of `f` on the triangle created by the
set of points `XY`, `Y` using a Gauss-Legendre quadrature.
`f` can return a vector.
"""
function get∫fdV(XY1, XY2, XY3, f)
	J =	(XY3[1] - XY2[1])*(XY1[2] - XY2[2]) - (XY1[1] - XY2[1])*(XY3[2] - XY2[2])
	x(u, v)	= XY2[1] + u*(XY3[1] - XY2[1]) + v*(XY1[1] - XY2[1])
	y(u, v)	= XY2[2] + u*(XY3[2] - XY2[2]) + v*(XY1[2] - XY2[2])
	Is = gausslegendretriangle(f, x, y)
	return abs(J).*Is
end

"""
    cubicspline(z, X, U)

Computes the cubic spline at point `z ∈ [0, 1]`, build from
interpolation of points `(X, U)` as well as its derivative.
"""
function cubicspline(z, X, U)
	@assert length(X) == length(U) "X and U must have the same length!"
	@assert 0 ≤ z ≤ 1 "z must be on the interval [0, 1]"

	n = length(U)-1
	dv = 4ones(n+1)
	dv[1] = 2
	dv[end] = 2
	ev = ones(n)
	A = SymTridiagonal(dv, ev)
	x = map((a,b)->3(b-a), U, circshift(U, -2))
	x[1] = 3(U[2] - U[1])
	x[end] = 3(U[end] - U[end-1])
	b = (A\x)
	a = U[1:end-1]
	c = map((w, x, y, z)->3(x-w) - 2y - z, U, circshift(U, -1), b, circshift(b, -1))
	d = map((w, x, y, z)->2(w-x) + y + z, U, circshift(U, -1), b, circshift(b, -1))
	b = b[1:end-1]
	c = c[1:end-1]
	d = d[1:end-1]

	x = (z*ones(n+1) - X) ./ (circshift(X, -1) - X)
	i = findfirst(0 .≤ x[1:end] .≤ 1)
	return a[i] + b[i]*x[i] + c[i]*x[i]^2 + d[i]*x[i]^3

end

"""
    getbodymass(XYm, XYt, XYb)

Integrate the constant and uniform unit density on the body.
Integral computed with 25-point Gauss quadrature.
"""
function getbodymass(XYm, XYt, XYb)
	I12 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->1), XYm[1:end-1], XYm[2:end], XYt[2:end]))
	I22 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->1), XYm[1:end-1], XYm[2:end], XYb[2:end]))
	I32 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->1), XYm[1:end-1], XYt[1:end-1], XYt[2:end]))
	I42 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->1), XYm[1:end-1], XYb[1:end-1], XYb[2:end]))
	return (I12 + I22 + I32 + I42)
end

"""
    getinertia(XYm, XYt, XYb)

Return the moment of inertia of the body through 25-point Gauss quatradure on triangles.
"""
function getinertia(XYm, XYt, XYb)
	I1 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->norm2(x, y)), XYm[1:end-1], XYm[2:end], XYt[2:end]))
	I2 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->norm2(x, y)), XYm[1:end-1], XYm[2:end], XYb[2:end]))
	I3 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->norm2(x, y)), XYm[1:end-1], XYt[1:end-1], XYt[2:end]))
	I4 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->norm2(x, y)), XYm[1:end-1], XYb[1:end-1], XYb[2:end]))
	return I1 + I2 + I3 + I4
end

"""
    getcenterofmass(xm, ym, xb, yb, xt, yt)

Compute the center of mass of the body.
Integrals computed with 25-point Gauss quadrature.
"""
function getcenterofmass(xm, ym, xb, yb, xt, yt)
	XYm = map((x,y)->[x, y], xm, ym)
	XYt = map((x,y)->[x, y], xt, yt)
	XYb = map((x,y)->[x, y], xb, yb)
	I1 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->[x, y]), XYm[1:end-1], XYm[2:end], XYt[2:end]))
	I2 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->[x, y]), XYm[1:end-1], XYm[2:end], XYb[2:end]))
	I3 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->[x, y]), XYm[1:end-1], XYt[1:end-1], XYt[2:end]))
	I4 = sum(map((XY1, XY2, XY3)->get∫fdV(XY1, XY2, XY3, (x, y)->[x, y]), XYm[1:end-1], XYb[1:end-1], XYb[2:end]))
	mass = getbodymass(XYm, XYt, XYb)
	return  [((I1 .+ I2 .+ I3 .+ I4) ./ mass)..., mass]
end

"""
    getcenterofmassvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)

Compute the velocity of the center of mass of the body.
Integrals computed with bi-linear assumption on quadrilaterals.
"""
function getcenterofmassvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)
	XYm = map((x,y)->[x, y], xm, ym)
	XYt = map((x,y)->[x, y], xt, yt)
	XYb = map((x,y)->[x, y], xb, yb)
	diags1 = (XYb[2:end] - XYt[1:end-1])
	diags2 = (XYb[1:end-1] - XYt[2:end])
	Areas = .5abs.(map((xy1, xy2)->xy1[1]*xy2[2] - xy1[2]*xy2[1], diags1, diags2))
	I1 = sum(.25Areas.*(ut[1:end-1] .+ ut[2:end] .+ ub[1:end-1] .+ ub[2:end]))
	I2 = sum(.25Areas.*(vt[1:end-1] .+ vt[2:end] .+ vb[1:end-1] .+ vb[2:end]))
	return  [I1; I2] ./ getbodymass(XYm, XYt, XYb)
end

"""
    correctvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)

Compute the angular velocity (ratio between angular impulse and
scalar inertia moment).
Integrals computed with bi-linear assumption on quadrilaterals,
and 25-point Gauss quadrature.
"""
function correctvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)
	XYm = map((x,y)->[x, y], xm, ym)
	XYt = map((x,y)->[x, y], xt, yt)
	XYb = map((x,y)->[x, y], xb, yb)
	diags1 = (XYb[2:end] - XYt[1:end-1])
	diags2 = (XYb[1:end-1] - XYt[2:end])
	Areas = .5abs.(map((xy1, xy2)->xy1[1]*xy2[2] - xy1[2]*xy2[1], diags1, diags2))
	height1 = map((x, y, u, v)->x*v - u*y, xt, yt, ut, vt)
	height3 = map((x, y, u, v)->x*v - u*y, xb, yb, ub, vb)
	I = sum(.25Areas.*(height1[1:end-1] .+ height1[2:end] .+ height3[1:end-1] .+ height3[2:end]))
	J = getinertia(XYm, XYt, XYb)
	return I/J, J
end
