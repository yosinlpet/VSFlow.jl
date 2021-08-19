#!/usr/bin/env julia
# File              : VortexElement.jl
# Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
# Date              : 13.03.2019
# Last Modified Date: 04.03.2021
include("Utils.jl")
using LinearAlgebra
#= using FMMLIB2D =#

const pi2_inv = .5/pi
const sq3_inv = 1/sqrt(3.)
const w14 = (18 - sqrt(30))/36
const w23 = (18 + sqrt(30))/36
const X14 = sqrt((15 + 2sqrt(30))/35)
const X23 = sqrt((15 - 2sqrt(30))/35)

mutable struct VortexPoint
	X::Float64
	Y::Float64
	Γ::Float64
	dv::Vector{Float64}
end

VortexPoint(X, Y, Γ) = VortexPoint(X, Y, Γ, [0., 0.])

function blobkernel(X, Y, vpX, vpY, δ)
	## Gaussian regularization kernel
	x2 = norm2(X - vpX, Y - vpY)
    abs(x2) < 1e-12 && return zeros(2)
	coeff = pi2_inv*(1-exp(-x2/δ^2))/x2

	## Low order algebraic i.e. Blob
	#= coeff = pi2_inv/norm2(X - vpX, Y - vpY, δ) =#
	return [(vpY - Y)coeff, (X - vpX)coeff]

end

function getinducedvelocity(vp::VortexPoint, X, Y, δ)
	return vp.Γ*blobkernel(X, Y, vp.X, vp.Y, δ)
end

"""
vortexcloudvelocity
	This function uses anti-symmetry to obtain the velocities induced by
	the vortex cloud on itself
"""
function getcloudvelocities(vortex_points, δ)
	Nv = length(vortex_points)
	U = zeros(Nv, Nv)
	V = zeros(Nv, Nv)
	for i = 1:Nv
		@inbounds vi = vortex_points[i]
		@inbounds for j = i+1:Nv
			@inbounds vj = vortex_points[j]
	 		@inbounds U[i, j], V[i, j] = blobkernel(vi.X, vi.Y, vj.X, vj.Y, δ)
			@inbounds U[j, i], V[j, i] = -U[i, j], -V[i, j]
		end
	end
	Γs = map(vp->vp.Γ, vortex_points)
	return map((u, v)->[u, v], U*Γs, V*Γs)
end

"""getvelocitygradient(vp::VortexPoint, X, Y, θ)
	Computes the gradient of induced velocity at position 'X, Y'
    due to a change of position of 'vp'

    Input :
        -vp : sensitive point vortex
        -X  : X position where velocity is computed
        -Y  : Y position where velocity is computed
        -θ  : rotation angle
        -δ  : Gaussian kernel cut-off width
"""
function getvelocitygradient(vp::VortexPoint, X, Y, θ, δ)
	coeff = -pi2_inv/norm2(X - vp.X, Y - vp.Y, δ)^2
	bx = (((X - vp.X)^2 - norm2(Y - vp.Y, δ))cos(θ) - 2(vp.X - X)*(vp.Y - Y)sin(θ))coeff
	by = ((norm2(X - vp.X, δ) - (Y - vp.Y)^2)sin(θ) + 2(vp.X - X)*(vp.Y - Y)cos(θ))coeff
	return [bx, by]
end

abstract type Panel
end

function BuildPanel(T, X1, Y1, X2, Y2, profID)
	Xc = .5 * (X1+X2)
	Yc = .5 * (Y1+Y2)
	b = .5norm([X2-X1, Y2-Y1])
	θ = atan(Y2 - Y1, X1 - X2)
	return T(X1, Y1, X2, Y2, Xc, Yc, b, θ, profID)
end

function inpanelframe(p::Panel, X, Y)
	X, Y = rotate(X - p.Xc, Y - p.Yc, π+p.θ)
end

function inreferenceframe(p::Panel, x, y)
	X, Y = rotate(x, y, -(π+p.θ))
	return X + p.Xc, Y + p.Yc
end

function cropedge!(p::Panel, ϵ, isend)
	p.b -= .5ϵ
	if isend
		p.X2 = p.X1 + 2(p.b * cos(pi - p.θ))
		p.Y2 = p.Y1 + 2(p.b * sin(pi - p.θ))
	else
		p.X1 = p.X2 + 2(p.b * cos(p.θ))
		p.Y1 = p.Y2 - 2(p.b * sin(p.θ))
	end
	p.Xc = .5(p.X1+p.X2)
	p.Yc = .5(p.Y1+p.Y2)
end

function setangle!(p::Panel, α, pivotlocation)
	x1 = p.X1 - pivotlocation[1]
	y1 = p.Y1 - pivotlocation[2]
	x2 = p.X2 - pivotlocation[1]
	y2 = p.Y2 - pivotlocation[2]
	p.X1, p.Y1 = pivotlocation + rotate(x1, y1, -α)
	p.X2, p.Y2 = pivotlocation + rotate(x2, y2, -α)
	p.Xc = .5(p.X1+p.X2)
	p.Yc = .5(p.Y1+p.Y2)
	p.θ = atan(p.Y2 - p.Y1, p.X1 - p.X2)
end

function translate!(p::Panel, Xt, Yt)
	p.X1 += Xt
	p.Y1 += Yt
	p.X2 += Xt
	p.Y2 += Yt
	p.Xc += Xt
	p.Yc += Yt
end


mutable struct ConstantPanel <: Panel
	X1::Float64
	Y1::Float64
	X2::Float64
	Y2::Float64
	Xc::Float64
	Yc::Float64
	b::Float64
	θ::Float64
	profID::String
end

function getparams(p::ConstantPanel, X, Y)
	x, y = inpanelframe(p, X, Y)
	u = (atan((x-p.b)/y) - atan((x+p.b)/y))pi2_inv
	v = .5pi2_inv*log(norm2(x+p.b, y)/norm2(x-p.b, y))
	return u, v
end

mutable struct LinearPanel <: Panel
	X1::Float64
	Y1::Float64
	X2::Float64
	Y2::Float64
	Xc::Float64
	Yc::Float64
	b::Float64
	θ::Float64
	profID::String
end

function LinearPanel(X1, Y1, X2, Y2, Xc, Yc, b, θ, profID)
	return LinearPanel(X1, Y1, X2, Y2, Xc, Yc, b, θ, profID)
end

"""
function getparams(p::LinearPanel, X, Y)
	Obtain the parameters to multiply with γL, γR in order to get the
	velocity at point X, Y.
"""
function getparams(p::LinearPanel, X, Y)
	(p.b == 0) && return [zeros(2), zeros(2)]
	x, y = inpanelframe(p, X, Y)
	b = p.b
	xb = x/b
	if abs(y) != 0
		y4b = y/4b
		logstuff = log(norm2(x-b, y) / norm2(x+b, y))
		atanstuff =  atan((x-b)/y) - atan((x+b)/y)
		u1 = (y4b*logstuff - .5atanstuff*(xb - 1.))pi2_inv
		u2 = -(y4b*logstuff - .5atanstuff*(xb + 1.))pi2_inv
		v1 = (1. + 2y4b*atanstuff + .25(xb - 1.)*logstuff)pi2_inv
		v2 = -(1. + 2y4b*atanstuff + .25(xb + 1.)*logstuff)pi2_inv
	else
		logstuff = log((b-x)^2/(b+x)^2)
		v1 = (1. + .25(xb - 1.)logstuff)pi2_inv
		v2 = -(1. + .25(xb + 1.)logstuff)pi2_inv
		u1 = 0
		u2 = 0
	end
	return [[u1, u2], [v1, v2]]
end

"""
function getmotionkelvin(p::LinearPanel, ω, Ut, Vt, X0, Y0)
	Get the circulation induced by the panel motion.

	Input:
		-p  : the linear panel to integrate on
		-ω  : angular velocity in the z direction
		-Ut : translational horizontal body velocity
		-Vt : translational vertical body velocity
		-X0 : position X of the center of rotation of the body
		-Y0 : position Y of the center of rotation of the body
"""
function getmotionkelvin(p::LinearPanel, ω, Ut, Vt, X0, Y0)
	x0, y0 = inpanelframe(p, X0, Y0)
	ut, vt = rotate(Ut, Vt, π+p.θ)
	return 2p.b*(ut + ω*y0)
end

"""
function getmotionparams(p::LinearPanel, X, Y)
	Get the velocity parameters at point X, Y due to the motion of the panel.
	Integrals are performed in the panel frame of reference then re-cast
	in the inertial frame through a rotation.

	Input:
		-p  : the linear panel to integrate on
		-X  : position X where we want the velocity induced by the body motion
		-Y  : position Y where we want the velocity induced by the body motion
		-ω  : angular velocity in the z direction
		-Ut : translational horizontal body velocity
		-Vt : translational vertical body velocity
		-X0 : position X of the center of rotation of the body
		-Y0 : position Y of the center of rotation of the body
		-vs : stretching velocity
"""
function getmotionparams(p::LinearPanel, X, Y, ω, Ut, Vt, X0, Y0, vs=0)
	(p.b == 0) && return zeros(2)
	x, y = inpanelframe(p, X, Y)
	x0, y0 = inpanelframe(p, X0, Y0)
	ut, vt = rotate(Ut, Vt, π+p.θ)
	k = ut + ω*(y + y0) + .5vs*x/p.b
	l = vt + ω*(x - x0) - .5vs*y/p.b
	atanstuff = 0
    logstuff = .5log(norm2(x-p.b, y)/norm2(x+p.b, y))
    y != 0 && (atanstuff = atan((x-p.b)/y) - atan((p.b+x)/y))
	u = (atanstuff*k + l*logstuff + 2ω*p.b)pi2_inv
	v = (l*atanstuff - k*logstuff)pi2_inv
	return rotate(u, v, -π-p.θ)
end

"""
function getinducedvelocity(p::LinearPanel, XY, k, γL, γR, ω, Ut, Vt, Xp, Yp, vs)
    Computes the velocity due to pane `p` at location `XY`
    k specifies which integral is required:
        0 full integral; vortex
        1 full integral; source
        2 half-left integral; vortex
        3 half-left integral; source
        4 half-right integral; vortex
        5 half-right integral; source
"""
function getinducedvelocity(p::LinearPanel, XY, k, γL, γR, ω, Ut, Vt, Xp, Yp, vs=0)
    @assert 0 ≤ k < 6
	x, y = inpanelframe(p, XY...)
	xp, yp = inpanelframe(p, Xp, Yp)
	ut, vt = rotate(Ut, Vt, π+p.θ)
    u = ut + ω*yp
    v = vt - ω*xp
    function p1(α, β, ξ)
        atanstuff = 0
        y != 0 && (atanstuff = atan((x-ξ)/y))
        logstuff = log(norm2(x-ξ, y))
        return [(α*x+β)*atanstuff - .5α*y*logstuff;
                -.5(α*x+β)*logstuff + α*((x-ξ) - y*atanstuff)].*pi2_inv
    end
    iseven(k) && begin
        α = .5(γR - γL + vs)/p.b
        β = .5(γR + γL + 2u)
        k == 0 && return rotate((p1(α, β, p.b) - p1(α, β, -p.b))..., -π-p.θ)
        x0 = -β/α
        k == 2 && return rotate((p1(α, β, x0) - p1(α, β, -p.b))..., -π-p.θ)
        k == 4 && return rotate((p1(α, β, p.b) - p1(α, β, x0))..., -π-p.θ)
    end
    p2(α, β, ξ) = [0 -1;1 0]*p1(α, β, ξ)
    k == 1 && return rotate((p2(ω, v, p.b) - p2(ω, v, -p.b))..., -π-p.θ)
    x0 = -v/ω
    k == 3 && return rotate((p2(ω, v, x0) - p2(ω, v, -p.b))..., -π-p.θ)
    k == 5 && return rotate((p2(ω, v, p.b) - p2(ω, v, x0))..., -π-p.θ)
end

"""
function getimpulse(p::LinearPanel, γ1, γ2)
	Compute the impulse by integrating X^γ and X^X^γ analytically on the linear panel.
	You can specify if you need to add the integral of X^(n^Vb) where n
	points towards the fluid.
	The integral is performed in the panel frame and the result is then
	rotated in the inertial frame of reference.
	Inertial frame: X
	Panel frame:   -x + x'

	Input:
		-p  : the linear panel to integrate on
		-γL : vortex strength on the left side of the panel (in its frame)
		-γR : vortex strength on the right side of the panel
		-ω  : angular velocity in the z direction
		-Ut : translational horizontal body velocity
		-Vt : translational vertical body velocity
		-X0 : position X of origin of the inertial frame
		-Y0 : position Y of origin of the inertial frame
		-Xp : position X of the center of rotation of the body
		-Yp : position Y of the center of rotation of the body
		-nc : need body velocity contribution?
		-vs : stretching velocity
"""
function getimpulse(p::LinearPanel, γL, γR, ω, Ut, Vt, X0, Y0, Xp, Yp, nc, vs=0)
	Θ = π+p.θ
	x, y = inpanelframe(p, X0, Y0)
	_, yp = inpanelframe(p, Xp, Yp)
	ut, _ = rotate(Ut, Vt, Θ)

	b3 = p.b/3
	xy2 = norm2(x, y)
	Σ = p.b*(γR + γL)
	Δ = p.b*(γR - γL)
	m = -.5b3*(p.b*Σ - 2x*Δ) - .5Σ*xy2
	Δ *= b3
	!nc && return [rotate(-y*Σ, x*Σ + Δ, -Θ)..., m]

	vb = (ut + ω*yp)*p.b
	b23 = p.b*b3
	vsb3 = vs*b23
	Σ += 2vb
	return [rotate(-y*Σ, x*Σ - vsb3 + Δ, -Θ)..., m - vb*(xy2 + b23) - x*vsb3]
end

"""
function getmonopolevortex()
    Compute the monopole of vorticity of a linear panel and its location in the reference frame.
    If vortex strength crosses the x axis, then the panel is split into to smaller
    panels of positive and negative strengths.
    In the case it is split, then the id is modified to be able to retrace it later.

    Input:
		-p  : the linear panel to integrate on
        -id : the id of the panel
		-γL : vortex strength on the left side of the panel (in its frame)
		-γR : vortex strength on the right side of the panel
"""
function getmonopolevortex(p::LinearPanel, id, γL, γR, ω, Ut, Vt, Xp, Yp, vs=0)
	Θ = π+p.θ
	yp, yp = inpanelframe(p, Xp, Yp)
	ut, _ = rotate(Ut, Vt, Θ)
    u = ut + ω*yp
    γL += u - .5vs
    γR += u + .5vs

    γL*γR ≥ 0 && return [[getcentroid(p, γL, γR)..., 1im*p.b*(γR + γL), id]]
    #= println("vortex panel $id is split") =#

    x0 = p.b*(γL + γR)/(γL - γR)
    xc1 = (x0^3 - 3x0*p.b^2 - 2p.b^3)/3(p.b+x0)^2
    xc2 = (x0^3 - 3x0*p.b^2 + 2p.b^3)/3(p.b-x0)^2
    m1 = [inreferenceframe(p, xc1, 0)..., .5im*γL*(p.b+x0), -id]
    m2 = [inreferenceframe(p, xc2, 0)..., .5im*γR*(p.b-x0), -id*1im]
    return [m1, m2]
end

"""
function getmonopolesource()
    Compute the monopole of source of a linear panel and its location in the reference frame.
    If vortex strength crosses the x axis, then the panel is split into to smaller
    panels of positive and negative strengths.
    In the case it is split, then the id is modified to be able to retrace it later.

    Input:
		-p  : the linear panel to integrate on
        -id : the id of the panel
"""
function getmonopolesource(p::LinearPanel, id, ω, Ut, Vt, Xp, Yp)
    @assert id < 1e9
	Θ = π+p.θ
	xp, _ = inpanelframe(p, Xp, Yp)
	_, vt = rotate(Ut, Vt, Θ)
    v = vt - ω*xp
    x0 = -v/ω

    abs(x0) ≥ p.b &&
        return [[inreferenceframe(p, ω*p.b^2/3v , 0)..., 2p.b*v, id]]

    #= println("source panel $id is split") =#
    xc1 = (x0^3 - 3x0*p.b^2 - 2p.b^3)/3(p.b+x0)^2
    xc2 = (x0^3 - 3x0*p.b^2 + 2p.b^3)/3(p.b-x0)^2
    m1 = [inreferenceframe(p, xc1, 0)..., .5*(p.b+x0)*(v-ω*p.b), -id]
    m2 = [inreferenceframe(p, xc2, 0)..., .5*(p.b-x0)*(v+ω*p.b), -id*1im]
    return [m1, m2]
end

"""
function getcentroid(p::LinearPanel)
	returns the absolute position of the vorticity centroid of the panel.
"""
function getcentroid(p::LinearPanel, γL, γR)
	xc = p.b*(γR - γL)/3(γR + γL)
	return inreferenceframe(p, xc, 0)
end

"""
function getpressureforce(p::LinearPanel)
"""
function getpressureforce(p::LinearPanel, p1, p2, Xp, Yp)
	x, _ = inpanelframe(p, Xp, Yp)
	[rotate(0, (p1+p2)*p.b, -π-p.θ)..., p.b*((p2-p1)*p.b/3 - (p1+p2)*x)]
end

"""
function getpotential(p::LinearPanel, γL, γR, ω, Ut, Xp, Yp)
	Compute the potential through the integral of tangent velocity
	along the panel "p".
"""
function getpotential(p::LinearPanel, γL, γR, ω, Ut, Vt, Xp, Yp)
	_, yp = inpanelframe(p, Xp, Yp)
	ut, _ = rotate(Ut, Vt, π+p.θ)

	vv = ut + ω*yp
	return p.b*((γL+γR) + 2vv)
end

function getenergyflux(p::LinearPanel, γL, γR, ω, Ut, Vt, X0, Y0, Xp, Yp, vs=0)
	xp, yp = inpanelframe(p, Xp, Yp)
	x, y = inpanelframe(p, X0, Y0)
	ut, vt = rotate(Ut, Vt, π+p.θ)

	αb = .5(γR - γL + vs)
	uu = .5(γL + γR) + ut + ω*yp
	vv = vt - ω*xp
	ωb = ω*p.b

	αωb2 = 2αb*ωb
	Δb = (ωb^2 - αb^2)
	uv2 = 2vv*uu

	I1 = p.b*(αωb2/3 + uv2)
	I2 = p.b*(Δb/3 + (vv + uu)*(vv - uu))
	I3 = p.b*((αωb2*y- Δb*x + 2p.b*(ωb*vv - αb*uu))/3 + y*uv2 - x*(vv + uu)*(vv - uu))
	return [rotate(I1, I2, -π-p.θ)..., I3]
end

function getangularimpulsecorrection(p::LinearPanel, ω, Up, Vp, Ut, Vt, X0, Y0, Xp, Yp)
	xp, _ = inpanelframe(p, Xp, Yp)
	x, y = inpanelframe(p, X0, Y0)
	_, vp = rotate(Up, Vp, π+p.θ)
	ut, vt = rotate(Ut, Vt, π+p.θ)

	vv = vp - ω*xp
	return -2p.b*(p.b*p.b*vt*ω/3 + (vt*x-y*ut)*vv)
end

"""
function getbodymomentofinertia(p::LinearPanel, X0, Y0)
	returns the contribution of panel "p" on the moment of inertia of
	the body.
"""
function getbodymomentofinertia(p::LinearPanel, X0, Y0)
	x, y = inpanelframe(p, X0, Y0)
	return .5*y*p.b*((x^2 + y^2) + p.b^2/3)
end

"""
function getbodyvolume(p::LinearPanel, X0, Y0)
	returns the contribution of panel "p" on the body volume.
"""
function getbodyvolume(p::LinearPanel, X0, Y0)
	_, y = inpanelframe(p, X0, Y0)
	return y*p.b
end

"""
function getbodycentroid(p::LinearPanel, X0, Y0)
	returns the contribution of panel "p" on the position of the
	centroid of the body with respect to position "X0, Y0",
	multiplied by the body volume.
"""
function getbodycentroid(p::LinearPanel, X0, Y0)
	x, y = inpanelframe(p, X0, Y0)
	return rotate(0, -p.b*((x^2 + y^2) + p.b^2/3), -π-p.θ)
end

function getlinearspeedcoeffs(p::LinearPanel, o::LinearPanel, vertical)
	u, v = getparams(o, p.Xc, p.Yc)
	uv1 = rotate(u[1], v[1], p.θ - o.θ)
	uv2 = rotate(u[2], v[2], p.θ - o.θ)
	vertical && return uv1[2], uv2[2]
	return uv1[1], uv2[1]
end

function getconstantspeedcoeffs(p::LinearPanel, o::ConstantPanel, vertical)
	u, v = getparams(o, p.Xc, p.Yc)
	uv = rotate(u, v, p.θ - o.θ)
	vertical && return uv[2]
	return uv[1]
end

"""
function getnormal(p::LinearPanel)
	Get the normal of the panel in its frame of reference
	(inside of the body)
"""
function getnormal(p::LinearPanel)
	return [-sin(p.θ), -cos(p.θ)]
end

function getγ(p::LinearPanel, γ1, γ2, X, Y)
	x, y = inpanelframe(p, X, Y)
	@assert	abs(y) < 1e-12
	return (γ2 - γ1)*x/2p.b + .5(γ2+γ1)
end

"""
function get∫fdV(pa::LinearPanel, XY, f)
	Computes the integral of 'f' on the triangle created by the
	panel and point XY using a Gauss-Legendre quadrature.
	f can return a vector.
"""
function get∫fdV(pa::LinearPanel, XY, f)
	J =	(XY[1] - pa.X2)*(pa.Y1 - pa.Y2) - (pa.X1 - pa.X2)*(XY[2] - pa.Y2)
	x(u, v)	= pa.X2 + u*(XY[1] - pa.X2) + v*(pa.X1 - pa.X2)
	y(u, v)	= pa.Y2 + u*(XY[2] - pa.Y2) + v*(pa.Y1 - pa.Y2)
	Is = gausslegendretriangle(f, x, y)
	return abs(J).*Is
end
