module VSFlow

using LinearAlgebra
using Printf
using Logging
using Random
using Plots

include("VortexElement.jl")

const EPS = eps()
Base.show(io::IO, x::Union{Float64,Float32}) = Base.Grisu._show(io, x, Base.Grisu.SHORTEST, 0, true, false)
"""
PROFILE object
"""
mutable struct Profile
	profile::Function
	N::Int64
	panels::Array{LinearPanel, 1}
	V::Float64
	vortex_points::Array{VortexPoint, 1}
	constant_panel::ConstantPanel
	old_α::Float64
	ω::Float64
	γs::Vector{Float64}
	γg::Float64
	average_length::Float64
	dt::Float64
	T::Float64
	δ::Float64
	ϵ::Float64
	timelapse::Float64
	is_init::Bool
	θg::Float64
	θ0::Float64
	θ1::Float64
	θ2::Float64
	old_impulse::Array{Array{Float64, 1}}
	old_nocaimpulse::Array{Array{Float64, 1}}
	old_ϕs::Array{Array{Float64, 1}}
	old_panel::ConstantPanel
	pivot_location::Array{Float64, 1}
	Ut::Float64
	Vt::Float64
	force_control::Array{Float64, 1}
	trajectoryx::Array{Float64, 1}
	trajectoryy::Array{Float64, 1}
	trailing::Array{Float64, 1}
	fname::String
	Ainv::Array{Float64, 2}
	cps::Array{Float64, 1}
	pcps::Array{Float64, 1}
	Bs::Vector{Float64}
	#Lumping
	η::Float64
	τ::Int64
	Tmin::Int64
	v_corr::Vector{Float64}
	sheet_size::Int64
	last_vp_ind::Int64
	rng::MersenneTwister
	profID::String
end

"""
Constructor
Profile(profile, N, position, args...)
"""
function Profile(id, profile::Function, N, position,
				 dt, T, δ, ϵ, rng, is_comb, lump, args...)
	iseven(N) || error("Need an even number of panels!")
	profID = "profile-$id"
	η = 0
	Tmin = 0
	sheet_size = 0
	if lump
		η = args[1]
		Tmin = args[2]
		sheet_size = args[3]
	end
	position[1] -= .25
	panels = genprofile(profile, N, profID, position[1:end-1], args[end])
	γs = zeros(N + 1)
	average_length = sum(2 .* map(p->p.b, panels))/N

	vortex_points = dipolecombgenerator(rng, T, δ, is_comb)
	trajectoryx = [position[1] + .25]
	trajectoryy = [position[2]]
	fname = "naca00"*string(args[end])*"_np"*string(N)*"_dt"*string(dt)*"_T"*string(T)*"_dv"*string(δ)*"_eps"*string(ϵ)
	lump && (fname = fname*"_lump_errMax"*string(η)*"_Tmin"*string(Tmin)*"_ss"*string(sheet_size))
	fname = replace(fname, "."=>"")

	pivot_location = [position[1] + .25, position[2]]
	V = sum(map(pa->getbodyvolume(pa, pivot_location...), panels))
	for p in panels
		setangle!(p, position[end], pivot_location)
	end

	θ0 = panels[1].θ - panels[end].θ + sign(panels[end].θ)*pi
	θ0 = abs(atan(sin(θ0), cos(θ0)))
	θ1 = .5θ0
	θ2 = .5θ0
	θg = panels[1].θ - θ1

	X1 = .5(panels[1].X1+panels[end].X2)
	Y1 = .5(panels[1].Y1+panels[end].Y2)
	trailing = [X1, Y1]

	#Constant panel positioning
	b = dt - ϵ*average_length
	X2 = X1 + 2b*cos(θg)
	Y2 = Y1 - 2b*sin(θg)
	constant_panel = BuildPanel(ConstantPanel, X1, Y1, X2, Y2, profID)
	cropedge!(constant_panel, ϵ*average_length, false)
	cropedge!(panels[1], ϵ*average_length, false)
	cropedge!(panels[end], ϵ*average_length, true)
	p0 = [zeros(3), zeros(3), zeros(3), zeros(3)]
	n0 = [zeros(3), zeros(3), zeros(3), zeros(3)]
	ϕs0 = [zeros(N+1), zeros(N+1), zeros(N+1), zeros(N+1)]
	return Profile(profile, N, panels, V, vortex_points, constant_panel, position[end],
				   0, γs, 0, average_length, dt, T, δ, ϵ, 0,
				   true, θg, θ0, θ1, θ2, p0, n0, ϕs0, constant_panel, pivot_location, 0, 0,
				   zeros(3), trajectoryx, trajectoryy, trailing, fname, zeros(N+2, N+2),
				   zeros(N+1), zeros(N), zeros(N), η, 0, Tmin, zeros(2), sheet_size, 1, rng, profID)
end

function copier!(p::Profile, p2::Profile)
	p.profile = p2.profile
	p.N = p2.N
	p.panels = p2.panels
	p.V = p2.V
	p.vortex_points = p2.vortex_points
	p.constant_panel = p2.constant_panel
	p.old_α = p2.old_α
	p.ω = p2.ω
	p.γs = p2.γs
	p.γg = p2.γg
	p.average_length = p2.average_length
	p.dt = p2.dt
	p.T = p2.T
	p.δ = p2.δ
	p.ϵ = p2.ϵ
	p.timelapse = p2.timelapse
	p.is_init = p2.is_init
	p.θg = p2.θg
	p.θ0 = p2.θ0
	p.θ1 = p2.θ1
	p.θ2 = p2.θ2
	p.old_impulse = p2.old_impulse
	p.old_nocaimpulse = p2.old_nocaimpulse
	p.old_ϕs = p2.old_ϕs
	p.old_panel = p2.old_panel
	p.pivot_location = p2.pivot_location
	p.Ut = p2.Ut
	p.Vt = p2.Vt
	p.force_control = p2.force_control
	p.trajectoryx = p2.trajectoryx
	p.trajectoryy = p2.trajectoryy
	p.trailing = p2.trailing
	p.fname = p2.fname
	p.Ainv = p2.Ainv
	p.cps = p2.cps
	p.pcps = p2.pcps
	p.η = p2.η
	p.τ = p2.τ
	p.Tmin = p2.Tmin
	p.v_corr = p2.v_corr
	p.sheet_size = p2.sheet_size
	p.last_vp_ind = p2.last_vp_ind
	p.rng = p2.rng
	p.profID = p2.profID
	return nothing
end

"""
    genprofile(profile, N, position, args...)

Generate a profile with a wedged trailing edge.
"""
function genprofile(profile::Function, N::Int64, profID, position::Array, args...)
	θ = range(0, pi, length=N÷2 + 1)
	x = @. .5cos(θ) + .5

	if length(args) == 1
		ye, yi = profile(x, args[1])
	else
		ye, yi = profile(x)
	end

	L = []
	for i in 1:length(x)-1
		@inbounds	push!(L, BuildPanel(LinearPanel, x[i] + position[1],
										  ye[i] + position[2], x[i+1] + position[1], ye[i+1] + position[2], profID))
	end
	for i in length(x)-1:-1:1
		@inbounds push!(L, BuildPanel(LinearPanel, x[i+1] + position[1],
										 yi[i+1] + position[2], x[i] + position[1], yi[i] + position[2], profID))
	end
	return L
end

"""
getprofilelengths(p::Profile)
"""
function getprofilelengths(p::Profile)
	λ = 2map(pr->pr.b, p.panels)
	pushfirst!(λ, 0)
	return cumsum(λ)
end


"""
inbodyframe(p::Profile, X, Y)
"""
function inbodyframe(p::Profile, X, Y)
	X0, Y0 = p.pivot_location
	x, y = rotate(X-X0, Y-Y0, p.old_α)
	return [x, y]
end

"""
getparamsinrefframe(panel, X, Y)
"""
function getparamsinrefframe(panel::LinearPanel, X, Y)
	u, v = getparams(panel, X, Y)
	uv1 = rotate(u[1], v[1], -(π+panel.θ))
	uv2 = rotate(u[2], v[2], -(π+panel.θ))
	U = [uv1[1], uv2[1]]
	V = [uv1[2], uv2[2]]
	return [U, V]
end

"""
getparamsinrefframe(panel, X, Y)
"""
function getparamsinrefframe(panel::ConstantPanel, X, Y)
	u, v = getparams(panel, X, Y)
	return rotate(u, v, -(π+panel.θ))
end

"""
getboundpanelsinducedvel(profile, X, Y)
"""
function getboundpanelsinducedvel(p::Profile, X, Y)
	X0, Y0 = p.pivot_location
	U = V = 0
	for (i, panel) in enumerate(p.panels)
		u, v = getparamsinrefframe(panel, X, Y)
		@inbounds U += u[1]*p.γs[i] + u[2]*p.γs[i+1]
		@inbounds V += v[1]*p.γs[i] + v[2]*p.γs[i+1]
	end
	int = sum(map((pa)->getmotionparams(pa, X, Y, p.ω, p.Ut, p.Vt, X0, Y0), p.panels))

	return [U, V] + int
end

"""
motionupdate!(p::Profile, X, Ẋ, t)
"""
function motionupdate!(p::Profile, X, Ẋ, t)
	X0, Y0 = p.pivot_location
	Xt, Yt, θ = X
	α = -θ
	p.Ut, p.Vt, p.ω = Ẋ
	for panel in p.panels
		translate!(panel, Xt-X0, Yt-Y0)
		setangle!(panel, α-p.old_α, [Xt, Yt])
	end
	p.old_α = α
	p.pivot_location = [Xt, Yt]
	X1 = .5(p.panels[1].X1 + p.panels[end].X2)
	Y1 = .5(p.panels[1].Y1 + p.panels[end].Y2)
	p.trailing = [X1, Y1]
	return nothing
end

"""
dipolecombgenerator(rng::MersenneTwister, T, δ, threshold=1e-2)
"""
function dipolecombgenerator(rng::MersenneTwister, T, δ, is_comb, threshold=5e-3)
	dipoles = []
	!is_comb && return dipoles
	n_dip = Random.rand(rng, 1:floor(UInt64, .4*T))

	d = .4 + .8Random.rand(rng)
	v_center = .05 + .25Random.rand(rng)
	Γ = .5pi*d*v_center
	lim = sqrt(.5d*Γ/(pi*threshold) + (.5d)^2)
	sn1 = Random.rand(rng, [-1., 1.])
	sn2 = Random.rand(rng, [-1., 1.])
	X = -(lim + (.5T - lim)Random.rand(rng))
	Y = .3 + .8Random.rand(rng) - .5Γ*X/(pi*d)
	push!(dipoles, VortexPoint(X + .5d, sn1*Y, -sn2*Γ))
	push!(dipoles, VortexPoint(X - .5d, sn1*Y, sn2*Γ))
	last_d = d

	for _ in 1:n_dip
		d = .4 + .8Random.rand(rng)
		v_center = .05 + .25Random.rand(rng)
		Γ = .5pi*d*v_center
		vp = dipoles[end]
		old_lim = sqrt(.5last_d*abs(vp.Γ)/(pi*threshold) + (.5d)^2)
		lim = sqrt(.5d*Γ/(pi*threshold) + (.5d)^2)
		sn1 = Random.rand(rng, [-1., 1.])
		sn2 = Random.rand(rng, [-1., 1.])
		last_X = -vp.X + .5last_d + max(lim, old_lim)
		X = -(last_X + (.5T - last_X)Random.rand(rng))
		Y = .3 + .8Random.rand(rng) - .5Γ*X/(pi*d)
		if T - last_X >= 0
			push!(dipoles, VortexPoint(X + .5d, sn1*Y, -sn2*Γ))
			push!(dipoles, VortexPoint(X - .5d, sn1*Y, sn2*Γ))
			last_d = d
		end
	end
	return dipoles
end

"""
globalvelocity(p::Profile, XY)
"""
function globalvelocity(p::Profile, XY, incloud=false)
	UV = getboundpanelsinducedvel(p, XY[1], XY[2])
	UV += p.γg*getparamsinrefframe(p.constant_panel, XY[1], XY[2])
	!incloud && !isempty(p.vortex_points) && (UV += sum(map(vp->getinducedvelocity(vp, XY[1], XY[2], p.δ), p.vortex_points)))
	return UV
end

"""
setforcecontrol!(p::Profile, fx, fy, mz)
"""
function setforcecontrol!(p::Profile, fx, fy, mz)
	p.force_control = [fx, fy, -mz]
	return nothing
end

"""
setinitvelocity(p::Profile, u, v, α̇)
"""
function setinitvelocity(p::Profile, u, v, α̇)
	p.Ut = u
	p.Vt = v
	p.ω = -α̇
	return nothing
end

"""
gettrailingedgevel!(p::Profile)
"""
function gettrailingedgevel!(p::Profile)
	v1 = p.γs[1]
	v2 = -p.γs[end]
    β = atan(tan(.5p.θ0)*(v2-v1)/(v2+v1))
    p.θ1 = .5p.θ0 + β
    p.θ2 = .5p.θ0 - β
	p.θg = p.panels[1].θ - p.θ1
    u3 = .5*(v1*cos(p.θ1) + v2*cos(p.θ2))
    p.Bs[5] = β

	return [-u3*cos(p.θg), u3*sin(p.θg)]
end

"""
placeconstantpanel!(p::Profile, X2, Y2)
"""
function placeconstantpanel!(p::Profile, X2, Y2)
	X1 = p.trailing[1] + p.ϵ*p.average_length*cos(p.θg)
	Y1 = p.trailing[2] - p.ϵ*p.average_length*sin(p.θg)
	p.constant_panel = BuildPanel(ConstantPanel, X1, Y1, X2, Y2, p.profID)
	return nothing
end

"""
step!(p::Profile)
"""
function step!(p::Profile, is4thorder=true, need_reset=true)
	X0, Y0 = p.pivot_location

	X1 = p.trailing[1] + p.ϵ*p.average_length*cos(p.θg)
	Y1 = p.trailing[2] - p.ϵ*p.average_length*sin(p.θg)

	XY0 = map(vp->[vp.X, vp.Y], p.vortex_points)
	XA0 = [X0, Y0, -p.old_α]
	ẊA0 = [p.Ut, p.Vt, p.ω]
	XP0 = [X1, Y1]

	order = 1 + is4thorder*3
	β = [0, .5, .5, 1]
	ζ = is4thorder*[.5, 1., 1., .5]/3. + !is4thorder*[1., 0, 0, 0]
	dXY = XY0
	dXA = XA0
	dXP = XP0
	dẊA = ẊA0
	XY = XY0
	XA = XA0
	XP = XP0
	ẊA = ẊA0
	XYL = XY0
	XAL = XA0
	XPL = XP0
	ẊAL = ẊA0
	for k in 1:order
		if k > 1
			@inbounds XYL = XY0 + β[k]*dXY
			@inbounds XAL = XA0 + β[k]*dXA
			@inbounds XPL = XP0 + β[k]*dXP
			@inbounds ẊAL = ẊA0 + β[k]*dẊA
			for (i, vp) in enumerate(p.vortex_points)
				@inbounds vp.X = XYL[i][1]
				@inbounds vp.Y = XYL[i][2]
			end
			@inbounds motionupdate!(p, XAL, ẊAL, β[k]*p.dt)
			placeconstantpanel!(p, XPL...)
			setγs!(p)
		end
		#= dXY = p.dt*map(x->globalvelocity(p, x), XYL) =#
		dXY = p.dt*(getcloudvelocities(p.vortex_points, p.δ) .+ map(x->globalvelocity(p, x, true), XYL))

		!isempty(p.vortex_points) && (dXY[p.last_vp_ind] += p.vortex_points[p.last_vp_ind].dv*p.dt)
		dXP = p.dt*gettrailingedgevel!(p)
		dẊA = p.dt*p.force_control
		dXA = p.dt*ẊAL
		@inbounds XY += dXY*ζ[k]
		@inbounds XP += dXP*ζ[k]
		@inbounds ẊA += dẊA*ζ[k]
		@inbounds XA += dXA*ζ[k]
	end

	for (i, vp) in enumerate(p.vortex_points)
		@inbounds vp.X = XY[i][1]
		@inbounds vp.Y = XY[i][2]
	end
	motionupdate!(p, XA, ẊA, p.dt)
	placeconstantpanel!(p, XP...)
	setγs!(p)
	gettrailingedgevel!(p)

	push!(p.trajectoryx, X0)
	push!(p.trajectoryy, Y0)

	need_reset && resetconstantpanel!(p)

	p.timelapse += p.dt
	p.τ += 1
	return nothing
end

"""
function resetconstantpanel!(p)
	Transform the constant panel into a point vortex.
	Then recompute the vortex strength distribution on the body.
"""
function resetconstantpanel!(p)
	getnoca!(p, true)
	X = .5(p.constant_panel.X1+p.constant_panel.X2)
	Y = .5(p.constant_panel.Y1+p.constant_panel.Y2)
	Γ = 2p.γg*p.constant_panel.b

	if Γ != 0
		new_vp = VortexPoint(X, Y, Γ)
		push!(p.vortex_points, new_vp)
	end

	X1 = p.trailing[1] + p.ϵ*p.average_length*cos(p.θg)
	Y1 = p.trailing[2] - p.ϵ*p.average_length*sin(p.θg)
	p.old_panel = p.constant_panel
	placeconstantpanel!(p, X1, Y1)
	p.Bs[1] = p.γg

	p.γg = 0
	setγs!(p)
	return nothing
end

"""
lumpvortices!(p::Profile)
"""
function lumpvortices!(p::Profile)
	if p.last_vp_ind > length(p.vortex_points) - p.sheet_size + 1
		step!(p)
		p.τ = 0
		return 0
	end

	t_index = p.last_vp_ind
	s_index = p.last_vp_ind + 1
	if sign(p.vortex_points[t_index].Γ) != sign(p.vortex_points[s_index].Γ)
		p.τ = 0
		p.vortex_points[t_index].dv = zeros(2)
		p.last_vp_ind += 1
		step!(p)
		return 0
	end

	γ̇ = p.vortex_points[s_index].Γ/p.dt
	s_vp = p.vortex_points[s_index]
	t_vp = p.vortex_points[t_index]
	s_imp = getvorteximpulse(p, s_vp)
	t_imp = getvorteximpulse(p, t_vp)

	grads = zeros(p.N+1, 2)
	for (i, panel) in enumerate(p.panels)
		@inbounds grads[i, :] = getvelocitygradient(t_vp, panel.Xc, panel.Yc, π+panel.θ, p.δ)
	end
	∇γ = (p.Ainv*grads)'

	κ = map(pa->[pa.Y1, -pa.X1], p.panels)
	push!(κ, [p.panels[end].Y2, -p.panels[end].X2])
	κ = hcat(κ...)
	A = mapslices(fook->[mapslices(foog->trapz(fook.*foog, getprofilelengths(p)), ∇γ, dims=2)], κ, dims=2)

	∇p = [0 1; -1 0] + hcat(A...)'
	p.v_corr = (∇p \ (s_imp - t_imp))γ̇/t_vp.Γ
	η = getforceerror!(p, s_index, t_index)
	return η
end

"""
getforceerror!(p::Profile, s_index, t_index)
"""
function getforceerror!(p::Profile, s_index, t_index)
	no_transfer_p = deepcopy(p)
	no_transfer_p.τ = 0
	no_transfer_p.vortex_points[t_index].dv = zeros(2)
	no_transfer_p.last_vp_ind += 1
	step!(no_transfer_p)
	no_transfer_imp = getimpulse(no_transfer_p)[1:2]

	transfer!(p, s_index, t_index)
	step!(p)
	imp = getimpulse(p)[1:2]

	η = norm(imp - no_transfer_imp)/p.dt
	if η > p.η && p.τ >= p.Tmin
		@info "No Lumping allowed"
		copier!(p, no_transfer_p)
		η = 0.
	end
	no_transfer_p = nothing
	return η
end

"""
getvorteximpulse(p, vp)
"""
function getvorteximpulse(p::Profile, vp)
	b = zeros(p.N + 1)
	for (m, panel) in enumerate(p.panels)
		U, V = getinducedvelocity(vp, panel.Xc, panel.Yc, p.δ)
		_, V = rotate(U, V, π+panel.θ)
		@inbounds b[m] = -V
	end
	b[end] = -vp.Γ
	γs = p.Ainv * b

	X0, Y0 = p.trajectoryx[1], p.trajectoryy[1]
	Xp, Yp = p.pivot_location
	panels_impulse = sum(map((pa, γL, γR)->getimpulse(pa, γL, γR, p.ω, p.Ut, p.Vt, X0, Y0, Xp, Yp, false)[1:2], p.panels, γs[1:end-1], γs[2:end]))
	return [vp.Y-Y0, X0-vp.X] + panels_impulse/vp.Γ
end

"""
transfer!(p, s_index, t_index)
"""
function transfer!(p::Profile, s_index, t_index)
	p.vortex_points[t_index].Γ += p.vortex_points[s_index].Γ
	p.vortex_points[t_index].dv = p.v_corr
	deleteat!(p.vortex_points, s_index)
	return nothing
end


"""
setγs!(p::Profile)
"""
function setγs!(p::Profile)
	X0, Y0 = p.pivot_location
	A = zeros(p.N + 2, p.N + 2)
	b = zeros(p.N + 2)

	m::UInt64 = 1
	for panel in p.panels

		V = - .5rotate(p.Ut - p.ω*(panel.Yc-Y0), p.Vt + p.ω*(panel.Xc-X0), π+panel.θ)[2]
		for vp in p.vortex_points
			Up, Vp = getinducedvelocity(vp, panel.Xc, panel.Yc, p.δ)
			Up, Vp = rotate(Up, Vp, π+panel.θ)
			V += Vp
		end
		Ub, Vb = sum(map((pa)->getmotionparams(pa, panel.Xc, panel.Yc, p.ω, p.Ut, p.Vt, X0, Y0), p.panels))

		@inbounds b[m] -= V + rotate(Ub, Vb, π+panel.θ)[2]

		n::UInt64 = 1
		for other in p.panels
			v1, v2 = getlinearspeedcoeffs(panel, other, true)
			@inbounds A[m, n] += v1
			@inbounds A[m, n+1] += v2
			n += 1
		end
		v = getconstantspeedcoeffs(panel, p.constant_panel, true)
		@inbounds A[m, n+1] = v

		#Kelvin
		@inbounds A[p.N+1, m] += panel.b
		@inbounds A[p.N+1, m+1] += panel.b
		m += 1
	end
	@inbounds A[p.N+1, m+1] += 2p.constant_panel.b

	int = -sum(map(pa->getmotionkelvin(pa, p.ω, p.Ut, p.Vt, X0, Y0), p.panels))
	@inbounds b[p.N+1] = -sum(map(vp->vp.Γ, p.vortex_points)) + int

	#Unsteady Kutta
	@inbounds A[p.N+2, 1] = cos(p.θ1)
	@inbounds A[p.N+2, p.N+1] = cos(p.θ2)
	@inbounds A[p.N+2, p.N+2] = -1

	#= #Giesing-Maskell =#
	#= bo = (p.Bs[1] < 0) =#
	#= A[p.N+2, 1] = 1 =#
	#= A[p.N+3, p.N+1] = 1 =#
	#= A[p.N+2, p.N+2] = -1*bo =#
	#= A[p.N+3, p.N+2] = -1*(!bo) =#

	if p.constant_panel.b == 0
		@inbounds p.Ainv = inv(A[1:p.N+1, 1:p.N+1])
		@inbounds p.γs = p.Ainv * b[1:p.N+1]
	else
		p.Ainv = inv(A)
		#= p.Ainv = pinv(A) =#
		p.γs = p.Ainv * b
	end

	if length(p.γs) == p.N + 2
		p.γg = p.γs[end]
		pop!(p.γs)
	end

	p.is_init = false
	return nothing
end

"""
getimpulse(p::Profile)
"""
function getimpulse(p::Profile)
	X0, Y0 = p.trajectoryx[1], p.trajectoryy[1]
	Xp, Yp = p.pivot_location
	panels_impulse = sum(map((pa, γL, γR)->getimpulse(pa, γL, γR, p.ω, p.Ut, p.Vt, X0, Y0, Xp, Yp, true), p.panels, p.γs[1:end-1], p.γs[2:end]))

	if !isempty(p.vortex_points)
		points_impulse = sum(map(vp->vp.Γ*[vp.Y-Y0, -vp.X+X0, -.5norm2(vp.Y - Y0, vp.X - X0)], p.vortex_points))
	else
		points_impulse = zeros(3)
	end
	return panels_impulse + points_impulse
end


"""
getcoeffimpulse!()
"""
function getcoeffimpulse!(p::Profile, update=true, β=1.)
	X0, Y0 = p.trajectoryx[1], p.trajectoryy[1]
	xyz_imp = getimpulse(p::Profile)

	dt = β*p.dt
	f = -backwarddifference(xyz_imp, p.old_impulse, dt)

	if update
		pop!(p.old_impulse)
		pushfirst!(p.old_impulse, xyz_imp)
	end

	cd, cl, cm = 2f
	cm -= dot(p.pivot_location - [X0, Y0], [cl, -cd])
	return [cd, cl, cm]
end

"""
"""
function getnoca!(p::Profile, update=true)
	X0, Y0 = p.trajectoryx[1], p.trajectoryy[1]
	Xp, Yp = p.pivot_location

	u = 2p.constant_panel.b/p.dt
	panels_impulse = sum(map((pa, γL, γR)->getimpulse(pa, γL, γR, p.ω, p.Ut, p.Vt, Xp, Yp, Xp, Yp, true), p.panels, p.γs[1:end-1], p.γs[2:end]))
	energyflux = sum(map((pa, γL, γR)->getenergyflux(pa, γL, γR, p.ω, 0, 0, Xp, Yp, Xp, Yp), p.panels, p.γs[1:end-1], p.γs[2:end]))

	vorticityflux = p.γg*u*([p.trailing[2]-Yp, Xp - p.trailing[1], -.5norm2.((p.trailing - p.pivot_location)...)])
	Δm = sum(map(pa->getangularimpulsecorrection(pa, p.ω, 0, 0, p.Ut, p.Vt, Xp, Yp, Xp, Yp), p.panels))

	f = energyflux - backwarddifference(panels_impulse, p.old_nocaimpulse, p.dt) - vorticityflux

	if update
		pop!(p.old_nocaimpulse)
		pushfirst!(p.old_nocaimpulse, panels_impulse)
	end

	cd, cl, cm = 2f
	p.Bs[2] = cd
	p.Bs[3] = cl
	p.Bs[4] = cm + Δm
end

"""
getboundpotential(p::Profile, nn=20)
	Integrates the velocity on the profile to obtain the potential.
	Integrating on a path from upstream is required to get the exact
	pressure. If the goal is to get the loads, then the upstream
	integration is unnecessary.
"""
function getboundpotential(p::Profile, nn=20)
	X0, Y0 = p.pivot_location
	ϕs = zeros(p.N+1)
	tmp = 0
	for i in 1:p.N
		@inbounds tmp += getpotential(p.panels[i], p.γs[i], p.γs[i+1], p.ω, p.Ut, p.Vt, X0, Y0)
		@inbounds ϕs[i+1] = tmp
	end
	return ϕs
end

"""
setcps
"""
function setcps!(p::Profile, update=true, β=1.)
	ϕs = getboundpotential(p)
	dt = β*p.dt

	X0, Y0 = p.pivot_location
	Us = map(pa->p.Ut - p.ω*(pa.Y1 - Y0), p.panels)
	Vs = map(pa->p.Vt + p.ω*(pa.X1 - X0), p.panels)
	push!(Us, p.Ut - p.ω*(p.panels[end].Y2 - Y0))
	push!(Vs, p.Vt + p.ω*(p.panels[end].X2 - X0))

	p.cps = 1. .+ norm2.(Us, Vs) .- p.γs.^2 .- 2backwarddifference(ϕs, p.old_ϕs, dt)
	if update
		pop!(p.old_ϕs)
		pushfirst!(p.old_ϕs, ϕs)
	end
	return nothing
end

"""
getcoeffpotential!
"""
function getcoeffpotential!(p::Profile)
	setcps!(p)
	cd, cl, cm = sum(map((pa, p1, p2)->getpressureforce(pa, p1, p2, p.pivot_location...), p.panels, p.cps[1:end-1], p.cps[2:end]))
	return cd, cl, cm
end

"""
getstate(p, n)
	Used for the gym model
"""
function getstate(p, n)
	skip = div(p.N, n)
	m = div(mod(p.N, n), 2)
	setcps!(p)
	pressures = p.cps[skip+m:skip:end-m]

	#= L = [p.pivot_location[2], p.Vt, p.old_α, p.ω] =#
	L = [p.old_α, p.ω]
	io = open("pressures_"*p.fname*".txt", "a")
	write(io, string(pressures)*"\n")
	close(io)
	return vcat(L, pressures)
end

getaoa(p) = rad2deg(atan(-p.Vt, -p.Ut) + p.old_α)

"""
logresults2(p, is_init, t, mem)
	Results for the RL controller
"""
function logresults2(p::Profile, is_init, t=0., mem=0)
	if is_init
		io = open(p.fname*".txt", "w")
		write(io,"time\tcd_i\tcl_i\tcm_i\tpx\tpy\tpz\tΓ_b\tθg\tγg\tα\tCPU_T\tBytes\tError\n")
		close(io)
		return nothing
	else
		cd_i, cl_i, cm_i = getcoeffimpulse!(p)
		imp = getimpulse(p)

		Γb = 2p.γg*p.constant_panel.b + sum(map(vp->vp.Γ, p.vortex_points))
		α = getaoa(p)

		io = open(p.fname*".txt", "a")
		s = @sprintf "%.2f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%d\n" p.timelapse cd_i cl_i cm_i imp[1] imp[2] imp[3] Γb (p.θg-p.old_α)/p.θ0 p.γg α t mem
		write(io, s)
		close(io)
		return [cd_i, cl_i, cm_i]
	end
end

"""
logresults(p, is_init, t, mem)
"""
function logresults(p::Profile, is_init, t=0., mem=0, η=0.)
	if is_init
		io = open(p.fname*".txt", "w")
		write(io,"time\tcd_i\tcl_i\tcm_i\tcd_p\tcl_p\tcm_p\tcd_b\tcl_b\tcm_b\tΓ_b\tθg\tγg\tα\tCPU_T\tBytes\tError\n")
		close(io)
		return nothing
	else
		cd_b, cl_b, cm_b = p.Bs[2], p.Bs[3], p.Bs[4]
		cd_i, cl_i, cm_i = getcoeffimpulse!(p)
		cd_p, cl_p, cm_p = getcoeffpotential!(p)

		Γb = 2p.γg*p.constant_panel.b + sum(map(vp->vp.Γ, p.vortex_points))
		α = getaoa(p)
        θ = p.Bs[5]/p.θ0#p.panels[1].θ - .5p.θ0

		io = open(p.fname*".txt", "a")
		s = @sprintf "%.2f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%f\n" p.timelapse cd_i cl_i cm_i cd_p cl_p cm_p cd_b cl_b cm_b Γb θ p.Bs[1] α t mem η
		write(io, s)
		close(io)
	end
end

"""
"""
function profilerun(p::Profile, is_write, is_lumped, is4thorder, forcefun, motion_args, show=false)
	run(`echo "==================================================================="`)
	run(`figlet -f larry3d Unsteady Panels`)
	run(`echo "==================================================================="`)
	run(pipeline(`git log`, `head -6`))
	run(`echo "==================================================================="`)
	println(" Params:")
	println("---------")
	println("Profile          :"*string(p.profile))
	println("Panel Number     :"*string(p.N))
	println("Timestep         :"*string(p.dt))
	println("Horizon Time     :"*string(p.T))
	println("Blob δ           :"*string(p.δ))
	println("-----------------------------------------------------------")

	is_write && logresults(p, true)
	setγs!(p)

	count::UInt64 = 0
	if show
		anim = Animation()
	end
	while(p.timelapse <= p.T + p.dt)
		if count%10 == 0
			@info @sprintf "%.3f/%.3f" p.timelapse p.T
		end

		η = 0.
		setforcecontrol!(p, forcefun(p.timelapse, motion_args...)[3]...)
		if is_lumped
			η, t, bytes, _, _ = @timed lumpvortices!(p)
		else
			_, t, bytes, _, _ = @timed step!(p, is4thorder)
		end
		is_write && logresults(p, false, t, bytes, η)
		count += 1

		skip=1
		##Animation
		if show && count%skip ==0
			x = map(vp->vp.X, p.vortex_points)
			y = map(vp->vp.Y, p.vortex_points)
			Γ = map(vp->vp.Γ, p.vortex_points)
			traj = plot(p.trajectoryx, p.trajectoryy, line=(1, :lightgray), aspect_ratio=1, xlim=(-.5-p.T, 1.2), ylim=(-3, 3), framestyle=:none, grid=false, ticks=false, colorbar=false, dpi=400)
			pp = scatter!(x, y, zcolor=Γ, markercolor=:coolwarm, markersize=3, markerstrokewidth=0, clims=(-.2,.2))
			xp = map(panel->panel.X1, p.panels)
			yp = map(panel->panel.Y1, p.panels)
			push!(xp, p.panels[end].X2)
			push!(yp, p.panels[end].Y2)
			plot!(pp, xp, yp, line=(2, :black), fill=(0, :black), legend=false)

			frame(anim)
		end
	end

	if show
		gif(anim, p.fname*".gif", fps=20)
	end

	println("===================================================================")
	println("DONE")
	println("===================================================================")
end

end
