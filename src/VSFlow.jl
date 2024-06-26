#!/usr/bin/env julia
# File              : src/VSFlow.jl
# Author            : Denis Dumoulin <denis.dumoulin@uclouvain.be>
# Date              : 30.08.2021
# Last Modified Date: 28.09.2021
module VSFlow

using LinearAlgebra
using Printf
using DelimitedFiles
using ForwardDiff
using Logging
using Plots
gr()

include("VortexElement.jl")
include("Airfoil.jl")
include("History.jl")

export uniform, heavepitch, circularmotion, heavepitch2
export naca00, gaw1, circle, ellipse
export Profile, profilerun
export History, getimpulse, getboundpotential, plotaero, plotcps, plotϕs

"""
PROFILE object
"""
mutable struct Profile
  profileshape::Function
  N::Int64
  panels::Array{LinearPanel,1}
  V::Float64
  vortex_points::Array{VortexPoint,1}
  constant_panel::ConstantPanel
  old_α::Float64
  ω::Float64
  γs::Vector{Float64}
  γg::Float64
  oldγg::Float64
  average_length::Float64
  dt::Float64
  T::Float64
  δ::Float64
  ϵ::Float64
  timelapse::Float64
  is_init::Bool
  θg::Float64
  θ0::Float64
  β::Float64
  θ1::Float64
  θ2::Float64
  old_panel::ConstantPanel
  pivot_location::Array{Float64,1}
  Ut::Float64
  Vt::Float64
  force_control::Array{Float64,1}
  trailing::Array{Float64,1}
  fname::String
  Ainv::Array{Float64,2}

  #Lumping
  η::Float64
  τ::Int64
  Tmin::Int64
  v_corr::Vector{Float64}
  sheet_size::Int64
  last_vp_ind::Int64
  islumped::Bool

  profID::String
  history::History
  isgaussian::Bool
end

"""
    Profile(; id, profileshape::Function, x0, ẋ0, N, dt, T, δ = 1e-2, ϵ = 1e-2, (η, Tmin, Smin) = zeros(3))

Constructs a `Profile` object.
If `eta`, `Tmin`, `Smin` are all `0`, then no lumping operation occurs.

# Keyword Arguments
 - `id`: `String` used to identify the profile.
 - `profileshape`: `function` determining the profile shape.
 - `N`: number of panels on the surface of the profile.
 - `x0`: initial position `[X0, Y0, α0]` of the quarter chord of the profile.
 - `ẋ0`: initial velocity `[Ẋ0, Ẏ0, α̇0]` of the quarter chord of the profile.
 - `dt`: timestep of the simulation.
 - `T`: horizon time.
 - `isgaussian`: `` determining if the vortex kernel is Gaussian or Low-order Algebraic.
 - `δ = 1e-2`: kernel cut-off width.
 - `ϵ = 1e-2`: percentage of the average panel length to chop the trailing edge.
 - `η = 0`: maximum error due to lumping.
 - `Tmin = 0`: minimum time between two successive active vortices.
 - `Smin = 0`: minimum vortex sheet length in the wake.

"""
function Profile(; id, profileshape::Function, x0, ẋ0, N, dt, T,
  isgaussian=false,
  δ=1e-2,
  ϵ=1e-2,
  lumpargs=zeros(3))
  iseven(N) || error("Need an even number of panels!")
  islumped = (lumpargs != zeros(3))
  η, Tmin, sheet_size = lumpargs

  # b = 1 / 3.0
  b = 0.25
  x0[1] -= b
  panels = genprofile(profileshape, N, id, x0[1:end-1])
  γs = zeros(N + 1)
  average_length = sum(2 .* map(p -> p.b, panels)) / N
  # average_length = panels[1].b+panels[end].b

  vortex_points = []
  fname = id * "_np" * string(N) * "_dt" * string(dt) * "_T" * string(T) * "_dv" * string(δ) * "_eps" * string(ϵ)
  islumped && (fname = fname * "_lump_errMax" * string(η) * "_Tmin" * string(Tmin) * "_ss" * string(sheet_size))
  fname = replace(fname, "." => "")

  pivot_location = [x0[1] + b, x0[2]]
  V = sum(map(pa -> getbodyvolume(pa, pivot_location...), panels))
  for p in panels
    setangle!(p, x0[end], pivot_location)
  end
  U0, V0, α̇0 = ẋ0

  θ0 = panels[1].θ - panels[end].θ + sign(panels[end].θ) * pi
  θ0 = abs(atan(sin(θ0), cos(θ0)))
  θ1 = 0.5θ0
  θ2 = 0.5θ0
  θg = panels[1].θ - θ1

  X1 = 0.5(panels[1].X1 + panels[end].X2)
  Y1 = 0.5(panels[1].Y1 + panels[end].Y2)
  trailing = [X1, Y1]

  #Constant panel positioning
  b = dt - ϵ * average_length
  X2 = X1 + 2b * cos(θg)
  Y2 = Y1 - 2b * sin(θg)
  constant_panel = BuildPanel(ConstantPanel, X1, Y1, X2, Y2, id)
  cropedge!(constant_panel, ϵ * average_length, false)
  cropedge!(panels[1], ϵ * average_length, false)
  cropedge!(panels[end], ϵ * average_length, true)
  h = History(N, dt, T)
  return Profile(profileshape, N, panels, V, vortex_points, constant_panel, x0[end],
    -α̇0, γs, 0, 0, average_length, dt, T, δ, ϵ, 0,
    true, θg, θ0, 0, θ1, θ2, constant_panel, pivot_location, U0, V0,
    zeros(3), trailing, fname, zeros(N + 2, N + 2),
    η, 0, Tmin, zeros(2), sheet_size, 1, islumped, id, h, isgaussian)
end

function copier!(p::Profile, p2::Profile)
  p.profileshape = p2.profileshape
  p.N = p2.N
  p.panels = p2.panels
  p.V = p2.V
  p.vortex_points = p2.vortex_points
  p.constant_panel = p2.constant_panel
  p.old_α = p2.old_α
  p.ω = p2.ω
  p.γs = p2.γs
  p.γg = p2.γg
  p.oldγg = p2.oldγg
  p.average_length = p2.average_length
  p.dt = p2.dt
  p.T = p2.T
  p.δ = p2.δ
  p.ϵ = p2.ϵ
  p.timelapse = p2.timelapse
  p.is_init = p2.is_init
  p.θg = p2.θg
  p.θ0 = p2.θ0
  p.β = p2.β
  p.θ1 = p2.θ1
  p.θ2 = p2.θ2
  p.old_panel = p2.old_panel
  p.pivot_location = p2.pivot_location
  p.Ut = p2.Ut
  p.Vt = p2.Vt
  p.force_control = p2.force_control
  p.trailing = p2.trailing
  p.fname = p2.fname
  p.Ainv = p2.Ainv
  p.η = p2.η
  p.τ = p2.τ
  p.Tmin = p2.Tmin
  p.v_corr = p2.v_corr
  p.sheet_size = p2.sheet_size
  p.last_vp_ind = p2.last_vp_ind
  p.islumped = p2.islumped
  p.profID = p2.profID
  p.history = p2.history
  p.isgaussian = p2.isgaussian
  return nothing
end

"""
    genprofile(profile, N, position, args...)

Generate a profile with a wedged trailing edge.
"""
function genprofile(profileshape::Function, N::Int64, profID, position::Array)
  #θ = range(0, pi, length=N÷2 + 1)
  #x = @. .5cos(θ) + .5
  x = uniformdiscretization(profileshape, N)[end:-1:1]
  ye, yi = profileshape(x)

  L = []
  for i in 1:length(x)-1
    @inbounds push!(L, BuildPanel(LinearPanel, x[i] + position[1],
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

Return an array containing lengths of all panels.
"""
function getprofilelengths(p::Profile)
  λ = 2map(pr -> pr.b, p.panels)
  pushfirst!(λ, 0)
  return cumsum(λ)
end


"""
    inbodyframe(p::Profile, X, Y)

Return the vector `[X, Y]` in a frame attached to the pivot point of the body.
"""
function inbodyframe(p::Profile, X, Y)
  X0, Y0 = p.pivot_location
  x, y = rotate(X - X0, Y - Y0, p.old_α)
  return [x, y]
end

"""
    getparamsinrefframe(panel, X, Y)

Return the velocity parameters from `panel` at point `X`, `Y` in the inertial frame.
"""
function getparamsinrefframe(panel::LinearPanel, X, Y)
  u, v = getparams(panel, X, Y)
  uv1 = rotate(u[1], v[1], -(π + panel.θ))
  uv2 = rotate(u[2], v[2], -(π + panel.θ))
  U = [uv1[1], uv2[1]]
  V = [uv1[2], uv2[2]]
  return [U, V]
end

"""
    getparamsinrefframe(panel, X, Y)

Return the velocity parameters from `panel` at point `X`, `Y` in the inertial frame.
"""
function getparamsinrefframe(panel::ConstantPanel, X, Y)
  u, v = getparams(panel, X, Y)
  return rotate(u, v, -(π + panel.θ))
end

"""
    getboundpanelsinducedvel(profile, X, Y)

Return the velocity at point `X`, `Y` due to the body panels.
"""
function getboundpanelsinducedvel(p::Profile, X, Y)
  X0, Y0 = p.pivot_location
  U = V = 0
  for (i, panel) in enumerate(p.panels)
    u, v = getparamsinrefframe(panel, X, Y)
    @inbounds U += u[1] * p.γs[i] + u[2] * p.γs[i+1]
    @inbounds V += v[1] * p.γs[i] + v[2] * p.γs[i+1]
  end
  int = sum(map((pa) -> getmotionparams(pa, X, Y, p.ω, p.Ut, p.Vt, X0, Y0), p.panels))

  return [U, V] + int
end

"""
    motionupdate!(p::Profile, X, Ẋ)

Update the position, AoA and velocities of the body.
"""
function motionupdate!(p::Profile, X, Ẋ)
  X0, Y0 = p.pivot_location
  Xt, Yt, θ = X
  α = -θ
  p.Ut, p.Vt, p.ω = Ẋ
  for panel in p.panels
    translate!(panel, Xt - X0, Yt - Y0)
    setangle!(panel, α - p.old_α, [Xt, Yt])
  end
  p.old_α = α
  p.pivot_location = [Xt, Yt]
  X1 = 0.5(p.panels[1].X1 + p.panels[end].X2)
  Y1 = 0.5(p.panels[1].Y1 + p.panels[end].Y2)
  p.trailing = [X1, Y1]
  return nothing
end

"""
    globalvelocity(p::Profile, XY, incloud=false)

Return the velocity at point `XY`.

# Keyword Arguments
 - `incloud = false`: flag indicating whether `XY` is part of the cloud of point vortices.
"""
function globalvelocity(p::Profile, XY, incloud=false)
  UV = getboundpanelsinducedvel(p, XY[1], XY[2])
  UV += p.γg * getparamsinrefframe(p.constant_panel, XY[1], XY[2])
  !incloud && !isempty(p.vortex_points) &&
    (UV += sum(map(vp -> getinducedvelocity(vp, XY[1], XY[2], p.δ, p.isgaussian), p.vortex_points)))
  return UV
end

"""
    setforcecontrol!(p::Profile, fx, fy, mz)

Update forces acting on the body pivot point.
"""
function setforcecontrol!(p::Profile, fx, fy, mz)
  p.force_control = [fx, fy, -mz]
  return nothing
end

"""
    gettrailingedgevel!(p::Profile)

Return the fluid velocity at the trailing edge expressed in the inertial frame.
"""
function gettrailingedgevel!(p::Profile)
  v1 = p.γs[1]
  v2 = -p.γs[end]
  p.β = atan(tan(0.5p.θ0) * (v2 - v1) / (v2 + v1))
  p.θ1 = 0.5p.θ0 + p.β
  p.θ2 = 0.5p.θ0 - p.β
  p.θg = p.panels[1].θ - p.θ1
  u3 = 0.5 * (v1 * cos(p.θ1) + v2 * cos(p.θ2))

  return [-u3 * cos(p.θg), u3 * sin(p.θg)]
end

"""
    placeconstantpanel!(p::Profile, X2, Y2)

Place the uniform panel at location `X2`, `Y2`.
"""
function placeconstantpanel!(p::Profile, X2, Y2)
  X1 = p.trailing[1] + p.ϵ * p.average_length * cos(p.θg)
  Y1 = p.trailing[2] - p.ϵ * p.average_length * sin(p.θg)
  p.constant_panel = BuildPanel(ConstantPanel, X1, Y1, X2, Y2, p.profID)
  return nothing
end

"""
    step!(p::Profile)

Update the simulation in time.

#Keyword Arguments
 - `is4thorder = true`: indicates if a RK4 or a ForwardEuler scheme is used.
 - `need_reset = true`: indicates whether the uniform panel needs to be turned
    into a point vortex.
"""
function step!(p::Profile, accfunc, is4thorder=true, need_reset=true)
  X0, Y0 = p.pivot_location
  X1 = p.trailing[1] + p.ϵ * p.average_length * cos(p.θg)
  Y1 = p.trailing[2] - p.ϵ * p.average_length * sin(p.θg)

  function updatesystem!(p, XV, XPL, XAL, ẊAL)
    for (i, vp) ∈ enumerate(p.vortex_points)
      @inbounds vp.X = XV[i][1]
      @inbounds vp.Y = XV[i][2]
    end
    motionupdate!(p, XAL, ẊAL)
    placeconstantpanel!(p, XPL...)
    setγs!(p)
    return nothing
  end


  XY0 = map(vp -> [vp.X, vp.Y], p.vortex_points)
  XA0 = [X0, Y0, -p.old_α]
  ẊA0 = [p.Ut, p.Vt, p.ω]
  XP0 = [X1, Y1]

  ##EULER
  #kx1 = p.force_control
  #kv1 = getcloudvelocities(p.vortex_points, p.δ, p.isgaussian) .+ map(x->globalvelocity(p, x, true), XY0)
  #!isempty(p.vortex_points) && (kv1[p.last_vp_ind] += p.vortex_points[p.last_vp_ind].dv)
  #kp1 = gettrailingedgevel!(p)
  #XY = XY0 + p.dt*kv1
  #XA = XA0 + p.dt*ẊA0
  #ẊA = ẊA0 + p.dt*kx1
  #XP = XP0 + p.dt*kp1
  ##END EULER

  #RK4
  βt = 0.5p.dt
  kx1 = p.force_control
  kv1 = getcloudvelocities(p.vortex_points, p.δ, p.isgaussian) .+ map(x -> globalvelocity(p, x, true), XY0)
  !isempty(p.vortex_points) && (kv1[p.last_vp_ind] += p.vortex_points[p.last_vp_ind].dv)
  kp1 = gettrailingedgevel!(p)
  XY1 = XY0 + βt * kv1
  XA1 = XA0 + βt * ẊA0
  ẊA1 = ẊA0 + βt * kx1
  XP1 = XP0 + βt * kp1
  updatesystem!(p, XY1, XP1, XA1, ẊA1)
  setforcecontrol!(p, accfunc(p.timelapse + βt)[3]...)
  kx2 = p.force_control
  kv2 = getcloudvelocities(p.vortex_points, p.δ, p.isgaussian) .+ map(x -> globalvelocity(p, x, true), XY1)
  !isempty(p.vortex_points) && (kv2[p.last_vp_ind] += p.vortex_points[p.last_vp_ind].dv)
  kp2 = gettrailingedgevel!(p)
  XY2 = XY0 + βt * kv2
  XA2 = XA0 + βt * ẊA1
  ẊA2 = ẊA0 + βt * kx2
  XP2 = XP0 + βt * kp2
  updatesystem!(p, XY2, XP2, XA2, ẊA2)
  βt = p.dt
  kx3 = kx2
  kv3 = getcloudvelocities(p.vortex_points, p.δ, p.isgaussian) .+ map(x -> globalvelocity(p, x, true), XY2)
  !isempty(p.vortex_points) && (kv3[p.last_vp_ind] += p.vortex_points[p.last_vp_ind].dv)
  kp3 = gettrailingedgevel!(p)
  XY3 = XY0 + βt * kv3
  XA3 = XA0 + βt * ẊA2
  ẊA3 = ẊA0 + βt * kx3
  XP3 = XP0 + βt * kp3
  updatesystem!(p, XY3, XP3, XA3, ẊA3)
  setforcecontrol!(p, accfunc(p.timelapse + βt)[3]...)
  kx4 = p.force_control
  kv4 = getcloudvelocities(p.vortex_points, p.δ, p.isgaussian) .+ map(x -> globalvelocity(p, x, true), XY3)
  !isempty(p.vortex_points) && (kv4[p.last_vp_ind] += p.vortex_points[p.last_vp_ind].dv)
  kp4 = gettrailingedgevel!(p)
  coeff = p.dt / 6.0
  XY = XY0 + coeff * (kv1 + 2kv2 + 2kv3 + kv4)
  XP = XP0 + coeff * (kp1 + 2kp2 + 2kp3 + kp4)
  XA = XA0 + coeff * (ẊA0 + 2ẊA1 + 2ẊA2 + ẊA3)
  ẊA = ẊA0 + coeff * (kx1 + 2kx2 + 2kx3 + kx4)
  #END RK4

  updatesystem!(p, XY, XP, XA, ẊA)
  gettrailingedgevel!(p)

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
  X = 0.5(p.constant_panel.X1 + p.constant_panel.X2)
  Y = 0.5(p.constant_panel.Y1 + p.constant_panel.Y2)
  Γ = 2p.γg * p.constant_panel.b

  if Γ != 0
    new_vp = VortexPoint(X, Y, Γ)
    push!(p.vortex_points, new_vp)
  end

  X1 = p.trailing[1] + p.ϵ * p.average_length * cos(p.θg)
  Y1 = p.trailing[2] - p.ϵ * p.average_length * sin(p.θg)
  p.old_panel = deepcopy(p.constant_panel)
  placeconstantpanel!(p, X1, Y1)

  p.oldγg = p.γg
  p.γg = 0
  setγs!(p)
  return nothing
end

"""
    lumpvortices!(p::Profile)

Lump the last shed point vortex into the active vortex if lumping criteria are met.
Used instead of `step!` if lumping is desired.
"""
function lumpvortices!(p::Profile, accfunc)
  if p.last_vp_ind > length(p.vortex_points) - p.sheet_size + 1
    step!(p, accfunc)
    p.τ = 0
    return 0
  end

  t_index = p.last_vp_ind
  s_index = p.last_vp_ind + 1
  if sign(p.vortex_points[t_index].Γ) != sign(p.vortex_points[s_index].Γ)
    p.τ = 0
    p.vortex_points[t_index].dv = zeros(2)
    p.last_vp_ind += 1
    step!(p, accfunc)
    return 0
  end

  γ̇ = p.vortex_points[s_index].Γ / p.dt
  s_vp = p.vortex_points[s_index]
  t_vp = p.vortex_points[t_index]
  s_imp = getvorteximpulse(p, s_vp)
  t_imp = getvorteximpulse(p, t_vp)

  grads = zeros(p.N + 1, 2)
  for (i, panel) in enumerate(p.panels)
    @inbounds grads[i, :] = getvelocitygradient(t_vp, panel.Xc, panel.Yc, π + panel.θ, p.δ)
  end
  ∇γ = (p.Ainv * grads)'

  κ = map(pa -> [pa.Y1, -pa.X1], p.panels)
  push!(κ, [p.panels[end].Y2, -p.panels[end].X2])
  κ = hcat(κ...)
  A = mapslices(fook -> [mapslices(foog -> trapz(fook .* foog, getprofilelengths(p)), ∇γ, dims=2)], κ, dims=2)

  ∇p = [0 1; -1 0] + hcat(A...)'
  p.v_corr = (∇p \ (s_imp - t_imp))γ̇ / t_vp.Γ
  η = getforceerror!(p, accfunc, s_index, t_index)
  return η
end

"""
    getforceerror!(p::Profile, s_index, t_index)

Return the error in impulse due to the lumping of vortex with index `s_index`.
"""
function getforceerror!(p::Profile, accfunc, s_index, t_index)
  no_transfer_p = deepcopy(p)
  no_transfer_p.τ = 0
  no_transfer_p.vortex_points[t_index].dv = zeros(2)
  no_transfer_p.last_vp_ind += 1
  step!(no_transfer_p, accfunc)
  no_transfer_imp = getimpulse(no_transfer_p)[1:2]

  transfer!(p, s_index, t_index)
  step!(p, accfunc)
  imp = getimpulse(p)[1:2]

  η = norm(imp - no_transfer_imp) / p.dt
  if η > p.η && p.τ >= p.Tmin
    @info "No Lumping allowed"
    copier!(p, no_transfer_p)
    η = 0.0
  end
  no_transfer_p = nothing
  return η
end

"""
    getvorteximpulse(p, vp)

Return the impulse linked with vortex `vp` and its image on the body.
"""
function getvorteximpulse(p::Profile, vp)
  b = zeros(p.N + 1)
  for (m, panel) in enumerate(p.panels)
    U, V = getinducedvelocity(vp, panel.Xc, panel.Yc, p.δ, p.isgaussian)
    _, V = rotate(U, V, π + panel.θ)
    @inbounds b[m] = -V
  end
  b[end] = -vp.Γ
  γs = p.Ainv * b

  X0, Y0 = p.history.X[1], p.history.X[2]
  Xp, Yp = p.pivot_location
  panels_impulse = sum(map((pa, γL, γR) -> getimpulse(pa, γL, γR, p.ω, p.Ut, p.Vt, X0, Y0, Xp, Yp, false)[1:2], p.panels, γs[1:end-1], γs[2:end]))
  return [vp.Y - Y0, X0 - vp.X] + panels_impulse / vp.Γ
end

"""
    transfer!(p, s_index, t_index)

Lump vortex of index `s_index` into active vortex of index `t_index`.
"""
function transfer!(p::Profile, s_index, t_index)
  p.vortex_points[t_index].Γ += p.vortex_points[s_index].Γ
  p.vortex_points[t_index].dv = p.v_corr
  deleteat!(p.vortex_points, s_index)
  return nothing
end


"""
    setγs!(p::Profile)

Compute the vortex strength distribution around the profile.
"""
function setγs!(p::Profile)
  X0, Y0 = p.pivot_location
  A = zeros(p.N + 2, p.N + 2)
  b = zeros(p.N + 2)

  m::UInt64 = 1
  for panel in p.panels

    V = -0.5rotate(p.Ut - p.ω * (panel.Yc - Y0), p.Vt + p.ω * (panel.Xc - X0), π + panel.θ)[2]
    for vp in p.vortex_points
      Up, Vp = getinducedvelocity(vp, panel.Xc, panel.Yc, p.δ, p.isgaussian)
      Up, Vp = rotate(Up, Vp, π + panel.θ)
      V += Vp
    end
    Ub, Vb = sum(map((pa) -> getmotionparams(pa, panel.Xc, panel.Yc, p.ω, p.Ut, p.Vt, X0, Y0), p.panels))

    @inbounds b[m] -= V + rotate(Ub, Vb, π + panel.θ)[2]

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

  int = -sum(map(pa -> getmotionkelvin(pa, p.ω, p.Ut, p.Vt, X0, Y0), p.panels))
  @inbounds b[p.N+1] = -sum(map(vp -> vp.Γ, p.vortex_points)) + int

  #Unsteady Kutta
  @inbounds A[p.N+2, 1] = cos(p.θ1)
  @inbounds A[p.N+2, p.N+1] = cos(p.θ2)
  @inbounds A[p.N+2, p.N+2] = -1

  if p.constant_panel.b == 0
    @inbounds p.Ainv = inv(A[1:p.N+1, 1:p.N+1])
    @inbounds p.γs = p.Ainv * b[1:p.N+1]
  else
    p.Ainv = inv(A)
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

Return the linear/angular impulse exerted by the fluid on the body.
"""
function getimpulse(p::Profile)
  X0, Y0 = p.history.X[1, 1], p.history.X[1, 2]
  Xp, Yp = p.pivot_location
  panels_impulse = sum(map((pa, γL, γR) -> getimpulse(pa, γL, γR, p.ω, p.Ut, p.Vt, X0, Y0, Xp, Yp, true), p.panels, p.γs[1:end-1], p.γs[2:end]))

  if !isempty(p.vortex_points)
    points_impulse = sum(map(vp -> vp.Γ * [vp.Y - Y0, -vp.X + X0, -0.5norm2(vp.Y - Y0, vp.X - X0)], p.vortex_points))
  else
    points_impulse = zeros(3)
  end
  return panels_impulse + points_impulse
end


"""
    setcoeffimpulse!(p::Profile, idx, bdorder)

Set aerodynamic coefficients computed with impulse conservation in the flow.
"""
function setcoeffimpulse!(p::Profile, idx, bdorder)
  X0, Y0 = p.history.X[1, 1], p.history.X[1, 2]
  xyz_imp = getimpulse(p::Profile)
  oldimpulse = getlastvalues(p.history, :P, idx - 1, bdorder)

  cd, cl, cm = -2backwarddifference(xyz_imp, oldimpulse, p.dt)
  cm -= dot(p.pivot_location - [X0, Y0], [cl, -cd])

  p.history.P[idx, :] = xyz_imp
  p.history.aci[idx, :] = [cd, cl, cm]
end

"""
    setcoeffnoca!(p::Profile, idx, bdorder)

Set aerodynamic coefficients computed with a control volume approach.
"""
function setcoeffnoca!(p::Profile, idx, bdorder)
  X0, Y0 = p.history.X[1, 1], p.history.X[1, 2]
  Xp, Yp = p.pivot_location

  panels_impulse = sum(map((pa, γL, γR) -> getimpulse(pa, γL, γR, p.ω, p.Ut, p.Vt, Xp, Yp, Xp, Yp, true), p.panels, p.γs[1:end-1], p.γs[2:end]))
  oldnocaimpulse = getlastvalues(p.history, :NP, idx - 1, bdorder)

  u = 2p.old_panel.b / p.dt
  energyflux = sum(map((pa, γL, γR) -> getenergyflux(pa, γL, γR, p.ω, 0, 0, Xp, Yp, Xp, Yp), p.panels, p.γs[1:end-1], p.γs[2:end]))
  vorticityflux = p.oldγg * u * ([p.trailing[2] - Yp, Xp - p.trailing[1], -0.5norm2.((p.trailing - p.pivot_location)...)])
  Δm = sum(map(pa -> getangularimpulsecorrection(pa, p.ω, p.Ut, p.Vt, p.Ut, p.Vt, Xp, Yp, Xp, Yp), p.panels))

  f = energyflux - backwarddifference(panels_impulse, oldnocaimpulse, p.dt) - vorticityflux
  cd, cl, cm = f

  p.history.NP[idx, :] = panels_impulse
  p.history.acn[idx, :] = 2 * [cd, cl, cm + Δm]
end

"""
    getboundpotential(p::Profile)

Integrates the velocity on the profile to obtain the velocity potential `ϕ`.
"""
function getboundpotential(p::Profile)
  X0, Y0 = p.pivot_location
  ϕs = zeros(p.N + 1)
  tmp = 0
  for i in 1:p.N
    @inbounds tmp += getpotential(p.panels[i], p.γs[i], p.γs[i+1], p.ω, p.Ut, p.Vt, X0, Y0)
    @inbounds ϕs[i+1] = tmp
  end
  return ϕs
end

"""
    setcoeffpotential!(p::Profile, idx, bdorder)

Set aerodynamic coefficients using pressure integration.
"""
function setcoeffpotential!(p::Profile, idx, bdorder)
  ϕs = getboundpotential(p)
  oldϕs = getlastvalues(p.history, :ϕs, idx - 1, bdorder)

  X0, Y0 = p.pivot_location
  Us = map(pa -> p.Ut - p.ω * (pa.Y1 - Y0), p.panels)
  Vs = map(pa -> p.Vt + p.ω * (pa.X1 - X0), p.panels)
  push!(Us, p.Ut - p.ω * (p.panels[end].Y2 - Y0))
  push!(Vs, p.Vt + p.ω * (p.panels[end].X2 - X0))

  cps = 1.0 .+ norm2.(Us, Vs) .- p.γs .^ 2 .- 2backwarddifference(ϕs, oldϕs, p.dt)
  cd, cl, cm = sum(map((pa, p1, p2) -> getpressureforce(pa, p1, p2, p.pivot_location...), p.panels, cps[1:end-1], cps[2:end]))

  p.history.ϕs[idx, :] = ϕs
  p.history.cps[idx, :] = cps
  p.history.acp[idx, :] = [cd, cl, cm]
end

"""
    profilerun(p::Profile, accfunc; animate=false, is4thorder=true, bdorder=4)

Simulate the testcase.

# Arguments
 - `p`: `Profile` simulation to run.
 - `accfunc`: `Function` giving the force applied on the pivot point of the body.

# Keyword Arguments
 - `animate = false`: `Bool` indicating whether an animation has to be generated.
 - `is4thorder = true`: `Bool` indicating if the order of the Runge-Kutta integrator is 4.
                        If `false`, then a Forward Euler method is used.
 - `bdorder = 4`: Order of the backward difference schemes for time derivatives.
"""
function profilerun(p::Profile, accfunc; animate=false, is4thorder=true, bdorder=4)
  println(" Params:")
  println("---------")
  println("Panel Number     :" * string(p.N))
  println("Timestep         :" * string(p.dt))
  println("Horizon Time     :" * string(p.T))
  println("Blob δ           :" * string(p.δ))
  println("-----------------------------------------------------------")

  setγs!(p)

  count = 0
  if animate
    anim = Animation()
  end
  while (count <= 1 + Int(p.T / p.dt))
    if count % 10 == 0
      @info @sprintf "%.3f/%.3f" p.timelapse p.T
    end

    η = 0.0
    setforcecontrol!(p, accfunc(p.timelapse)[3]...)
    if p.islumped
      η, t, bytes, _, _ = @timed lumpvortices!(p, accfunc)
    else
      _, t, bytes, _, _ = @timed step!(p, accfunc, is4thorder)
    end
    count += 1

    ##Update history
    θ = p.β / p.θ0
    Γb = sum(map(vp -> vp.Γ, p.vortex_points))
    setcoeffpotential!(p, count, bdorder)
    setcoeffimpulse!(p, count, bdorder)
    setcoeffnoca!(p, count, bdorder)
    updatehistory!(p.history, count, p.timelapse,
      [p.pivot_location[1], p.pivot_location[2], p.old_α],
      [p.Ut, p.Vt, p.ω], Γb, θ)

    ##Animation
    scale = 1.2
    skip = 1
    if animate && count % skip == 0
      x = map(vp -> vp.X, p.vortex_points)
      y = map(vp -> vp.Y, p.vortex_points)
      Γ = map(vp -> vp.Γ, p.vortex_points)
      # traj = plot(p.history.X[1:count, 1], p.history.X[1:count, 2],
      # line=(1, :lightgray), aspect_ratio=1, xlim=(-0.5 - p.T, 1.2),
      # ylim=(-3, 3), framestyle=:none, grid=false, ticks=false,
      # colorbar=false, dpi=400)
      pp = scatter(x, y, zcolor=Γ, markercolor=palette(:balance), markersize=1.8,
        markerstrokewidth=0, aspect_ratio=1, xlim=(-0.5 - scale * p.T, scale),
        ylim=(-2scale, 2scale),
        framestyle=:none,
        grid=false,
        ticks=false,
        colorbar=false,
        dpi=200,
        clims=(-0.01, 0.01))
      xp = map(panel -> panel.X1, p.panels)
      yp = map(panel -> panel.Y1, p.panels)
      push!(xp, p.panels[end].X2)
      push!(yp, p.panels[end].Y2)
      plot!(pp, xp, yp, line=(1, :black), legend=false)

      frame(anim)

      ix = findfirst(x -> abs(x) < 0.5p.dt, p.timelapse .- [1, 2, 5, 10])
      istime = !isnothing(ix)
      animate && istime && savefig(p.fname * replace("$ix", "." => "") * ".svg")
    end
  end

  dumphistory(p.history, p.fname * ".txt")
  dumppressure(p.history, uniformdiscretization(p.profileshape, p.N), p.fname * "_pressure.txt")
  dumppotential(p.history, uniformdiscretization(p.profileshape, p.N), p.fname * "_potential.txt")

  animate && gif(anim, p.fname * ".gif", fps=20)

  p.timelapse = 0

  println("===================================================================")
  println("DONE")
  println("===================================================================")
end

end
