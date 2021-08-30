#!/usr/bin/env julia
# File              : src/History.jl
# Author            : Denis Dumoulin <denis.dumoulin@uclouvain.be>
# Date              : 30.08.2021
# Last Modified Date: 30.08.2021
"""
    History

A structure allowing to record useful data.

# Fields
 - `t`: timelapse.
 - `X`: position `[X, Y, α]`.
 - `Ẋ`: velocity `[U, V, ω]`.
 - `acp`: aerodynamic coefficients obtained with pressure integration.
 - `acn`: aerodynamic coefficients obtained with control volume approach.
 - `aci`: aerodynamic coefficients obtained with momentum conservation.
 - `ϕs`: velocity potential on body surface.
 - `cps`: pressure coefficients on body surface.
 - `P`: impulse `[px, py, pm]`.
 - `Γ`: total body circulation.
 - `θ`: shedding angle (0 is the bissector).
"""
mutable struct History
    t::Array{Float64, 1}
    X::Array{Float64, 2}
    Ẋ::Array{Float64, 2}
    acp::Array{Float64, 2}
    acn::Array{Float64, 2}
    aci::Array{Float64, 2}
    ϕs::Array{Float64, 2}
    cps::Array{Float64, 2}
    P::Array{Float64, 2}
    NP::Array{Float64, 2}
    Γ::Array{Float64, 1}
    θ::Array{Float64, 1}
end

function History(N, dt, T)
    n = 2 + Int(T/dt)
    return History(zeros(n), zeros(n, 3), zeros(n, 3), zeros(n, 3), zeros(n, 3), zeros(n, 3),
                   zeros(n, N+1), zeros(n, N+1), zeros(n, 3), zeros(n, 3), zeros(n), zeros(n))
end

"""
    updatehistory!(h::History, i, t, X, Ẋ, acp, acn, aci, ϕs, cps, Γ, θ)

Update history of the airfoil.
"""
function updatehistory!(h::History, i, t, X, Ẋ, Γ, θ)
    h.t[i] = t
    h.X[i, :] = X
    h.Ẋ[i, :] = Ẋ
    h.Γ[i] = Γ
    h.θ[i] = θ
end

"""
    updatehistory!(h::History, sym::Symbol, i, val)

Update field `sym` of the history.
"""
function updatehistory!(h::History, sym::Symbol, i, val)
    v = getfield(h, sym)
    v[idx, :] = val
    setfield!(h, sym, v)
end

"""
    getlastvalues(h::History, sym::Symbol, i, n)

Return the `n` last values of `h.sym` up to index `i`.
Values with negative indices are replaced with zeros.
"""
function getlastvalues(h::History, sym::Symbol, i, n)
    v = getfield(h, sym)
    i < n &&  return [circshift(v, (-i, 0))[end-n+1:end, :][j, :] for j in n:-1:1]
    return [v[i-n+1:i, :][j, :] for j in n:-1:1]
end

"""
    plotaero(h::history, save=false)
"""
function plotaero(h::History, save=false, ylim = (-.5, 1))
    p1 = plot(h.t[6:end], h.acn[6:end, 1])
    plot!(p1, h.t[6:end], h.aci[6:end, 1])
    plot!(p1, h.t[6:end], h.acp[6:end, 1],
          title = "Cd", legend = false)

    p2 = plot(h.t[6:end], h.acn[6:end, 2])
    plot!(p2, h.t[6:end], h.aci[6:end, 2])
    plot!(p2, h.t[6:end], h.acp[6:end, 2],
          title = "Cl", legend = false)

    p3 = plot(h.t[6:end], h.acn[6:end, 3], label="Control volume")
    plot!(p3, h.t[6:end], h.aci[6:end, 3], label="Impulse conservation")
    plot!(p3, h.t[6:end], h.acp[6:end, 3],
          title = "Cm", label="Cp integration")

    plot(p1, p2, p3, layout = (1, 3), size=(900, 300))
end
