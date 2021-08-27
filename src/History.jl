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
    Γ::Array{Float64, 1}
    θ::Array{Float64, 1}
end

function History(N, dt, T)
    n = 2 + Int(T/dt)
    return History(zeros(n), zeros(n, 3), zeros(n, 3), zeros(n, 3), zeros(n, 3), zeros(n, 3),
                   zeros(n, N+1), zeros(n, N+1), zeros(n), zeros(n))
end

"""
    updatehistory!(h::History, i, t, X, Ẋ, acp, acn, aci, ϕs, cps, Γ, θ)

Update history of the airfoil.
"""
function updatehistory!(h::History, i, t, X, Ẋ, acp, acn, aci, ϕs, cps, Γ, θ)
    h.t[i] = t
    h.X[i, :] = X
    h.Ẋ[i, :] = Ẋ
    h.acp[i, :] = acp
    h.acn[i, :] = acn
    h.aci[i, :] = aci
    h.ϕs[i, :] = ϕs
    h.cps[i, :] = cps
    h.Γ[i] = Γ
    h.θ[i] = θ
end
