export p_bl_d, penrose_walker, screen_polarisation, PowerLaw, evpa, jac_zamo_d_bl_u, jac_bl_d_zamo_u, jac_zamo_u_bl_d, jac_bl_u_zamo_d, jac_fluid_u_zamo_d
##----------------------------------------------------------------------------------------------------------------------
#Polarization stuff
##----------------------------------------------------------------------------------------------------------------------
"""
    Returns the momentum form in the Boyer-Lindquist basis.
"""
function p_bl_d(metric::Kerr{T}, r, θ, η, λ, νr::Bool, νθ::Bool) where {T}
    return @SVector[
        -one(T),
        (νr ? one(T) : -one(T))* √max(zero(T), r_potential(metric, η, λ, r)) / Δ(metric, r),
        (νθ ? one(T) : -one(T))* √max(zero(T), θ_potential(metric, η, λ, θ)),
        λ
    ]
end

"""
Returns the Jacobian which converts a Boyer-Lindquist covector to ZAMO basis.
"""
function jac_zamo_d_bl_u(metric::Kerr{T}, r, θ) where {T}
    a = metric.spin
    # coords = {t, r, θ, ϕ}
    Σt = Σ(metric, r, θ)
    Δt = Δ(metric, r)
    At = A(metric, r, θ)

    return @SMatrix [# Eq 3.1 1972ApJ...178..347B
        √(At / (Σt * Δt)) zero(T) zero(T) 2*a*r/√(At * Δt * Σt)
        zero(T) √(Δt / Σt) zero(T) zero(T)
        zero(T) zero(T) zero(T) √(Σt / At)*csc(θ)
        zero(T) zero(T) -inv(√Σt) zero(T)
    ]
end

"""
Jacobian which converts ZAMO covector to a Boyer-Lindquist basis
"""
function jac_bl_d_zamo_u(metric::Kerr{T}, r, θ) where {T}
    a = metric.spin
    Σt = Σ(metric, r, θ)
    Δt = Δ(metric, r)
    At = A(metric, r, θ)

    return @SMatrix [
        # coords = {t, r, ϕ, θ}
        √((Σt * Δt) / At) zero(T) -2*a*r*sin(θ)/√(At * Σt) zero(T)
        zero(T) √(Σt / Δt) zero(T) zero(T)
        zero(T) zero(T) zero(T) -√Σt
        zero(T) zero(T) √(At / Σt)*sin(θ) zero(T)
    ]
end

"""
Jacobian which converts Boyer-Lindquist vector to a ZAMO basis
"""
function jac_zamo_u_bl_d(metric::Kerr{T}, r, θ) where {T}
    a = metric.spin
    Σt = Σ(metric, r, θ)
    Δt = Δ(metric, r)
    At = A(metric, r, θ)

    return @SMatrix [#  Eq 3.2 1972ApJ...178..347B
        # coords = {t, r, ϕ, θ}
        √((Σt * Δt) / At) zero(T) zero(T) zero(T)
        zero(T) √(Σt / Δt) zero(T) zero(T)
        -(2a * r * sin(θ))/√(At * Σt) zero(T) zero(T) √(At / Σt)*sin(θ)
        zero(T) zero(T) -√Σt zero(T)
    ]
end

"""
Jacobian which converts ZAMO vector to a Boyer-Lindquist basis
"""
function jac_bl_u_zamo_d(metric::Kerr{T}, r, θ) where {T}
    a = metric.spin
    Σt = Σ(metric, r, θ)
    Δt = Δ(metric, r)
    At = A(metric, r, θ)

    return @SMatrix [
        # coords = {t, r, ϕ, θ}
        √(At / (Σt * Δt)) zero(T) zero(T) zero(T)
        zero(T) √(Δt / Σt) zero(T) zero(T)
        zero(T) zero(T) zero(T) -inv(√Σt)
        2a*r/√(Σt * Δt * At) zero(T) √(Σt / At)*csc(θ) zero(T)
    ]
end

"""
Jacobian which expreases ZAMO vector in the fluid frame
"""
function jac_fluid_u_zamo_d(::Kerr{T}, β, θ, φ) where {T}
    γ = inv(√(1 - β^2))
    sinφ = sin(φ)
    cosφ = cos(φ)
    sinθ = sin(θ)
    cosθ = cos(θ)

    return @SMatrix [
        γ -β*γ*cosφ*sinθ -β*γ*sinφ*sinθ -β*γ*cosθ
        -β*γ*cosφ*sinθ cosθ^2*cosφ^2+γ*cosφ^2*sinθ^2+sinφ^2 (γ-1)*cosφ*sinθ^2*sinφ (γ-1)*cosθ*cosφ*sinθ
        -β*γ*sinθ*sinφ (γ-1)*cosφ*sinθ^2*sinφ cosφ^2+(cosθ^2+γ*sinθ^2)*sinφ^2 (γ-1)*cosθ*sinθ*sinφ
        -β*γ*cosθ (γ-1)*cosθ*cosφ*sinθ (γ-1)*cosθ*sinθ*sinφ γ*cosθ^2+sinθ^2
    ]
end

"""
Returns the Penrose walker constant for a photon with momentum p_u emitted from a fluid particle with momentum f_u.
"""
function penrose_walker(metric::Kerr{T}, r, θ, p_u::AbstractVector, f_u::AbstractVector) where {T}# Eq 6 arXiv:2001.08750v1
    a = metric.spin
    pt, pr, pϕ, pθ = p_u
    ft, fr, fϕ, fθ = f_u
    sinθ = sin(θ)
    cosθ = cos(θ)

    A = pt * fr - pr * ft + a * sinθ * sinθ * (pr * fϕ - pϕ * fr)
    B = ((r * r + a * a) * (pϕ * fθ - pθ * fϕ) - a * (pt * fθ - pθ * ft)) * sinθ
    return A * r - B * a * cosθ, -(A * a * cosθ - B * r)
end

