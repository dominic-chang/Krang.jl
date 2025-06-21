export p_bl_d,
    penrose_walker,
    screen_polarization,
    PowerLaw,
    evpa
    
##----------------------------------------------------------------------------------------------------------------------
#Polarization stuff
##----------------------------------------------------------------------------------------------------------------------

"""
Returns the Penrose walker constant for a photon with momentum p_u emitted from a fluid particle with momentum f_u.
"""
function penrose_walker(
    metric::Kerr{T},
    r,
    θ,
    p_u::AbstractVector,
    f_u::AbstractVector,
) where {T}# Eq 6 arXiv:2001.08750v1
    a = metric.spin
    pt, pr, pϕ, pθ = p_u
    ft, fr, fϕ, fθ = f_u
    sinθ = sin(θ)
    cosθ = cos(θ)

    A = pt * fr - pr * ft + a * sinθ * sinθ * (pr * fϕ - pϕ * fr)
    B = ((r * r + a * a) * (pϕ * fθ - pθ * fϕ) - a * (pt * fθ - pθ * ft)) * sinθ
    return A * r - B * a * cosθ, -(A * a * cosθ + B * r)
end
