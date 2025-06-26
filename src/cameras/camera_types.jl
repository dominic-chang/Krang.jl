"""
    $TYPEDEF

Abstract Observer Type
"""
abstract type AbstractCamera end

"""
    $TYPEDEF

Abstract Screen Type
"""
abstract type AbstractScreen end

"""
    $TYPEDEF

Abstract Pixel Type
"""
abstract type AbstractPixel{T} end

"""
    screen_coordinate(pix::AbstractPixel)

Get the screen coordinates of a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel for which to get the screen coordinates.

# Returns
- The screen coordinates of the pixel.
"""
function screen_coordinate(pix::AbstractPixel)
    return pix.screen_coordinate
end

"""
    metric(pix::AbstractPixel)

Get the space time metric.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The space time metric
"""
function metric(pix::AbstractPixel)
    return pix.metric
end

"""
    inclination(pix::AbstractPixel)

Get the inclination angle (θo) of a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The inclination angle (θo) of the pixel.
"""
function inclination(pix::AbstractPixel)
    return pix.θo
end

"""
    η(pix::AbstractPixel)

Calculate the η value for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The η value of the pixel.
"""
function η(pix::AbstractPixel)
    @warn "η is not stored in pixel, calculating η under infinite observer assumption instead"

    met = metric(pix)
    α, β = screen_coordinate(pix)
    θo = inclination(pix)
    return Krang.η(met, α, β, θo)
end

"""
    λ(pix::AbstractPixel)

Calculate the λ value for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The λ value of the pixel.
"""
function λ(pix::AbstractPixel)
    @warn "λ is not stored in pixel, calculating λ under infinite observer assumption instead"
    met = metric(pix)
    α, _ = screen_coordinate(pix)
    θo = inclination(pix)
    return λ(met, α, θo)
end

"""
    roots(pix::AbstractPixel)

Calculate the radial roots for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The radial roots of the pixel.
"""
function roots(pix::AbstractPixel)
    return get_radial_roots(metric(pix), η(pix), λ(pix))
end

"""
    I0_inf(pix::AbstractPixel)

Calculate the I0 infinity value for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The I0 infinity value of the pixel.
"""
function I0_inf(pix::AbstractPixel)
    return Ir_inf(metric(pix), roots(pix))
end


"""
    I0_o(pix::AbstractPixel)

Calculate the I0 anti-derivative value for a given pixel at radius ro.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The I0 anti-derivative of the pixel.
"""
function I0_o(pix::AbstractPixel)
    if isinf(pix.ro)
        return I0_inf(pix)
    end
    return Ir_s(metric(pix), pix.ro, roots(pix), true)
end

"""
    total_mino_time(pix::AbstractPixel)

Return the total possible Mino time for a ray associated with a pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The total possible Mino time for a ray associated with the pixel.
"""
function total_mino_time(pix::AbstractPixel)
    return total_mino_time(metric(pix), roots(pix))
end
"""
    Ir_inf(pix::AbstractPixel)

Calculate the Ir infinity value for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The Ir infinity value of the pixel.
"""
function Ir_inf(pix::AbstractPixel)
    return I0_inf(pix)
end

"""
    I1_inf_m_I0_terms(pix::AbstractPixel)

Calculate the I1 infinity minus I0 terms for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The I1 infinity minus I0 terms of the pixel.
"""
function I1_inf_m_I0_terms(pix::AbstractPixel)
    @warn "I1_inf_m_I0_terms are not stored in pixel, calculating I1_inf_integrals under infinite observer assumption instead"
    return radial_inf_integrals(metric(pix), roots(pix))[1]
end

"""
    I1_o_m_I0_terms(pix::AbstractPixel)

Calculate the I1 anti-derivative minus I0 terms for a given pixel at radius ro.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The I1 anti-derivative minus I0 terms of the pixel.
"""
function I1_o_m_I0_terms(pix::AbstractPixel)
    @warn "I1_o_m_I0_terms are not stored in pixel, calculating I1_o_integrals under infinite observer assumption instead"
    return radial_o_integrals(metric(pix), roots(pix))[1]
end


"""
    I2_inf_m_I0_terms(pix::AbstractPixel)

Calculate the I2 infinity minus I0 terms for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The I2 infinity minus I0 terms of the pixel.
"""
function I2_inf_m_I0_terms(pix::AbstractPixel)
    @warn "I2_inf_m_I0_terms are not stored in pixel, calculating I2_inf_integrals under infinite observer assumption instead"
    return radial_inf_integrals(metric(pix), roots(pix))[2]
end

"""
    I2_o_m_I0_terms(pix::AbstractPixel)

Calculate the I2 anti-derivative minus I0 terms for a given pixel at radius ro.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The I2 anti-derivative minus I0 terms of the pixel.
"""
function I2_o_m_I0_terms(pix::AbstractPixel)
    @warn "I2_o_m_I0_terms are not stored in pixel, calculating I2_o_integrals under infinite observer assumption instead"
    return radial_o_integrals(metric(pix), roots(pix))[2]
end


"""
    Ip_inf_m_I0_terms(pix::AbstractPixel)

Calculate the Ip infinity minus I0 terms for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The Ip infinity minus I0 terms of the pixel.
"""
function Ip_inf_m_I0_terms(pix::AbstractPixel)
    return radial_inf_integrals(metric(pix), roots(pix))[3]
end

"""
    Ip_o_m_I0_terms(pix::AbstractPixel)

Calculate the Ip anti-derivative minus I0 terms for a given pixel at radius ro.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The Ip anti-derivative minus I0 terms of the pixel.
"""
function Ip_o_m_I0_terms(pix::AbstractPixel)
    return radial_o_integrals(metric(pix), roots(pix))[3]
end


"""
    Im_inf_m_I0_terms(pix::AbstractPixel)

Calculate the Im infinity minus I0 terms for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The Im infinity minus I0 terms of the pixel.
"""
function Im_inf_m_I0_terms(pix::AbstractPixel)
    return radial_inf_integrals(metric(pix), roots(pix))[4]
end

"""
    Im_o_m_I0_terms(pix::AbstractPixel)

Calculate the Im anti-derivative minus I0 terms for a given pixel at radius ro.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The Im anti-derivative minus I0 terms of the pixel.
"""
function Im_o_m_I0_terms(pix::AbstractPixel)
    return radial_o_integrals(metric(pix), roots(pix))[4]
end


"""
    radial_inf_integrals_m_I0_terms(pix::AbstractPixel)

Calculate the radial infinity minus I0 terms for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The radial infinity minus I0 terms of the pixel.
"""
function radial_inf_integrals_m_I0_terms(pix::AbstractPixel)
    return radial_inf_integrals(metric(pix), roots(pix))
end

"""
    radial_o_integrals(pix::AbstractPixel)

Calculate the radial anti-derivative minus I0 terms for a given pixel at radius ro.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The radial anti-derivative minus I0 terms of the pixel.
"""
function radial_o_integrals(pix::AbstractPixel)
    return radial_o_integrals_m_I0_terms(metric(pix), roots(pix))
end


"""
    Iϕ_inf(pix::AbstractPixel)

Calculate the Iϕ infinity terms for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The Iϕ infinity terms of the pixel.
"""
function Iϕ_inf(pix::AbstractPixel)
    return Iϕ_inf(metric(pix), roots(pix), λ(pix))
end

"""
    It_inf(pix::AbstractPixel)

Calculate the It infinity terms for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The It infinity terms of the pixel.
"""
function It_inf(pix::AbstractPixel)
    return It_inf(metric(pix), roots(pix), λ(pix))
end

"""
    absGθo_Gθhat(pix::AbstractPixel)

Calculate the absolute value of Gθo divided by Gθhat for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The absolute value of Gθo divided by Gθhat of the pixel.
"""
function absGθo_Gθhat(pix::AbstractPixel)
    return _absGθo_Gθhat(metric(pix), inclination(pix), η(pix), λ(pix))
end

"""
    absGϕo_Gϕhat(pix::AbstractPixel)

Calculate the absolute value of Gϕo divided by Gϕhat for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The absolute value of Gϕo divided by Gϕhat of the pixel.
"""
function absGϕo_Gϕhat(pix::AbstractPixel)
    return _absGϕo_Gϕhat(metric(pix), inclination(pix), η(pix), λ(pix))
end

"""
    absGto_Gthat(pix::AbstractPixel)

Calculate the absolute value of Gto divided by Gthat for a given pixel.

# Arguments
- `pix::AbstractPixel`: The pixel of a screen.

# Returns
- The absolute value of Gto divided by Gthat of the pixel.
"""
function absGto_Gthat(pix::AbstractPixel)
    return _absGto_Gthat(metric(pix), inclination(pix), η(pix), λ(pix))
end
