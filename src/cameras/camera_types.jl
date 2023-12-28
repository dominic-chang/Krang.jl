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
abstract type AbstractPixel end

function screen_coordinate(pix::AbstractPixel) return pix.screen_coordinate end
function metric(pix::AbstractPixel) return pix.metric end
function inclination(pix::AbstractPixel) return pix.θo end

function η(pix::AbstractPixel) 
    met = metric(pix)
    α, β = screen_coordinate(pix)
    θo = inclination(pix)
    return Krang.η(met, α, β, θo)
end

function λ(pix::AbstractPixel)
    met = metric(pix)
    α, _ = screen_coordinate(pix)
    θo = inclination(pix)
    return λ(met, α, θo)
end
function roots(pix::AbstractPixel) return get_radial_roots(metric(pix), η(pix), λ(pix)) end
function I0_inf(pix::AbstractPixel) return Ir_inf(metric(pix), roots(pix)) end
function Ir_inf(pix::AbstractPixel) return I0_inf(pix) end
function I1_inf_m_I0_terms(pix::AbstractPixel) return radial_inf_integrals(metric(pix), roots(pix))[1] end
function I2_inf_m_I0_terms(pix::AbstractPixel) return radial_inf_integrals(metric(pix), roots(pix))[2] end
function Ip_inf_m_I0_terms(pix::AbstractPixel) return radial_inf_integrals(metric(pix), roots(pix))[3] end
function Im_inf_m_I0_terms(pix::AbstractPixel) return radial_inf_integrals(metric(pix), roots(pix))[4] end
function radial_inf_integrals_m_I0_terms(pix::AbstractPixel) return radial_inf_integrals(metric(pix), roots(pix)) end
function Iϕ_inf(pix::AbstractPixel) return Iϕ_inf(metric(pix), roots(pix), λ(pix)) end
function It_inf(pix::AbstractPixel) return It_inf(metric(pix), roots(pix), λ(pix)) end
function absGθo_Gθhat(pix::AbstractPixel) return _absGθo_Gθhat(metric(pix), inclination(pix), η(pix), λ(pix)) end
function absGϕo_Gϕhat(pix::AbstractPixel) return _absGϕo_Gϕhat(metric(pix), inclination(pix), η(pix), λ(pix)) end
function absGto_Gthat(pix::AbstractPixel) return _absGto_Gthat(metric(pix), inclination(pix), η(pix), λ(pix)) end