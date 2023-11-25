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

function η(pix::AbstractPixel) return pix.η end
function λ(pix::AbstractPixel) return pix.λ end
function roots(pix::AbstractPixel) return pix.roots end
function screen_coordinate(pix::AbstractPixel) return pix.screen_coordinate end
function metric(pix::AbstractPixel) return pix.metric end
function inclination(pix::AbstractPixel) return pix.θo end
function I0_inf(pix::AbstractPixel) return pix.I0_inf end
function Ir_inf(pix::AbstractPixel) return pix.I0_inf end
function I1_inf_m_I0_terms(pix::AbstractPixel) return pix.I1_inf_m_I0_terms end
function I2_inf_m_I0_terms(pix::AbstractPixel) return pix.I2_inf_m_I0_terms end
function Ip_inf_m_I0_terms(pix::AbstractPixel) return pix.Ip_inf_m_I0_terms end
function Im_inf_m_I0_terms(pix::AbstractPixel) return pix.Im_inf_m_I0_terms end
function radial_inf_integrals_m_I0_terms(pix::AbstractPixel) return I1_inf_m_I0_terms(pix), I2_inf_m_I0_terms(pix), Ip_inf_m_I0_terms(pix), Im_inf_m_I0_terms(pix) end
function Iϕ_inf(pix::AbstractPixel) return pix.Iϕ_inf end
function It_inf(pix::AbstractPixel) return pix.It_inf end
function absGθo_Gθhat(pix::AbstractPixel) return pix.absGθo_Gθhat end
function absGϕo_Gϕhat(pix::AbstractPixel) return pix.absGϕo_Gϕhat end
function absGto_Gthat(pix::AbstractPixel) return pix.absGto_Gthat end





