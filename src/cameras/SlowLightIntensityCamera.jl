export SlowLightIntensityCamera

"""
    $TYPEDEF

Intensity Pixel Type.
"""
struct SlowLightIntensityPixel{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17} <: AbstractPixel
    metric::Kerr{T1}
    "Pixel screen_coordinate"
    screen_coordinate::NTuple{2,T2}
    "Radial roots"
    roots::NTuple{4,T3}
    "Radial antiderivative"
    I0_inf::T4
    "Total possible Mino time"
    total_mino_time::T5
    "Radial phi antiderivative"
    Iϕ_inf::T6
    "Radial time antiderivative"
    It_inf::T7
    I1_inf_m_I0_terms::T8
    I2_inf_m_I0_terms::T9
    Ip_inf_m_I0_terms::T10
    Im_inf_m_I0_terms::T11
    "Angular antiderivative"
    absGθo_Gθhat::NTuple{2,T12}
    "Angular ϕ antiderivative"
    absGϕo_Gϕhat::NTuple{2,T13}
    "Angular t antiderivative"
    absGto_Gthat::NTuple{2,T14}
    "Half orbit of angular t antiderivative"
    θo::T15
    η::T16
    λ::T17
    
end

@doc """
    SlowLightIntensityPixel(met::Kerr{T}, α::T, β::T, θo::T) where {T}

Construct a `SlowLightIntensityPixel` object with the given Kerr metric, screen coordinates, and inclination.

# Arguments
- `met::Kerr{T}`: The Kerr metric.
- `α::T`: The Bardeen α value (screen coordinate).
- `β::T`: The Bardeen β value (screen coordinate).
- `θo::T`: The inclination angle.

# Returns
- A `SlowLightIntensityPixel` object initialized with the given parameters.

# Details
This function calculates the η and λ values using the provided Kerr metric and screen coordinates. 
It then computes the radial roots and adjusts them if necessary. 
It also calculates the radial and angular antiderivatives. 
Finally, it initializes a `SlowLightIntensityPixel` object with the calculated values and the provided parameters.
"""
function SlowLightIntensityPixel(met::Kerr{T}, α, β, θo) where {T}
    tempη = Krang.η(met, α, β, θo)
    tempλ = Krang.λ(met, α, θo)
    roots = Krang.get_radial_roots(met, tempη, tempλ)
    numreals = sum(_isreal2.(roots))
    if (numreals == 2) && (abs(imag(roots[4]) / real(roots[4])) < eps(T))
        roots = (roots[1], roots[4], roots[2], roots[3])
    end
    I1, I2, Ip, Im = radial_inf_integrals(met, roots)
    I0_inf = Krang.Ir_inf(met, roots)
    τ_total = total_mino_time(met, roots)
    Iϕ_inf_temp = Krang.Iϕ_inf(met, roots, tempλ)
    It_inf_temp = Krang.It_inf(met, roots, tempλ)
    Gθ_Gθhat_temp = Krang._absGθo_Gθhat(met, θo, tempη, tempλ)
    Gϕ_Gϕhat_temp = Krang._absGϕo_Gϕhat(met, θo, tempη, tempλ)
    Gt_Gthat_temp = Krang._absGto_Gthat(met, θo, tempη, tempλ)

    SlowLightIntensityPixel(
        met,
        (α, β),
        roots,
        I0_inf,
        τ_total,
        Iϕ_inf_temp,
        It_inf_temp,
        I1,
        I2,
        Ip,
        Im,
        Gθ_Gθhat_temp,
        Gϕ_Gϕhat_temp,
        Gt_Gthat_temp,
        θo,
        tempη,
        tempλ,
    )
end

"""
    $TYPEDEF

Screen made of `SlowLightIntensityPixel`s.
"""
struct SlowLightIntensityScreen{A<:AbstractMatrix} <: AbstractScreen
    "Minimum and Maximum Bardeen α values"
    αrange::NTuple{2}

    "Minimum and Maximum Bardeen β values"
    βrange::NTuple{2}

    "Data type that stores screen pixel information"
    pixels::A

    SlowLightIntensityScreen{A}(αrange::NTuple{2}, βrange::NTuple{2}, pixels::A) where {A<:AbstractMatrix} =
        new{A}(αrange, βrange, pixels)

    function SlowLightIntensityScreen(met::Kerr, αmin, αmax, βmin, βmax, θo, res) 
        screen = Matrix{SlowLightIntensityPixel}(undef, res, res)
        αvals = range(αmin, αmax, length=res)
        βvals = range(βmin, βmax, length=res)
        for (iα, α) in collect(enumerate(αvals))
            for (iβ, β) in enumerate(βvals)
                screen[iα, iβ] = SlowLightIntensityPixel(met, α, β, θo)
            end
        end
        new{typeof(screen)}((αmin, αmax), (βmin, βmax), screen)
    end
end

"""
    $TYPEDEF

Camera that caches slow light raytracing information for an observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct SlowLightIntensityCamera{A} <: AbstractCamera
    metric::Kerr
    "Data type that stores screen pixel information"
    screen::SlowLightIntensityScreen{A}
    "Observer screen_coordinate"
    screen_coordinate::NTuple{2}
    @doc """
        SlowLightIntensityCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res; A=Matrix) where {T}

    Constructs a `SlowLightIntensityCamera` object.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric.
    - `θo`: The Observer's inclination angle. θo ∈ (0, π).
    - `αmin`: Minimum α coordinate on the screen.
    - `αmax`: Maximum α coordinate on the screen.
    - `βmin`: Minimum β coordinate on the screen.
    - `βmax`: Maximum β coordinate on the screen.
    - `res`: Resolution of the screen (number of pixels along one dimension).
    - `A`: Data type that stores screen pixel information (default is `Matrix`). A GPUMatrix can be used for GPU computations.

    # Returns
    - A `SlowLightIntensityCamera` object.
    """
    function SlowLightIntensityCamera(
        met::Kerr{T},
        θo,
        αmin,
        αmax,
        βmin,
        βmax,
        res
    ) where {T}
        screen = SlowLightIntensityScreen(met, αmin, αmax, βmin, βmax, θo, res)
        new{typeof(screen.pixels)}(met, screen, (T(Inf), θo))
    end
end

@inline function η(pix::SlowLightIntensityPixel)
    return pix.η
end
@inline function λ(pix::SlowLightIntensityPixel)
    return pix.λ
end
@inline function roots(pix::SlowLightIntensityPixel)
    return pix.roots
end
@inline function screen_coordinate(pix::SlowLightIntensityPixel)
    return pix.screen_coordinate
end
@inline function inclination(pix::SlowLightIntensityPixel)
    return pix.θo
end
function I0_inf(pix::SlowLightIntensityPixel)
    return pix.I0_inf
end
function total_mino_time(pix::SlowLightIntensityPixel)
    return pix.total_mino_time
end
function Ir_inf(pix::SlowLightIntensityPixel)
    return pix.I0_inf
end
function I1_inf_m_I0_terms(pix::SlowLightIntensityPixel)
    return pix.I1_inf_m_I0_terms
end
function I2_inf_m_I0_terms(pix::SlowLightIntensityPixel)
    return pix.I2_inf_m_I0_terms
end
function Ip_inf_m_I0_terms(pix::SlowLightIntensityPixel)
    return pix.Ip_inf_m_I0_terms
end
function Im_inf_m_I0_terms(pix::SlowLightIntensityPixel)
    return pix.Im_inf_m_I0_terms
end
function radial_inf_integrals_m_I0_terms(pix::SlowLightIntensityPixel)
    return I1_inf_m_I0_terms(pix),
    I2_inf_m_I0_terms(pix),
    Ip_inf_m_I0_terms(pix),
    Im_inf_m_I0_terms(pix)
end
function Iϕ_inf(pix::SlowLightIntensityPixel)
    return pix.Iϕ_inf
end
function It_inf(pix::SlowLightIntensityPixel)
    return pix.It_inf
end
function absGθo_Gθhat(pix::SlowLightIntensityPixel)
    return pix.absGθo_Gθhat
end
function absGϕo_Gϕhat(pix::SlowLightIntensityPixel)
    return pix.absGϕo_Gϕhat
end
function absGto_Gthat(pix::SlowLightIntensityPixel)
    return pix.absGto_Gthat
end
