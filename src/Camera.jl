export AssymptoticSlowLightCamera

"""
    $TYPEDEF

Abstract Observer Type
"""
abstract type AbstractObserver end

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

"""
    $TYPEDEF

Basic Pixel Type.
"""
struct BasicPixel{T} <: AbstractPixel
    "Pixel location"
    location::NTuple{2, T}
    "Radial roots"
    roots::NTuple{4,Complex{T}}
    "Radial antiderivative"
    Ir::T
    "Radial phi antiderivative"
    Iϕ::T
    "Radial time antiderivative"
    It::T
    "Angular antiderivative"
    absGθ::T
    "Half orbit of angular antiderivative"
    Ghat::T
    "Angular ϕ antiderivative"
    absGϕ::T
    "Half orbit of angular ϕ antiderivative"
    Gϕhat::T
    "Angular t antiderivative"
    absGt::T
    "Half orbit of angular t antiderivative"
    Gthat::T
    I1o_m_I0_terms::T
    I2o_m_I0_terms::T
    Ipo_m_I0_terms::T
    Imo_m_I0_terms::T
    function BasicPixel(met::Kerr{T}, α, β, θo) where {T}
        tempη = Krang.η(met, α, β, θo)
        tempλ = Krang.λ(met, α, θo)
        roots = Krang.get_radial_roots(met, tempη, tempλ)
        numreals = sum(_isreal2.(roots))
        if (numreals == 2) && (abs(imag(roots[4])) < sqrt(eps(T)))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        I1, I2, Ip, Im = radial_inf_integrals(met, roots)
        new{T}(
            (α, β), 
            roots, 
            Krang.Ir_inf(met, roots), 
            Krang.Iϕ_inf(met, roots, tempλ), 
            Krang.It_inf(met, roots, tempλ), 
            Krang._absGθo_Gθhat(met, θo, tempη, tempλ)..., 
            Krang._absGϕo_Gϕhat(met, θo, tempη, tempλ)..., 
            Krang._absGto_Gthat(met, θo, tempη, tempλ)...,
            I1, I2, Ip, Im
        )
    end
end

"""
    $TYPEDEF

Screen made of Basic Pixels.
"""
struct BasicScreen{T} <: AbstractScreen
    "Minimum and Maximum Bardeen α values"
    αrange::NTuple{2, T}

    "Minimum and Maximum Bardeen β values"
    βrange::NTuple{2, T}

    "Data type that stores screen pixel information"
    pixels::Matrix{BasicPixel{T}}
    function BasicScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res) where {T}
        screen = Matrix{BasicPixel}(undef, res, res)
        αvals = range(αmin, αmax, length=res)
        βvals = range(βmin, βmax, length=res)
        
        for (iα, α) in enumerate(αvals)
            for (iβ, β) in enumerate(βvals)
                screen[iα, iβ] = BasicPixel(met, α, β, θo)
            end
        end
        new{T}((αmin, αmax), (βmin, βmax), screen)
    end
end

"""
    $TYPEDEF

Observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct AssymptoticSlowLightCamera{T} <: AbstractObserver
    metric::Kerr{T}
    "Data type that stores screen pixel information"
    screen::BasicScreen{T}
    "Observer location"
    location::SVector{2, T}
    function AssymptoticSlowLightCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res) where {T}
        new{T}(met, BasicScreen(met, αmin, αmax, βmin, βmax, θo, res), SVector{2, T}(T(Inf), θo))
    end
end