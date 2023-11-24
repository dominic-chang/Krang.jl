export BasicCamera
"""
    $TYPEDEF

Basic Pixel Type.
"""
struct BasicPixel{T} <: AbstractPixel
    metric::Kerr{T}
    screen_coordinate::NTuple{2, T}
    "Radial roots"
    roots::NTuple{4,Complex{T}}
    "Radial antiderivative"
    I0_inf::T
    "Angular antiderivative"
    absGθ::T
    "Half orbit of angular antiderivative"
    Ghat::T
    "Inclination"
    θo::T
    function BasicPixel(met::Kerr{T}, α, β, θo) where {T}
        tempη = Krang.η(met, α, β, θo)
        tempλ = Krang.λ(met, α, θo)
        roots = Krang.get_radial_roots(met, tempη, tempλ)
        numreals = sum(_isreal2.(roots))
        if (numreals == 2) && (abs(imag(roots[4])) < sqrt(eps(T)))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        new{T}(
            met,
            (α, β), 
            roots,
            Krang.Ir_inf(met, roots), 
            Krang._absGθo_Gθhat(met, θo, tempη, tempλ)..., 
            θo
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
struct BasicCamera{T} <: AbstractCamera
    metric::Kerr{T}
    "Data type that stores screen pixel information"
    screen::BasicScreen{T}
    "Observer screen_coordinate"
    screen_coordinate::NTuple{2, T}
    function BasicCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res) where {T}
        new{T}(met, BasicScreen(met, αmin, αmax, βmin, βmax, θo, res), (T(Inf), θo))
    end
end

function Iϕ_inf(pix::BasicPixel)
    λtemp = λ(pix.metric, pix.screen_coordinate[1], pix.θo)
    return Iϕ_inf(pix.metric, pix.roots, λtemp)
end

function It_inf(pix::BasicPixel)
    λtemp = λ(pix.metric, pix.screen_coordinate[1], pix.θo)
    return It_inf(pix.metric, pix.roots, λtemp)
end

function radial_inf_integrals_m_I0_terms(pix::BasicPixel)
    return radial_inf_integrals(metric(pix),roots(pix))
end
function absGθo_Gθhat(pix::BasicPixel)
    @unpack metric, screen_coordinate, θo =  pix
    α, β = screen_coordinate
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    return _absGθo_Gθhat(metric, θo, ηtemp, λtemp)
end
function absGϕo_Gϕhat(pix::BasicPixel)
    @unpack metric, screen_coordinate, θo =  pix
    α, β = screen_coordinate
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    return _absGϕo_Gϕhat(metric, θo, ηtemp, λtemp)
end
function absGto_Gthat(pix::BasicPixel)
    @unpack metric, screen_coordinate, θo =  pix
    α, β = screen_coordinate
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    return _absGto_Gthat(metric, θo, ηtemp, λtemp)
end