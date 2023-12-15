using Metal

export BasicCamera
"""
    $TYPEDEF

Basic Pixel Type.
"""
struct BasicGPUPixel{T} <: AbstractPixel
    metric::Kerr{T}
    screen_coordinate::NTuple{2, T}
    "Radial roots"
    roots::NTuple{4,Complex{T}}
    θo::T
    η::T
    λ::T
    function BasicGPUPixel(met::Kerr{T}, α, β, θo) where {T}
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
            θo, tempη, tempλ
        )
    end
end

"""
    $TYPEDEF

Screen made of Basic Pixels.
"""
struct BasicGPUScreen{T} <: AbstractScreen
    "Minimum and Maximum Bardeen α values"
    αrange::NTuple{2, T}

    "Minimum and Maximum Bardeen β values"
    βrange::NTuple{2, T}

    "Data type that stores screen pixel information"
    pixels::MtlMatrix{BasicGPUPixel{T}}
    function BasicGPUScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res) where {T}
        screen = Matrix{BasicGPUPixel}(undef, res, res)
        αvals = range(αmin, αmax, length=res)
        βvals = range(βmin, βmax, length=res)
        
        for (iα, α) in enumerate(αvals)
            for (iβ, β) in enumerate(βvals)
                screen[iα, iβ] = BasicGPUPixel(met, α, β, θo)
            end
        end
        new{T}((αmin, αmax), (βmin, βmax), MtlArray([i for i in screen]))
    end
end

"""
    $TYPEDEF

Observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct BasicGPUCamera{T} <: AbstractCamera
    metric::Kerr{T}
    "Data type that stores screen pixel information"
    screen::BasicGPUScreen{T}
    "Observer screen_coordinate"
    screen_coordinate::NTuple{2, T}
    function BasicGPUCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res) where {T}
        new{T}(met, BasicGPUScreen(met, αmin, αmax, βmin, βmax, θo, res), (T(Inf), θo))
    end
end

function I0_inf(pix::BasicGPUPixel) 
    return Krang.Ir_inf(pix.metric, pix.roots)
end

function Ir_inf(pix::BasicGPUPixel) 
    return Krang.Ir_inf(pix.metric, pix.roots)
end

function Iϕ_inf(pix::BasicGPUPixel)
    λtemp = λ(pix.metric, pix.screen_coordinate[1], pix.θo)
    return Iϕ_inf(pix.metric, pix.roots, λtemp)
end

function It_inf(pix::BasicGPUPixel)
    λtemp = λ(pix.metric, pix.screen_coordinate[1], pix.θo)
    return It_inf(pix.metric, pix.roots, λtemp)
end

function radial_inf_integrals_m_I0_terms(pix::BasicGPUPixel)
    return radial_inf_integrals(metric(pix),roots(pix))
end

function absGθo_Gθhat(pix::BasicGPUPixel)
    @unpack metric, screen_coordinate, θo =  pix
    α, β = screen_coordinate
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    return _absGθo_Gθhat(metric, θo, ηtemp, λtemp)
end
function absGϕo_Gϕhat(pix::BasicGPUPixel)
    @unpack metric, screen_coordinate, θo =  pix
    α, β = screen_coordinate
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    return _absGϕo_Gϕhat(metric, θo, ηtemp, λtemp)
end
function absGto_Gthat(pix::BasicGPUPixel)
    @unpack metric, screen_coordinate, θo =  pix
    α, β = screen_coordinate
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    return _absGto_Gthat(metric, θo, ηtemp, λtemp)
end