export SlowLightIntensityCamera

"""
    $TYPEDEF

Intensity Pixel Type.
"""
struct SlowLightIntensityPixel{T} <: AbstractPixel{T}
    metric::Kerr{T}
    "Pixel screen_coordinate"
    screen_coordinate::NTuple{2, T}
    "Radial roots"
    roots::NTuple{4,Complex{T}}
    "Radial antiderivative"
    I0_inf::T
    "Radial phi antiderivative"
    Iϕ_inf::T
    "Radial time antiderivative"
    It_inf::T
    I1_inf_m_I0_terms::T
    I2_inf_m_I0_terms::T
    Ip_inf_m_I0_terms::T
    Im_inf_m_I0_terms::T
    "Angular antiderivative"
    absGθo_Gθhat::NTuple{2,T}
    "Angular ϕ antiderivative"
    absGϕo_Gϕhat::NTuple{2,T}
    "Angular t antiderivative"
    absGto_Gthat::NTuple{2,T}
    "Half orbit of angular t antiderivative"
    θo::T
    η::T
    λ::T
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
    function SlowLightIntensityPixel(met::Kerr{T}, α::T, β::T, θo) where {T}
        tempη = Krang.η(met, α, β, θo)
        tempλ = Krang.λ(met, α, θo)
        roots = Krang.get_radial_roots(met, tempη, tempλ)
        numreals = sum(_isreal2.(roots))
        if (numreals == 2) && (abs(imag(roots[4])) < sqrt(eps(T)))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        I1, I2, Ip, Im = radial_inf_integrals(met, roots)
        new{T}(
            met,
            (α, β), 
            roots, 
            Krang.Ir_inf(met, roots), 
            Krang.Iϕ_inf(met, roots, tempλ), 
            Krang.It_inf(met, roots, tempλ), 
            I1, I2, Ip, Im, 
            Krang._absGθo_Gθhat(met, θo, tempη, tempλ), 
            Krang._absGϕo_Gϕhat(met, θo, tempη, tempλ), 
            Krang._absGto_Gthat(met, θo, tempη, tempλ),
            θo, tempη, tempλ
        )
    end
end

"""
    $TYPEDEF

Screen made of `SlowLightIntensityPixel`s.
"""
struct SlowLightIntensityScreen{T, A <:AbstractMatrix} <: AbstractScreen
    "Minimum and Maximum Bardeen α values"
    αrange::NTuple{2, T}

    "Minimum and Maximum Bardeen β values"
    βrange::NTuple{2, T}

    "Data type that stores screen pixel information"
    pixels::A

    @kernel function _generate_screen!(screen, met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res) where T
        I,J = @index(Global, NTuple)
        α = αmin + (αmax - αmin) * (T(I)-1) / (res-1)
        β = βmin + (βmax - βmin) * (T(J)-1) / (res-1)
        screen[I, J] = SlowLightIntensityPixel(met, α, β, θo)
    end
    @doc """
        SlowLightIntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res; A=Matrix) where {T}

    Construct a `SlowLightIntensityScreen` object.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric object.
    - `αmin`: Minimum Bardeen α value.
    - `αmax`: Maximum Bardeen α value.
    - `βmin`: Minimum Bardeen β value.
    - `βmax`: Maximum Bardeen β value.
    - `θo`: Observer's inclination angle. θo ∈ (0, π).
    - `res`: Resolution of the screen (number of pixels along one dimension).
    - `A=Matrix`: Data type that stores screen pixel information (default is `Matrix`). A GPUMatrix can be used for GPU computations.

    # Returns
    A `SlowLightIntensityScreen` object.
    """
    function SlowLightIntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res; A=Matrix) where {T}
        screen = A(Matrix{SlowLightIntensityPixel{T}}(undef, res, res))

        backend = get_backend(screen)

        _generate_screen!(backend)(screen, met, αmin, αmax, βmin, βmax, θo, res, ndrange = (res, res))
        
        new{T, typeof(screen)}((αmin, αmax), (βmin, βmax), screen)
    end
end

"""
    $TYPEDEF

Camera that caches slow light raytracing information for an observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct SlowLightIntensityCamera{T, A} <: AbstractCamera
    metric::Kerr{T}
    "Data type that stores screen pixel information"
    screen::SlowLightIntensityScreen{T, A}
    "Observer screen_coordinate"
    screen_coordinate::NTuple{2, T}
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
    function SlowLightIntensityCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res; A=Matrix) where {T}
        screen = SlowLightIntensityScreen(met, αmin, αmax, βmin, βmax, θo, res; A)
        new{T,typeof(screen.pixels)}(met, screen, (T(Inf), θo))
    end
end

function η(pix::SlowLightIntensityPixel) return pix.η end
function λ(pix::SlowLightIntensityPixel) return pix.λ end
function roots(pix::SlowLightIntensityPixel) return pix.roots end
function screen_coordinate(pix::SlowLightIntensityPixel) return pix.screen_coordinate end
function inclination(pix::SlowLightIntensityPixel) return pix.θo end
function I0_inf(pix::SlowLightIntensityPixel) return pix.I0_inf end
function Ir_inf(pix::SlowLightIntensityPixel) return pix.I0_inf end
function I1_inf_m_I0_terms(pix::SlowLightIntensityPixel) return pix.I1_inf_m_I0_terms end
function I2_inf_m_I0_terms(pix::SlowLightIntensityPixel) return pix.I2_inf_m_I0_terms end
function Ip_inf_m_I0_terms(pix::SlowLightIntensityPixel) return pix.Ip_inf_m_I0_terms end
function Im_inf_m_I0_terms(pix::SlowLightIntensityPixel) return pix.Im_inf_m_I0_terms end
function radial_inf_integrals_m_I0_terms(pix::SlowLightIntensityPixel) return I1_inf_m_I0_terms(pix), I2_inf_m_I0_terms(pix), Ip_inf_m_I0_terms(pix), Im_inf_m_I0_terms(pix) end
function Iϕ_inf(pix::SlowLightIntensityPixel) return pix.Iϕ_inf end
function It_inf(pix::SlowLightIntensityPixel) return pix.It_inf end
function absGθo_Gθhat(pix::SlowLightIntensityPixel) return pix.absGθo_Gθhat end
function absGϕo_Gϕhat(pix::SlowLightIntensityPixel) return pix.absGϕo_Gϕhat end
function absGto_Gthat(pix::SlowLightIntensityPixel) return pix.absGto_Gthat end
