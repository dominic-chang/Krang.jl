export IntensityCamera
"""
    $TYPEDEF

Intensity Pixel Type. 
Each Pixel is associated with a single ray, and caches some information about the ray.
"""
struct IntensityPixel{T} <: AbstractPixel{T}
    metric::Kerr{T}
    screen_coordinate::NTuple{2,T}
    "Radial roots"
    roots::NTuple{4,Complex{T}}
    "Radial antiderivative"
    I0_inf::T
    "Total possible Mino time"
    total_mino_time::T
    "Angular antiderivative"
    absGθo_Gθhat::NTuple{2,T}
    "Inclination"
    θo::T
    η::T
    λ::T
    @doc """
        IntensityPixel(met::Kerr{T}, α::T, β::T, θo::T) where {T}

    Construct an `IntensityPixel` object with the given Kerr metric, screen coordinates, and inclination.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric.
    - `α::T`: The Bardeen α value (screen coordinate).
    - `β::T`: The Bardeen β value (screen coordinate).
    - `θo::T`: The inclination angle.

    # Returns
    - An `IntensityPixel` object initialized with the given parameters.

    # Details
    This function calculates the η and λ values using the provided Kerr metric and screen coordinates. 
    It then computes the radial roots and adjusts them if necessary. 
    Finally, it initializes an `IntensityPixel` object with the calculated values and the provided parameters.
    """
    function IntensityPixel(met::Kerr{T}, α::T, β::T, θo::T) where {T}
        tempη = Krang.η(met, α, β, θo)
        tempλ = Krang.λ(met, α, θo)
        roots = Krang.get_radial_roots(met, tempη, tempλ)
        numreals = sum(_isreal2.(roots))
        if (numreals == 2) && (abs(imag(roots[4]) / real(roots[4])) < eps(T))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        I0_inf = Krang.Ir_inf(met, roots)
        new{T}(
            met,
            (α, β),
            roots,
            I0_inf,
            total_mino_time(met, roots),
            Krang._absGθo_Gθhat(met, θo, tempη, tempλ),
            θo,
            tempη,
            tempλ,
        )
    end
end

"""
    $TYPEDEF

Screen made of `IntensityPixel`s.
"""
struct IntensityScreen{T,A<:AbstractMatrix} <: AbstractScreen
    "Minimum and Maximum Bardeen α values"
    αrange::NTuple{2,T}

    "Minimum and Maximum Bardeen β values"
    βrange::NTuple{2,T}

    "Data type that stores screen pixel information"
    pixels::A

    function IntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res) where {T}
        screen = Matrix{IntensityPixel}(undef, res, res)
        αvals = range(αmin, αmax, length=res)
        βvals = range(βmin, βmax, length=res)
        for (iα, α) in enumerate(αvals)
            for (iβ, β) in enumerate(βvals)
                screen[iα, iβ] = IntensityPixel(met, α, β, θo)
            end
        end
        new{T, typeof(screen)}((αmin, αmax), (βmin, βmax), screen)
    end
end

"""
    $TYPEDEF

Camera that caches fast light raytracing information for an observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct IntensityCamera{T,A} <: AbstractCamera
    metric::Kerr{T}
    "Data type that stores screen pixel information"
    screen::IntensityScreen{T,A}
    "Observer screen_coordinate"
    screen_coordinate::NTuple{2,T}

    @doc """
        IntensityCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res::Int; A=Matrix) where {T}

    Constructor for an Intensity Camera.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric.
    - `θo`: Observer's inclination angle. θo ∈ (0, π).
    - `αmin`: Minimum value of α.
    - `αmax`: Maximum value of α.
    - `βmin`: Minimum value of β.
    - `βmax`: Maximum value of β.
    - `res::Int`: Resolution of the screen.
    - `A=Matrix`: Optional argument to specify the type of matrix to use. A GPUMatrix can be used for GPU computations.

    # Returns
    - `IntensityCamera{T, A}`: An intensity camera object.
    """
    function IntensityCamera(
        met::Kerr{T},
        θo,
        αmin,
        αmax,
        βmin,
        βmax,
        res::Int
    ) where {T}
        screen = IntensityScreen(met, αmin, αmax, βmin, βmax, θo, res)
        new{T,typeof(screen.pixels)}(
            met,
            screen,
            (T(Inf), θo),
        )
    end
end

function η(pix::IntensityPixel)
    return pix.η
end
function λ(pix::IntensityPixel)
    return pix.λ
end
function roots(pix::IntensityPixel)
    return pix.roots
end
function screen_coordinate(pix::IntensityPixel)
    return pix.screen_coordinate
end
function inclination(pix::IntensityPixel)
    return pix.θo
end
function I0_inf(pix::IntensityPixel)
    return pix.I0_inf
end
function total_mino_time(pix::IntensityPixel)
    return pix.total_mino_time
end
function Ir_inf(pix::IntensityPixel)
    return pix.I0_inf
end
function absGθo_Gθhat(pix::IntensityPixel)
    return pix.absGθo_Gθhat
end
