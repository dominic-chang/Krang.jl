export IntensityCamera
"""
    $TYPEDEF

Intensity Pixel Type. 
Each Pixel is associated with a single ray, and caches some information about the ray.
"""
struct IntensityPixel{T1, T2, T3, T4, T5, T6, T7, T8, T9} <: AbstractPixel
    metric::Kerr{T1}
    screen_coordinate::NTuple{2,T2}
    "Radial roots"
    roots::NTuple{4,T3}
    "Radial antiderivative"
    I0_inf::T4
    "Total possible Mino time"
    total_mino_time::T5
    "Angular antiderivative"
    absGθo_Gθhat::NTuple{2,T6}
    "Inclination"
    θo::T7
    η::T8
    λ::T9
end

@doc """
    IntensityPixel(met::Kerr{T}, α, β, θo) where {T}

Construct an `IntensityPixel` object with the given Kerr metric, screen coordinates, and inclination.

# Arguments
- `met::Kerr{T}`: The Kerr metric.
- `α`: The Bardeen α value (screen coordinate).
- `β`: The Bardeen β value (screen coordinate).
- `θo`: The inclination angle.

# Returns
- An `IntensityPixel` object initialized with the given parameters.

# Details
This function calculates the η and λ values using the provided Kerr metric and screen coordinates. 
It then computes the radial roots and adjusts them if necessary. 
Finally, it initializes an `IntensityPixel` object with the calculated values and the provided parameters.
"""
function IntensityPixel(met::Kerr{T}, α, β, θo) where {T}
    tempη = Krang.η(met, α, β, θo)
    tempλ = Krang.λ(met, α, θo)
    roots = Krang.get_radial_roots(met, tempη, tempλ)
    numreals = sum(_isreal2, roots)
    if (numreals == 2) && (abs(imag(roots[4]) / real(roots[4])) < eps(T))
        roots = (roots[1], roots[4], roots[2], roots[3])
    end
    I0_inf = Krang.Ir_inf(met, roots)
    τ_total = total_mino_time(met, roots)
    Gθo_Gθhat = Krang._absGθo_Gθhat(met, θo, tempη, tempλ)
    IntensityPixel(
        met,
        (α, β),
        roots,
        I0_inf,
        τ_total,
        Gθo_Gθhat,
        θo,
        tempη,
        tempλ,
    )
end

"""
    $TYPEDEF

Screen made of `IntensityPixel`s.
"""
struct IntensityScreen{A<:AbstractMatrix} <: AbstractScreen
    "Minimum and Maximum Bardeen α values"
    αrange::NTuple{2}

    "Minimum and Maximum Bardeen β values"
    βrange::NTuple{2}

    "Data type that stores screen pixel information"
    pixels::A

    IntensityScreen{A}(αrange::NTuple{2}, βrange::NTuple{2}, pixels::A) where {A<:AbstractMatrix} =
        new{A}(αrange, βrange, pixels)

    function IntensityScreen(met::Kerr, αmin, αmax, βmin, βmax, θo, res)
        screen = Matrix{IntensityPixel}(undef, res, res)
        αvals = range(αmin, αmax, length=res)
        βvals = range(βmin, βmax, length=res)
        for (iα, α) in enumerate(αvals)
            for (iβ, β) in enumerate(βvals)
                screen[iα, iβ] = IntensityPixel(met, α, β, θo)
            end
        end
        new{typeof(screen)}((αmin, αmax), (βmin, βmax), screen)
    end
end

"""
    $TYPEDEF

Camera that caches fast light raytracing information for an observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct IntensityCamera{A} <: AbstractCamera
    metric::Kerr
    "Data type that stores screen pixel information"
    screen::IntensityScreen{A}
    "Observer screen_coordinate"
    screen_coordinate::NTuple{2}

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
        new{typeof(screen.pixels)}(
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
