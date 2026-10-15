export IntensityCamera
"""
    $TYPEDEF

Intensity Pixel Type. 
Each Pixel is associated with a single ray, and caches some information about the ray.
"""
struct IntensityPixel{T1,T2,T3,T4,T5,T6,T7,T8,T9} <: AbstractPixel
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
It then computes the radial roots and the radial and angular antiderivatives.
Finally, it initializes an `IntensityPixel` object with the calculated values and the provided parameters.
"""
function IntensityPixel(met::Kerr{T}, α, β, θo) where {T}
    tempη = Krang.η(met, α, β, θo)
    tempλ = Krang.λ(met, α, θo)
    roots = Krang.get_radial_roots(met, tempη, tempλ)
    I0_inf = Krang.Ir_inf(met, roots)
    τ_total = total_mino_time(met, roots)
    Gθo_Gθhat = Krang._absGθo_Gθhat(met, θo, tempη, tempλ)
    IntensityPixel(met, (α, β), roots, I0_inf, τ_total, Gθo_Gθhat, θo, tempη, tempλ)
end

"""
    $TYPEDEF

Screen of fast-light pixels, stored as a struct of arrays: one array per cached geodesic
invariant rather than one `IntensityPixel` per screen position.

`metric` and the inclination `θo` are constant across the screen, so they are held once here
instead of being duplicated into every pixel. [`IntensityPixel`](@ref) is lazily constructed 
when `screen.pixels` is indexed. This layout is what allows the cached invariants to cross a 
`Reactant.@compile` boundary as data.
"""
struct IntensityScreen{T,M<:AbstractMatrix{T},MC<:AbstractMatrix{Complex{T}}} <:
       AbstractScreen
    "Minimum and Maximum Bardeen α values"
    αrange::NTuple{2,T}
    "Minimum and Maximum Bardeen β values"
    βrange::NTuple{2,T}
    "Spacetime metric, shared by every pixel"
    metric::Kerr{T}
    "Observer inclination, shared by every pixel"
    θo::T
    α::M
    β::M
    r1::MC
    r2::MC
    r3::MC
    r4::MC
    I0_inf::M
    total_mino_time::M
    Gθo::M
    Gθhat::M
    η::M
    λ::M
end

"Concrete `IntensityPixel` type produced by a screen whose element type is `T`."
@inline _intensity_pixel_type(::Type{T}) where {T} =
    IntensityPixel{T,T,Complex{T},T,T,T,T,T,T}

"""
    $TYPEDEF

Lazy `AbstractMatrix{<:IntensityPixel}` view over an [`IntensityScreen`](@ref). `getindex`
assembles an `IntensityPixel` from the screen's field arrays, so existing code that iterates,
indexes or broadcasts over `screen.pixels` keeps working unchanged.
"""
struct IntensityPixelArray{P,S<:IntensityScreen} <: AbstractMatrix{P}
    screen::S
end

@inline Base.size(pa::IntensityPixelArray) = size(pa.screen.α)
Base.IndexStyle(::Type{<:IntensityPixelArray}) = IndexCartesian()

@inline function Base.getindex(pa::IntensityPixelArray, i::Int, j::Int)
    s = pa.screen
    @inbounds IntensityPixel(
        s.metric,
        (s.α[i, j], s.β[i, j]),
        (s.r1[i, j], s.r2[i, j], s.r3[i, j], s.r4[i, j]),
        s.I0_inf[i, j],
        s.total_mino_time[i, j],
        (s.Gθo[i, j], s.Gθhat[i, j]),
        s.θo,
        s.η[i, j],
        s.λ[i, j],
    )
end

# `similar` must hand back real storage (the view is read-only), and must stay on whatever
# device the screen lives on.
@inline Base.similar(pa::IntensityPixelArray, ::Type{T}, dims::Dims) where {T} =
    similar(pa.screen.α, T, dims)

# Keeps `camera.screen.pixels` working across src/, ext/, the tests and the examples.
@inline function Base.getproperty(s::IntensityScreen{T}, name::Symbol) where {T}
    name === :pixels && return IntensityPixelArray{_intensity_pixel_type(T),typeof(s)}(s)
    return getfield(s, name)
end

"""
    IntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res) where {T}

Build a fast-light screen by evaluating [`IntensityPixel`](@ref) on a `res`×`res` grid of
Bardeen coordinates and scattering each cached invariant into its own array.
"""
function IntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res) where {T}
    θoT = T(θo)
    αvals = range(T(αmin), T(αmax), length = res)
    βvals = range(T(βmin), T(βmax), length = res)

    α = Matrix{T}(undef, res, res)
    β = similar(α)
    I0_inf = similar(α)
    total_mino_time = similar(α)
    Gθo = similar(α)
    Gθhat = similar(α)
    ηm = similar(α)
    λm = similar(α)
    r1 = Matrix{Complex{T}}(undef, res, res)
    r2 = similar(r1)
    r3 = similar(r1)
    r4 = similar(r1)

    for j = 1:res, i = 1:res
        pix = IntensityPixel(met, αvals[i], βvals[j], θoT)
        α[i, j], β[i, j] = pix.screen_coordinate
        r1[i, j], r2[i, j], r3[i, j], r4[i, j] = pix.roots
        I0_inf[i, j] = pix.I0_inf
        total_mino_time[i, j] = pix.total_mino_time
        Gθo[i, j], Gθhat[i, j] = pix.absGθo_Gθhat
        ηm[i, j] = pix.η
        λm[i, j] = pix.λ
    end

    return IntensityScreen(
        (T(αmin), T(αmax)),
        (T(βmin), T(βmax)),
        met,
        θoT,
        α,
        β,
        r1,
        r2,
        r3,
        r4,
        I0_inf,
        total_mino_time,
        Gθo,
        Gθhat,
        ηm,
        λm,
    )
end

"""
    $TYPEDEF

Camera that caches fast light raytracing information for an observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct IntensityCamera{S<:IntensityScreen} <: AbstractCamera
    metric::Kerr
    "Data type that stores screen pixel information"
    screen::S
    "Observer screen_coordinate"
    screen_coordinate::NTuple{2}

    @doc """
        IntensityCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res::Int) where {T}

    Constructor for an Intensity Camera.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric.
    - `θo`: Observer's inclination angle. θo ∈ (0, π).
    - `αmin`: Minimum value of α.
    - `αmax`: Maximum value of α.
    - `βmin`: Minimum value of β.
    - `βmax`: Maximum value of β.
    - `res::Int`: Resolution of the screen.
    # Returns
    - `IntensityCamera{S}`: An intensity camera object, where `S` is the screen type.
    """
    function IntensityCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res::Int) where {T}
        screen = IntensityScreen(met, αmin, αmax, βmin, βmax, θo, res)
        new{typeof(screen)}(met, screen, (T(Inf), θo))
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

# `Adapt` support: lets the struct-of-arrays screen (and the lazy pixel view over it) be moved
# onto a device and, critically, be converted to device-array types there -- an array of
# pixels needed none of this, a struct holding arrays does.
function Adapt.adapt_structure(to, s::IntensityScreen)
    return IntensityScreen(
        s.αrange,
        s.βrange,
        s.metric,
        s.θo,
        adapt(to, getfield(s, :α)),
        adapt(to, getfield(s, :β)),
        adapt(to, getfield(s, :r1)),
        adapt(to, getfield(s, :r2)),
        adapt(to, getfield(s, :r3)),
        adapt(to, getfield(s, :r4)),
        adapt(to, getfield(s, :I0_inf)),
        adapt(to, getfield(s, :total_mino_time)),
        adapt(to, getfield(s, :Gθo)),
        adapt(to, getfield(s, :Gθhat)),
        adapt(to, getfield(s, :η)),
        adapt(to, getfield(s, :λ)),
    )
end

function Adapt.adapt_structure(to, pa::IntensityPixelArray{P}) where {P}
    screen = adapt(to, pa.screen)
    return IntensityPixelArray{P,typeof(screen)}(screen)
end

# Collecting the view would otherwise scalar-index the backing arrays, which GPU array types
# forbid. Move the fields to the host once, then materialise the pixels.
function Base.Array(pa::IntensityPixelArray{P}) where {P}
    host = adapt(Array, pa)
    return P[host[i, j] for i in axes(host, 1), j in axes(host, 2)]
end
