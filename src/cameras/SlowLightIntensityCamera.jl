export SlowLightIntensityCamera

"""
    $TYPEDEF

Slow Light Intensity Pixel Type.
Each Pixel is associated with a single ray, and caches some information about the ray.
On top of what an [`IntensityPixel`](@ref) caches, a slow light pixel also caches the
antiderivatives needed to recover the azimuthal angle and the coordinate time of emission.
"""
struct SlowLightIntensityPixel{
    T1,
    T2,
    T3,
    T4,
    T5,
    T6,
    T7,
    T8,
    T9,
    T10,
    T11,
    T12,
    T13,
    T14,
    T15,
    T16,
    T17,
} <: AbstractPixel
    metric::Kerr{T1}
    screen_coordinate::NTuple{2,T2}
    "Radial roots"
    roots::NTuple{4,T3}
    "Radial antiderivative"
    I0_inf::T4
    "Total possible Mino time"
    total_mino_time::T5
    "Radial ϕ antiderivative"
    Iϕ_inf::T6
    "Radial t antiderivative"
    It_inf::T7
    "Linear radial antiderivative, with its I0 divergence subtracted"
    I1_inf_m_I0_terms::T8
    "Quadratic radial antiderivative, with its I0 divergence subtracted"
    I2_inf_m_I0_terms::T9
    "Outer horizon radial antiderivative, with its I0 divergence subtracted"
    Ip_inf_m_I0_terms::T10
    "Inner horizon radial antiderivative, with its I0 divergence subtracted"
    Im_inf_m_I0_terms::T11
    "Angular antiderivative"
    absGθo_Gθhat::NTuple{2,T12}
    "Angular ϕ antiderivative"
    absGϕo_Gϕhat::NTuple{2,T13}
    "Angular t antiderivative"
    absGto_Gthat::NTuple{2,T14}
    "Inclination"
    θo::T15
    η::T16
    λ::T17
end

@doc """
    SlowLightIntensityPixel(met::Kerr{T}, α, β, θo) where {T}

Construct a `SlowLightIntensityPixel` object with the given Kerr metric, screen coordinates, and inclination.

# Arguments
- `met::Kerr{T}`: The Kerr metric.
- `α`: The Bardeen α value (screen coordinate).
- `β`: The Bardeen β value (screen coordinate).
- `θo`: The inclination angle.

# Returns
- A `SlowLightIntensityPixel` object initialized with the given parameters.

# Details
This function calculates the η and λ values using the provided Kerr metric and screen coordinates.
It then computes the radial roots, the radial and angular antiderivatives, and the extra ϕ and t
antiderivatives that fast light does not need.
Finally, it initializes a `SlowLightIntensityPixel` object with the calculated values and the provided parameters.
"""
function SlowLightIntensityPixel(met::Kerr{T}, α, β, θo) where {T}
    tempη = Krang.η(met, α, β, θo)
    tempλ = Krang.λ(met, α, θo)
    roots = Krang.get_radial_roots(met, tempη, tempλ)
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

Screen of slow-light pixels, stored as a struct of arrays: one array per cached geodesic
invariant rather than one `SlowLightIntensityPixel` per screen position.

`metric` and the inclination `θo` are constant across the screen, so they are held once here
instead of being duplicated into every pixel. [`SlowLightIntensityPixel`](@ref) is lazily
constructed when `screen.pixels` is indexed. This layout is what allows the cached invariants
to cross a `Reactant.@compile` boundary as data.
"""
struct SlowLightIntensityScreen{T,M<:AbstractMatrix{T},MC<:AbstractMatrix{Complex{T}}} <:
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
    Iϕ_inf::M
    It_inf::M
    I1_inf_m_I0_terms::M
    I2_inf_m_I0_terms::M
    Ip_inf_m_I0_terms::M
    Im_inf_m_I0_terms::M
    Gθo::M
    Gθhat::M
    Gϕo::M
    Gϕhat::M
    Gto::M
    Gthat::M
    η::M
    λ::M
end

"Concrete `SlowLightIntensityPixel` type produced by a screen whose element type is `T`."
@inline _slow_light_intensity_pixel_type(::Type{T}) where {T} =
    SlowLightIntensityPixel{T,T,Complex{T},T,T,T,T,T,T,T,T,T,T,T,T,T,T}

"""
    $TYPEDEF

Lazy `AbstractMatrix{<:SlowLightIntensityPixel}` view over a
[`SlowLightIntensityScreen`](@ref). `getindex` assembles a `SlowLightIntensityPixel` from the
screen's field arrays, so existing code that iterates, indexes or broadcasts over
`screen.pixels` keeps working unchanged.
"""
struct SlowLightIntensityPixelArray{P,S<:SlowLightIntensityScreen} <: AbstractMatrix{P}
    screen::S
end

@inline Base.size(pa::SlowLightIntensityPixelArray) = size(pa.screen.α)
Base.IndexStyle(::Type{<:SlowLightIntensityPixelArray}) = IndexCartesian()

@inline function Base.getindex(pa::SlowLightIntensityPixelArray, i::Int, j::Int)
    s = pa.screen
    @inbounds SlowLightIntensityPixel(
        s.metric,
        (s.α[i, j], s.β[i, j]),
        (s.r1[i, j], s.r2[i, j], s.r3[i, j], s.r4[i, j]),
        s.I0_inf[i, j],
        s.total_mino_time[i, j],
        s.Iϕ_inf[i, j],
        s.It_inf[i, j],
        s.I1_inf_m_I0_terms[i, j],
        s.I2_inf_m_I0_terms[i, j],
        s.Ip_inf_m_I0_terms[i, j],
        s.Im_inf_m_I0_terms[i, j],
        (s.Gθo[i, j], s.Gθhat[i, j]),
        (s.Gϕo[i, j], s.Gϕhat[i, j]),
        (s.Gto[i, j], s.Gthat[i, j]),
        s.θo,
        s.η[i, j],
        s.λ[i, j],
    )
end

# `similar` must hand back real storage (the view is read-only), and must stay on whatever
# device the screen lives on.
@inline Base.similar(pa::SlowLightIntensityPixelArray, ::Type{T}, dims::Dims) where {T} =
    similar(pa.screen.α, T, dims)

# Keeps `camera.screen.pixels` working across src/, ext/, the tests and the examples.
@inline function Base.getproperty(s::SlowLightIntensityScreen{T}, name::Symbol) where {T}
    name === :pixels &&
        return SlowLightIntensityPixelArray{_slow_light_intensity_pixel_type(T),typeof(s)}(
            s,
        )
    return getfield(s, name)
end

"""
    SlowLightIntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res) where {T}

Build a slow-light screen by evaluating [`SlowLightIntensityPixel`](@ref) on a `res`×`res` grid
of Bardeen coordinates and scattering each cached invariant into its own array.
"""
function SlowLightIntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res) where {T}
    θoT = T(θo)
    αvals = range(T(αmin), T(αmax), length = res)
    βvals = range(T(βmin), T(βmax), length = res)

    α = Matrix{T}(undef, res, res)
    β = similar(α)
    I0_inf = similar(α)
    total_mino_time = similar(α)
    Iϕ_inf = similar(α)
    It_inf = similar(α)
    I1_inf_m_I0_terms = similar(α)
    I2_inf_m_I0_terms = similar(α)
    Ip_inf_m_I0_terms = similar(α)
    Im_inf_m_I0_terms = similar(α)
    Gθo = similar(α)
    Gθhat = similar(α)
    Gϕo = similar(α)
    Gϕhat = similar(α)
    Gto = similar(α)
    Gthat = similar(α)
    η = similar(α)
    λ = similar(α)
    r1 = Matrix{Complex{T}}(undef, res, res)
    r2 = similar(r1)
    r3 = similar(r1)
    r4 = similar(r1)

    for j = 1:res, i = 1:res
        pix = SlowLightIntensityPixel(met, αvals[i], βvals[j], θoT)
        α[i, j], β[i, j] = pix.screen_coordinate
        r1[i, j], r2[i, j], r3[i, j], r4[i, j] = pix.roots
        I0_inf[i, j] = pix.I0_inf
        total_mino_time[i, j] = pix.total_mino_time
        Iϕ_inf[i, j] = pix.Iϕ_inf
        It_inf[i, j] = pix.It_inf
        I1_inf_m_I0_terms[i, j] = pix.I1_inf_m_I0_terms
        I2_inf_m_I0_terms[i, j] = pix.I2_inf_m_I0_terms
        Ip_inf_m_I0_terms[i, j] = pix.Ip_inf_m_I0_terms
        Im_inf_m_I0_terms[i, j] = pix.Im_inf_m_I0_terms
        Gθo[i, j], Gθhat[i, j] = pix.absGθo_Gθhat
        Gϕo[i, j], Gϕhat[i, j] = pix.absGϕo_Gϕhat
        Gto[i, j], Gthat[i, j] = pix.absGto_Gthat
        η[i, j] = pix.η
        λ[i, j] = pix.λ
    end

    return SlowLightIntensityScreen(
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
        Iϕ_inf,
        It_inf,
        I1_inf_m_I0_terms,
        I2_inf_m_I0_terms,
        Ip_inf_m_I0_terms,
        Im_inf_m_I0_terms,
        Gθo,
        Gθhat,
        Gϕo,
        Gϕhat,
        Gto,
        Gthat,
        η,
        λ,
    )
end

"""
    $TYPEDEF

Camera that caches slow light raytracing information for an observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct SlowLightIntensityCamera{S<:SlowLightIntensityScreen} <: AbstractCamera
    metric::Kerr
    "Data type that stores screen pixel information"
    screen::S
    "Observer screen_coordinate"
    screen_coordinate::NTuple{2}
    @doc """
        SlowLightIntensityCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res) where {T}

    Constructs a `SlowLightIntensityCamera` object.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric.
    - `θo`: The Observer's inclination angle. θo ∈ (0, π).
    - `αmin`: Minimum α coordinate on the screen.
    - `αmax`: Maximum α coordinate on the screen.
    - `βmin`: Minimum β coordinate on the screen.
    - `βmax`: Maximum β coordinate on the screen.
    - `res`: Resolution of the screen (number of pixels along one dimension).
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
        res,
    ) where {T}
        screen = SlowLightIntensityScreen(met, αmin, αmax, βmin, βmax, θo, res)
        new{typeof(screen)}(met, screen, (T(Inf), θo))
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

# `Adapt` support: lets the struct-of-arrays screen (and the lazy pixel view over it) be moved
# onto a device and, critically, be converted to device-array types inside a KernelAbstractions
# kernel -- an array of pixels needed none of this, a struct holding arrays does.
function Adapt.adapt_structure(to, s::SlowLightIntensityScreen)
    return SlowLightIntensityScreen(
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
        adapt(to, getfield(s, :Iϕ_inf)),
        adapt(to, getfield(s, :It_inf)),
        adapt(to, getfield(s, :I1_inf_m_I0_terms)),
        adapt(to, getfield(s, :I2_inf_m_I0_terms)),
        adapt(to, getfield(s, :Ip_inf_m_I0_terms)),
        adapt(to, getfield(s, :Im_inf_m_I0_terms)),
        adapt(to, getfield(s, :Gθo)),
        adapt(to, getfield(s, :Gθhat)),
        adapt(to, getfield(s, :Gϕo)),
        adapt(to, getfield(s, :Gϕhat)),
        adapt(to, getfield(s, :Gto)),
        adapt(to, getfield(s, :Gthat)),
        adapt(to, getfield(s, :η)),
        adapt(to, getfield(s, :λ)),
    )
end

function Adapt.adapt_structure(to, pa::SlowLightIntensityPixelArray{P}) where {P}
    screen = adapt(to, pa.screen)
    return SlowLightIntensityPixelArray{P,typeof(screen)}(screen)
end

# Collecting the view would otherwise scalar-index the backing arrays, which GPU array types
# forbid. Move the fields to the host once, then materialise the pixels.
function Base.Array(pa::SlowLightIntensityPixelArray{P}) where {P}
    host = adapt(Array, pa)
    return P[host[i, j] for i in axes(host, 1), j in axes(host, 2)]
end
