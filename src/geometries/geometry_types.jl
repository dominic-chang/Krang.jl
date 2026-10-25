export ⊕, add, Scene
"""
    $TYPEDEF

Abstract Geometry
"""
abstract type AbstractGeometry end

struct Intersection{T}
    ts::T
    rs::T
    θs::T
    ϕs::T
    νr::Bool
    νθ::Bool
end

"""
    $TYPEDEF

The angular solve of a ray that is shared by all sub-images at a given `θs`.

# Fields

$TYPEDFIELDS
"""
struct AngularCache{T}
    "Emission inclination this cache was built for"
    θs::T
    "The `θs` angular antiderivative"
    Gθs::T
    "`|Gθo|`, carrying the sign the vortical branch gives it"
    Gθo::T
    "Half period `Ĝθ` between angular turning points"
    Gθhat::T
    "Whether the ray is vortical (`η < 0`)"
    isvortical::Bool
    "`|cos θs| < |cos θo|`, the test [`Gθ`](@ref) branches on"
    isincone::Bool
    "Whether the screen boundary rejects this ray; only computed for a `θs` outside the observer's cone, the one case a caller reads it"
    βreject::Bool
    "Whether the ray reaches `θs` at all"
    issuccess::Bool
end

"""
    AngularCache(pix::AbstractPixel, θs)

Build the [`Krang.AngularCache`](@ref) for the ray at `pix` and the cone of opening angle
`θs`.
"""
@inline function AngularCache(pix::AbstractPixel, θs::T) where {T}
    met = metric(pix)
    θo = inclination(pix)
    α, β = screen_coordinate(pix)

    cosθo = cos(θo)
    # The emission functions only consult the screen boundary for a cone that lies outside
    # the observer's own, so for `isincone` -- an equatorial disk, say -- it is never read
    # and never worth computing.
    βreject = false
    if !(θo ≤ θs ≤ (π - θo) || (π - θo) ≤ θs ≤ θo)
        αmin = @inline αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? @inline(βboundary(met, α, θo, θs)) : zero(T))
        βreject = (abs(β) + eps(T)) < βbound
    end

    # `cos(θs) > |cos(θo)|` implies the `isincone` test of every emission function fails,
    # so a rejected ray here is rejected for every sub-image and the angular solve below
    # would be thrown away.
    if βreject && (cos(θs) > abs(cosθo))
        z = zero(T)
        return AngularCache(θs, z, z, z, η(pix) < z, abs(cos(θs)) < abs(cosθo), true, false)
    end

    Go, Ghat = absGθo_Gθhat(pix)
    Gθs, Gθo, Gθhat, isvortical, isincone, issuccess =
        @inline Krang._Gθ_setup(met, θo, η(pix), λ(pix), Go, Ghat, θs)
    return AngularCache(θs, Gθs, Gθo, Gθhat, isvortical, isincone, βreject, issuccess)
end

"""
    $TYPEDEF

Pairs an `AbstractGeometry` with the `AbstractMaterial` that emits on it. A `Scene` is a
tuple of these, and rendering sums each mesh's contribution to every pixel.
"""
struct Mesh{G<:AbstractGeometry,M<:AbstractMaterial}
    geometry::G
    material::M
end

const Scene = NTuple{N,Mesh} where {N}
Scene() = NTuple{0,Mesh}()

function add(scene::Scene, mesh::Mesh)
    return (scene..., mesh)
end

Rotations.:⊕(scene::Scene, mesh::Mesh) = add(scene, mesh)

"""
    $TYPEDEF

Cone Geometry with half opening angle `opening_angle`.
"""
struct ConeGeometry{T,A} <: AbstractGeometry
    opening_angle::T
    attributes::A
    ConeGeometry(opening_angle::T) where {T} = new{T,Nothing}(opening_angle)
    ConeGeometry(opening_angle::T, attributes::A) where {T,A} =
        new{T,A}(opening_angle, attributes)
end

"""
    _intersection(::AbstractLightTrait, ::AbstractSymmetryTrait, pix, θs, isindir, n)
    _intersection(::AbstractLightTrait, ::AbstractSymmetryTrait, pix, cache::AngularCache, isindir, n)

Emission coordinates on the `θs` cone, dispatching on the material's
[`AbstractLightTrait`](@ref) and [`AbstractSymmetryTrait`](@ref). Fast-light materials
skip the emission time and axisymmetric ones also skip the emission azimuth; skipped
coordinates are reported as zero. The second form takes the sub-image independent part of
the angular solve from a [`Krang.AngularCache`](@ref) rather than recomputing it.
"""
@inline function _intersection(::FastLightTrait, ::AxisymmetricTrait, pix, θs, isindir, n)
    rs, νr, νθ, _, issuccess = @inline emission_radius(pix, θs, isindir, n)
    return Intersection(zero(rs), rs, θs, zero(rs), νr, νθ), issuccess
end

@inline function _intersection(
    ::FastLightTrait,
    ::AxisymmetricTrait,
    pix,
    cache::AngularCache,
    isindir,
    n,
)
    rs, νr, νθ, _, issuccess = @inline emission_radius(pix, cache, isindir, n)
    return Intersection(zero(rs), rs, cache.θs, zero(rs), νr, νθ), issuccess
end

@inline function _intersection(
    ::FastLightTrait,
    ::NonAxisymmetricTrait,
    pix,
    θs,
    isindir,
    n,
)
    rs, ϕs, νr, νθ, issuccess = @inline emission_coordinates_fast_light(pix, θs, isindir, n)
    return Intersection(zero(rs), rs, θs, ϕs, νr, νθ), issuccess
end

@inline function _intersection(
    ::FastLightTrait,
    ::NonAxisymmetricTrait,
    pix,
    cache::AngularCache,
    isindir,
    n,
)
    rs, ϕs, νr, νθ, issuccess =
        @inline emission_coordinates_fast_light(pix, cache, isindir, n)
    return Intersection(zero(rs), rs, cache.θs, ϕs, νr, νθ), issuccess
end

@inline function _intersection(
    ::SlowLightTrait,
    ::AbstractSymmetryTrait,
    pix,
    θs,
    isindir,
    n,
)
    ts, rs, ϕs, νr, νθ, issuccess = @inline emission_coordinates(pix, θs, isindir, n)
    return Intersection(ts, rs, θs, ϕs, νr, νθ), issuccess
end

@inline function _intersection(
    ::SlowLightTrait,
    ::AbstractSymmetryTrait,
    pix,
    cache::AngularCache,
    isindir,
    n,
)
    ts, rs, ϕs, νr, νθ, issuccess = @inline emission_coordinates(pix, cache, isindir, n)
    return Intersection(ts, rs, cache.θs, ϕs, νr, νθ), issuccess
end

@inline function _raytrace(
    observation,
    pix::AbstractPixel,
    mesh::Mesh{<:ConeGeometry{T,A},<:AbstractMaterial};
    res,
) where {T,A}
    geometry = mesh.geometry
    material = mesh.material
    θs = geometry.opening_angle
    subimgs = material.subimgs

    # The angular solve is the same for every sub-image and for both values of `isindir`,
    # so it is done once here rather than inside the loop.
    cache = @inline AngularCache(pix, θs)

    isindir = false
    for _ = 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        for n in subimgs
            #νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
            intersection, issuccess = @inline _intersection(
                lightTrait(material),
                symmetryTrait(material),
                pix,
                cache,
                isindir,
                n,
            )

            if issuccess && (horizon(metric(pix)) < intersection.rs < T(Inf))
                observation = _accumulate(observation, @inline(material(pix, intersection)))
            end
        end
    end

    return observation
end


Disk(; attributes = nothing) = ConeGeometry(π / 2, attributes)
