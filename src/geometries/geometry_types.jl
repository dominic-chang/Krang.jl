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

Emission coordinates on the `θs` cone, dispatching on the material's
[`AbstractLightTrait`](@ref) and [`AbstractSymmetryTrait`](@ref). Fast-light materials
skip the emission time and axisymmetric ones also skip the emission azimuth; skipped
coordinates are reported as zero.
"""
@inline function _intersection(::FastLightTrait, ::AxisymmetricTrait, pix, θs, isindir, n)
    rs, νr, νθ, _, issuccess = @inline emission_radius(pix, θs, isindir, n)
    return Intersection(zero(rs), rs, θs, zero(rs), νr, νθ), issuccess
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

    isindir = false
    for _ = 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        for n in subimgs
            #νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
            intersection, issuccess = @inline _intersection(
                lightTrait(material),
                symmetryTrait(material),
                pix,
                θs,
                isindir,
                n,
            )

            if issuccess && (horizon(metric(pix)) < intersection.rs < T(Inf))
                observation += @inline(material(pix, intersection))
            end
        end
    end

    return observation
end


Disk(; attributes = nothing) = ConeGeometry(π / 2, attributes)
