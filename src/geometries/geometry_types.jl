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

    isindir = true
    for n in subimgs
        #νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
        intersection, issuccess = if isFastLight(material)
            if isAxisymmetric(material)
                rs, νr, νθ, _, issuccess = @inline emission_radius(pix, θs, isindir, n)
                Intersection(zero(rs), rs, θs, zero(rs), νr, νθ), issuccess
            else
                rs, ϕs, νr, νθ, issuccess =
                    @inline emission_coordinates_fast_light(pix, θs, isindir, n)
                Intersection(zero(rs), rs, θs, ϕs, νr, νθ), issuccess
            end
        else
            ts, rs, ϕs, νr, νθ, issuccess = @inline emission_coordinates(pix, θs, isindir, n)
            Intersection(ts, rs, θs, ϕs, νr, νθ), issuccess
        end
        
        if issuccess && (horizon(metric(pix)) < rs < T(Inf))
            observation += @inline(material(pix, intersection))
        end
    end

    isindir = false
    for n in subimgs
        #νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
        intersection, issuccess = if isFastLight(material)
            if isAxisymmetric(material)
                rs, νr, νθ, _, issuccess = @inline emission_radius(pix, θs, isindir, n)
                Intersection(zero(rs), rs, θs, zero(rs), νr, νθ), issuccess
            else
                rs, ϕs, νr, νθ, issuccess =
                    @inline emission_coordinates_fast_light(pix, θs, isindir, n)
                Intersection(zero(rs), rs, θs, ϕs, νr, νθ), issuccess
            end
        else
            ts, rs, ϕs, νr, νθ, issuccess = @inline emission_coordinates(pix, θs, isindir, n)
            Intersection(ts, rs, θs, ϕs, νr, νθ), issuccess
        end

        if issuccess && (horizon(metric(pix)) < rs < T(Inf))
            observation += @inline(material(pix, intersection))
        end
    end

    return observation
end


Disk(; attributes = nothing) = ConeGeometry(π / 2; attributes = attributes)
