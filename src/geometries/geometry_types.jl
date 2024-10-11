export ⊕, add, Scene
"""
    $TYPEDEF

Abstract Geometry
"""
abstract type AbstractGeometry end

"""
    $TYPEDEF

Abstract Material
"""
abstract type AbstractMaterial end



"""
    $TYPEDEF

Cone Geometry with half opening angle `opening_angle`.
"""
struct ConeGeometry{T, A} <: AbstractGeometry
    opening_angle::T
    attributes::A
    ConeGeometry(opening_angle::T) where {T} = new{T, Nothing}(opening_angle)
    ConeGeometry(opening_angle::T, attributes::A) where {T,A} = new{T, A}(opening_angle,attributes)
end

Disk(;attributes=nothing) = ConeGeometry(π/2;attributes=attributes)

"""
    $TYPEDEF
"""
struct Mesh{G<:AbstractGeometry,M<:AbstractMaterial}
    geometry::G
    material::M
end

const Scene = NTuple{N, Mesh} where {N}
Scene() = NTuple{0, Mesh}()

function add(scene::Scene, mesh::Mesh)
    return (scene..., mesh)
end
