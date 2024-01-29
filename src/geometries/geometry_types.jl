export ⊕
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
    attriributes::A
    ConeGeometry(opening_angle::T) where {T} = new{T, Nothing}(opening_angle)
    ConeGeometry(opening_angle::T, attriributes::A) where {T,A} = new{T, A}(opening_angle,attriributes)
end

Disk(;attributes=nothing) = ConeGeometry(π/2;attributes=attributes)

"""
    $TYPEDEF

Geometry that is comprised of the union of two geometries.
"""
struct UnionGeometry{G1,G2} <: AbstractGeometry where {G1<:AbstractGeometry, G2<:AbstractGeometry}
    geometry1::G1
    geometry2::G2
end

function ⊕(geometry1::G1, geometry2::G2) where {G1<:AbstractGeometry, G2<:AbstractGeometry}
    return UnionGeometry(geometry1, geometry2)
end

"""
    $TYPEDEF
"""
struct Mesh{G<:AbstractGeometry,M<:AbstractMaterial}
    geometry::G
    material::M
end
