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

Disk() = ConeGeometry(π/2)

"""
    $TYPEDEF

Geometry that is comprised of the union of two geometries.
"""
struct UnionGeometry <: AbstractGeometry
    geometry1::AbstractGeometry
    geometry2::AbstractGeometry
end

function ⊕(geometry1::AbstractGeometry, geometry2::AbstractGeometry) 
    return UnionGeometry(geometry1, geometry2)
end

"""
    $TYPEDEF
"""
struct Mesh 
    geometry::AbstractGeometry
    material::AbstractMaterial
end

