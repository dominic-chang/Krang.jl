"""
    $TYPEDEF

Abstract Geometry
"""
abstract type AbstractGeometry end

"""
    $TYPEDEF

Cone Geometry with half opening angle `opening_angle`.
"""
struct ConeGeometry{T} <: AbstractGeometry
    opening_angle::T
end

Disk() = ConeGeometry(π/2)
