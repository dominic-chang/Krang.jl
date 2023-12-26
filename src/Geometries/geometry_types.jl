"""
    $TYPEDEF

Abstract Geometry
"""
abstract type AbstractGeometry end

"""
    $TYPEDEF

Empty Geometry
"""
abstract type EmptyGeometry end


"""
    $TYPEDEF

Abstract Composite Geometry
"""
abstract type AbstractCompositeGeometry{M1,M2} <: AbstractGeometry end

"""
    $TYPEDEF

Sum of 2 Geometries
"""
struct AddGeometry{M1,M2} <: AbstractCompositeGeometry{M1,M2} 
    m1::M1
    m2::M2
end


"""
    $TYPEDEF

Cone Geometry with half opening angle `opening_angle`.
"""
struct ConeGeometry{T} <: AbstractGeometry
    opening_angle::T
end

Disk() = ConeGeometry(π/2)


"""
    $TYPEDEF

Sphere Geometry with radius `radius`.
"""
struct SphereGeometry{T} <: AbstractGeometry
    radius::T
end


