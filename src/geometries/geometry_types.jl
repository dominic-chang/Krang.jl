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

#TODO: replace with trait
isOcclusive(material::AbstractMaterial) = false
#TODO: replace with trait
isAxisymmetric(material::AbstractMaterial) = false
#TODO: replace with trait
isFastLight(material::AbstractMaterial) = false

struct Intersection{T}
    ts::T
    rs::T
    θs::T
    ϕs::T
    νr::Bool
    νθ::Bool
end

yield(material::AbstractMaterial) = 0

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

function raytrace(pix::AbstractPixel{T}, mesh::Mesh{<:ConeGeometry{T,A}, <:AbstractMaterial}) where {T,A}
    #(;magnetic_field, fluid_velocity, spectral_index, R, p1, p2, subimgs) = linpol
    
    geometry = mesh.geometry
    material = mesh.material
    θs = geometry.opening_angle
    subimgs= material.subimgs

    observation = zero(T)#StokesParams(zero(T), zero(T), zero(T), zero(T))

    isindir = false
    for _ in 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        for n in subimgs
            #νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
            intersection, issuccess = if isFastLight(material)
                if isAxisymmetric(material)
                    rs, νr, νθ, _, issuccess = emission_radius(pix, θs, isindir, n)
                    Intersection(zero(rs), rs, θs, zero(rs), νr, νθ), issuccess
                else
                    rs, ϕs, νr, νθ, issuccess = emission_coordinates_fast_light(pix, θs, isindir, n)
                    Intersection(zero(rs), rs, θs, ϕs, νr, νθ), issuccess
                end
            else
                ts, rs, ϕs, νr, νθ, issuccess = emission_coordinates(pix, θs, isindir, n)
                Intersection(ts, rs, θs, ϕs, νr, νθ), issuccess
            end
            if issuccess && (horizon(metric(pix)) < rs < T(Inf))

                if isOcclusive(material)
                    observation = material(pix, intersection)
                else
                    observation += material(pix, intersection)
                end
            end
        end
    end
    return observation
end


Disk(;attributes=nothing) = ConeGeometry(π/2;attributes=attributes)

