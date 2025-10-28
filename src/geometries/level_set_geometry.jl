"""
    Level Sets should be a functor
"""
abstract type AbstractLevelSetGeometry{T} <: AbstractGeometry end

@inline function _raytrace(
    observation::A,
    pixel::AbstractPixel{T},
    mesh::Mesh{<:AbstractLevelSetGeometry,<:AbstractMaterial};
    res = 100,
) where {A,T}
    geometry = mesh.geometry
    material = mesh.material
    τf = Krang.total_mino_time(pixel)
    Δτ = τf / res - eps()
    origin_bl = 1e100,pixel.θo, 0.0, true, true
    for i in 1:res
        τref = Δτ*i
        τ = Krang.line_intersection(pixel, origin_bl, Δτ*(i-1), τref, geometry)
        rs, θs, ϕs, νr, νθ = Krang.emission_coordinates_fast_light(pixel, τ)

        if rs < Krang.horizon(pixel.metric)
            continue
        end

        intersection = Krang.Intersection(0.0, rs, θs, ϕs, νr, νθ)
        observation += τ == 0 ? T(0) : material(pixel, intersection)#didintersect ? 1 : 0
        origin_bl = Krang.emission_coordinates_fast_light(pixel, τref)

    end
    return observation
end

@inline function line_intersection(
    pixel::Krang.AbstractPixel,
    origin_bl,
    τi::T,
    τf::T,
    geometry::Krang.AbstractLevelSetGeometry{T},
) where {T}
    origin = Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(pixel.metric, origin_bl[1:3]...)
    line_point_2_bl = Krang.emission_coordinates_fast_light(pixel, τf)
    line_point_2 = Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(pixel.metric, line_point_2_bl[1:3]...) 
    geo_origin = geometry(origin...)
    didintersect = geo_origin * geometry(line_point_2...) <= zero(T)
    if didintersect
        t = Krang.Roots.find_zero((x) -> x == zero(T) ? geo_origin : geometry(Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(pixel.metric, Krang.emission_coordinates_fast_light(pixel, x)[1:3]...)...), (τi, τf))
        return t
    end
    return zero(T)
end