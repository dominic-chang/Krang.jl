"""
    Level Sets should be a functor
"""
abstract type AbstractLevelSetGeometry{T} <: AbstractGeometry end

@inline function _raytrace(
    observation::A,
    pixel::AbstractPixel,
    mesh::Mesh{<:AbstractLevelSetGeometry,<:AbstractMaterial};
    res = 100,
) where {A}
    geometry = mesh.geometry
    material = mesh.material
    T = typeof(pixel.metric.spin)
    ray = Vector{Intersection{T}}(undef, res)
    generate_ray!(ray, pixel, res)
    (; rs, θs, ϕs, νr, νθ) = ray[1]
    origin =
        boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(pixel.metric, rs, θs, ϕs)
    z = zero(A)
    for i = 1:res
        (; ts, rs, θs, ϕs, νr, νθ) = ray[i]
        if rs <= Krang.horizon(pixel.metric) || iszero(rs)
            continue
        end
        line_point_2 = boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
            pixel.metric,
            rs,
            θs,
            ϕs,
        )
        if isinf(line_point_2[1]) || isinf(line_point_2[2]) || isinf(line_point_2[3])
            continue
        end
        didintersect, point = line_intersection(origin, line_point_2, geometry)
        rs = sqrt(sum(point .^ 2))
        θs = acos(point[3] / rs)
        ϕs = ϕ_BL(pixel.metric, rs, atan(point[2], point[1]))

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