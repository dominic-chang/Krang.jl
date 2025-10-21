"""
    Level Sets should be a functor
"""
abstract type AbstractLevelSetGeometry{T} <: AbstractGeometry end

function _raytrace(
    observation::A,
    pixel::AbstractPixel{T},
    mesh::Mesh{<:AbstractLevelSetGeometry,<:AbstractMaterial};
    res = 100,
) where {A,T}
    geometry = mesh.geometry
    material = mesh.material
    ray = Vector{Intersection{T}}(undef, res)
    generate_ray!(ray, pixel, res)
    (; rs, θs, ϕs, νr, νθ) = ray[end]
    origin =
        boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(pixel.metric, rs, θs, ϕs)
    z = zero(A)
    for i = res:-1:1
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

        intersection = Intersection(ts, rs, θs, ϕs, νr, νθ)
        observation += !didintersect ? z : material(pixel, intersection)#didintersect ? 1 : 0
        origin = line_point_2
    end
    return observation
end

@inline function line_intersection(
    origin::NTuple{3,T},
    line_point_2,
    geometry::AbstractLevelSetGeometry{T},
) where {T}
    didintersect = geometry(origin...) * geometry(line_point_2...) <= zero(T)
    if didintersect
        direction = (
            line_point_2[1] - origin[1],
            line_point_2[2] - origin[2],
            line_point_2[3] - origin[3],
        )

        t = find_zero((x) -> geometry((origin .+ (direction .* x))...), (zero(T), one(T)))

        return true, origin .+ (direction .* t)
    end
    return false, (zero(T), zero(T), zero(T))

end
