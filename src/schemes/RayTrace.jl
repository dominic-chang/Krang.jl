export generate_ray, raytrace

struct RayTrace <: AbstractScheme end

function generate_ray!(
    ray::Vector{Intersection{T}},
    pixel::Krang.AbstractPixel,
    res::Int,
) where {T}
    actual_res = res
    τf = total_mino_time(pixel)

    Δτ = τf / actual_res
    for I in range(1, res)
        ts, rs, θs, ϕs, νr, νθ, _ = emission_coordinates(pixel, Δτ * I)
        ϕs = ϕs % T(2π)
        ray[I] = Intersection(ts, rs, θs, ϕs, νr, νθ)
    end
end

function generate_ray(pixel::AbstractPixel{T}, res::Int) where {T}
    ray = Vector{Intersection{T}}(undef, res)#zeros(T, 3, res)
    generate_ray!(ray, pixel, res)
    return ray
end

#TODO: debug this
@kernel function _generate_rays!(
    rays::AbstractArray{Intersection{T}},
    pixels::AbstractMatrix{P},
    res::Int,
) where {T,P<:AbstractPixel}
    (I, J, K) = @index(Global, NTuple)
    k = Int(K)
    pixel = pixels[I, J]
    actual_res =
        unsafe_trunc(Int, sum((Krang._isreal2.(Krang.roots(pixel))))) == 4 ? res + 1 : res
    τf = total_mino_time(pixel)

    Δτ = τf / actual_res

    ts, rs, θs, ϕs, νr, νθ = (2.0f0, 0.0f0, 0.0f0, 0.0f0, true, true)
    rays[I, J, K] = Intersection(ts, rs, θs, ϕs, νr, νθ)
    emission_coordinates(pixel, Δτ * k)


end
function generate_rays(
    pixels::AbstractMatrix{P},
    res::Int;
    A = Array,
) where {P<:AbstractPixel{T}} where {T}
    dims = (size(pixels)..., res)
    rays = A{Intersection{T}}(undef, dims...)
    backend = get_backend(rays)
    _generate_rays!(backend)(rays, pixels, res, ndrange = dims)
    return rays
end

function line_intersection(origin::NTuple{3,T}, line_point_2, t_a, t_b, t_c) where {T}
    e1 = (t_b[1] - t_a[1], t_b[2] - t_a[2], t_b[3] - t_a[3])
    e2 = (t_c[1] - t_a[1], t_c[2] - t_a[2], t_c[3] - t_a[3])
    s = (origin[1] - t_a[1], origin[2] - t_a[2], origin[3] - t_a[3])
    direction = (
        line_point_2[1] - origin[1],
        line_point_2[2] - origin[2],
        line_point_2[3] - origin[3],
    )
    dir_len = sqrt(dot(direction, direction))
    direction = direction ./ dir_len

    ray_cross_e2 = cross(direction, e2)
    det = dot(e1, ray_cross_e2)

    outarr = (zero(T), zero(T), zero(T))

    if (-eps(T) < det < eps(T))
        return false, outarr
    end#This ray is parallel to this triangle.

    inv_det = inv(det)
    u = inv_det * dot(s, ray_cross_e2)
    if !(one(T) > u > zero(T))
        return false, outarr
    end

    s_cross_e1 = cross(s, e1)
    v = inv_det * dot(direction, s_cross_e1)
    if ((v < zero(T)) || (u + v > one(T)))
        return false, outarr
    end
    # At this stage we can compute t to find out where the intersection point is on the line.
    t = inv_det * dot(e2, s_cross_e1)

    intersection_point = Tuple(origin .+ (direction .* t))
    #return Some(intersection_point);

    if t > zero(T)
        return t < dir_len, intersection_point
    end# ray intersection

    # This means that there is a line intersection but not a ray intersection.
    return false, outarr
end

function cross(a, b)
    return a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1]
end
function dot(a, b)
    return a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
end

function raytrace(pixel::AbstractPixel, mesh::Mesh; res = 100)
    return_trait = returnTrait(mesh.material)
    return raytrace(return_trait, pixel, mesh; res = res)
end

function raytrace(
    ::AbstractReturnTrait,
    pix::AbstractPixel{T},
    mesh::Mesh;
    res = 100,
) where {T}
    observation = zero(T)
    return _raytrace(observation, pix, mesh; res = res)
end

function raytrace(
    ::SimplePolarizationTrait,
    pix::AbstractPixel{T},
    mesh::Mesh;
    res = 100,
) where {T}
    observation = StokesParams(zero(T), zero(T), zero(T), zero(T))
    return _raytrace(observation, pix, mesh; res = res)
end

function raytrace(camera::AbstractCamera, mesh::Mesh{A}; res = 100) where {A}
    return_trait = returnTrait(mesh.material)
    return raytrace(return_trait, camera, mesh; res = res)
end

function raytrace(
    ::AbstractReturnTrait,
    camera::AbstractCamera,
    mesh::Mesh{A};
    res = 100,
) where {A}
    intersections =
        Array{typeof(camera.metric).parameters[1]}(undef, size(camera.screen.pixels)...)
    for I in CartesianIndices(camera.screen.pixels)
        intersections[I] = raytrace(camera.screen.pixels[I], mesh; res = res)
    end
    return intersections
end

function raytrace(
    ::SimplePolarizationTrait,
    camera::AbstractCamera,
    mesh::Mesh{A};
    res = 100,
) where {A}
    intersections = Array{StokesParams}(undef, size(camera.screen.pixels)...)
    for I in CartesianIndices(camera.screen.pixels)
        intersections[I] = raytrace(camera.screen.pixels[I], mesh; res = res)
    end
    return intersections
end



#https://journals.aps.org/prd/abstract/10.1103/PhysRevD.86.084049
function t_kerr_schild(metric::Kerr{T}, tBL, rBL) where {T}
    a = metric.spin
    temp = sqrt(1 - a^2)
    rp = 1 + temp
    rm = 1 - temp
    num = rBL - rp
    den = rBL - rm

    term1 = (rp^2 + a^2) * log(abs(num / den)) / (2temp)
    term2 = -2 * log(2 / (rBL - rm))
    return tBL + term1 + term2
end

function ϕ_kerr_schild(metric::Kerr{T}, rBL, ϕBL) where {T}
    a = metric.spin
    temp = sqrt(1 - a^2)
    rp = 1 + temp
    rm = 1 - temp
    num = rBL - rp + eps()
    den = rBL - rm + eps()

    term1 = a / (2temp) * log(abs(num / den))
    term2 = -atan(a, rBL)
    ans = ϕBL + term1 + term2
    if isinf(ans)
        @warn "ϕ_kerr_schild is inf at rs=$rBL. This usually happens if the ray intersects the horizon."
        return ϕBL + term2
    end
    return ans
end

function ϕ_BL(metric::Kerr{T}, rKS, ϕKS) where {T}
    a = metric.spin
    temp = sqrt(1 - a^2)
    rp = 1 + temp
    rm = 1 - temp
    num = rKS - rp + eps()
    den = rKS - rm + eps()

    term1 = -a / (2temp) * log(abs(num / den))
    term2 = atan(a, rKS)
    ans = ϕKS + term1 + term2
    if isinf(ans)
        @warn "ϕ_BL is inf at rs=$rKS. This usually happens if the ray intersects the horizon."
        return ϕKS + term2
    end
    return ans
end


"""    
    Transforms from Boyer-Lindquist to Kerr-Schild coordinates
"""
function boyer_lindquist_to_quasi_cartesian_kerr_schild(metric::Kerr, tBL, rBL, θBL, ϕBL)
    t = t_kerr_schild(metric, tBL, rBL)
    r = rBL
    θ = θBL
    ϕ = ϕ_kerr_schild(metric, rBL, ϕBL)
    sθ = sin(θ)
    sϕ = sin(ϕ)
    cϕ = cos(ϕ)
    return (t, sθ * (r * cϕ), sθ * (r * sϕ), r * cos(θ))
end

"""    
    Transforms from Boyer-Lindquist to Kerr-Schild coordinates ignoring time
"""
function boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
    metric::Kerr{T},
    rBL,
    θBL,
    ϕBL,
) where {T}
    a = metric.spin
    r = rBL
    θ = θBL
    ϕ = ϕ_kerr_schild(metric, rBL, ϕBL)
    sθ = sin(θ)
    sϕ = sin(ϕ)
    cϕ = cos(ϕ)

    return r .* (sθ * cϕ, sθ * sϕ, cos(θ))
end
