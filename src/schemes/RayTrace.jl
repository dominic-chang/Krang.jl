export generate_ray, raytrace

struct RayTrace <: AbstractScheme end

function generate_ray!(ray::Vector{Intersection{T}}, pixel::Krang.AbstractPixel, res::Int) where T
    actual_res = unsafe_trunc(Int, sum((Krang._isreal2.(Krang.roots(pixel))))) == 4 ? res+1 : res 
    τf = total_mino_time(pixel)

    Δτ = τf/actual_res
    for I in range(1, res)
        ts, rs, θs, ϕs, νr, νθ, _ = emission_coordinates(pixel, Δτ*I)
        ϕs = ϕs % T(2π)
        ray[I] = Intersection(ts, rs, θs, ϕs, νr, νθ)
    end
end

function generate_ray(pixel::AbstractPixel{T}, res::Int) where T
    ray = Vector{Intersection{T}}(undef,res)#zeros(T, 3, res)
    generate_ray!(ray, pixel, res)
    return ray
end

@kernel function _generate_rays!(rays::AbstractArray{Intersection{T}}, pixels::AbstractMatrix{P}, res::Int) where {T, P<:AbstractPixel}
    (I,J,K) = @index(Global, NTuple)
    k = Int(K)
    pixel = pixels[I,J]
    actual_res = unsafe_trunc(Int, sum((Krang._isreal2.(Krang.roots(pixel))))) == 4 ? res+1 : res 
    τf =  total_mino_time(pixel)

    Δτ = τf/actual_res

    ts, rs, θs, ϕs, νr, νθ = (2f0, 0f0, 0f0, 0f0, true, true)
    rays[I,J,K] = Intersection(ts, rs, θs, ϕs, νr, νθ)
    emission_coordinates(pixel, Δτ*k)

    #ts, rs, θs, ϕs, νr, νθ = (2f0, 0f0, 0f0, 0f0, true, true)
    #rays[I,J,K] = Intersection(ts, rs, θs, ϕs, νr, νθ)


end
function generate_rays(pixels::AbstractMatrix{P}, res::Int; A=Array)  where {P <: AbstractPixel{T}} where T
    dims = (size(pixels)..., res)
    rays = A{Intersection{T}}(undef, dims...)
    backend = get_backend(rays)
    _generate_rays!(backend)(rays, pixels, res, ndrange = dims)
    return rays
end

function line_intersection(origin::NTuple{3,T}, line_point_2, t_a, t_b, t_c) where T
    e1 = (t_b[1]-t_a[1], t_b[2]-t_a[2], t_b[3]-t_a[3])
    e2 = (t_c[1]-t_a[1], t_c[2]-t_a[2], t_c[3]-t_a[3])
	s = (origin[1] - t_a[1], origin[2] - t_a[2], origin[3] - t_a[3])
    direction = (line_point_2[1] - origin[1], line_point_2[2] - origin[2], line_point_2[3] - origin[3])
    dir_len = sqrt(dot(direction, direction))
    direction = direction ./ dir_len

	ray_cross_e2 = cross(direction, e2)
	det = dot(e1, ray_cross_e2)

    outarr = (zero(T), zero(T), zero(T))

	if ( -eps(T) < det < eps(T)) return false, outarr end#This ray is parallel to this triangle.

	inv_det = inv(det)
	u = inv_det * dot(s, ray_cross_e2)
	if !(one(T) > u > zero(T)) return false, outarr end

	s_cross_e1 = cross(s, e1)
	v = inv_det * dot(direction, s_cross_e1)
	if ((v < zero(T)) || (u + v > one(T) )) return false, outarr end
	# At this stage we can compute t to find out where the intersection point is on the line.
	t = inv_det * dot(e2, s_cross_e1)

	intersection_point = Tuple(origin .+ (direction .* t))
	#return Some(intersection_point);

	if t > zero(T) return t < dir_len, intersection_point end# ray intersection
        
	# This means that there is a line intersection but not a ray intersection.
	return false, outarr
end

function cross(a, b) 
    return a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1]
end
function dot(a,b)
    return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
end

