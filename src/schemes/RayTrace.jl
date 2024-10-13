export generate_ray, raytrace

struct RayTrace <: AbstractScheme end

function generate_ray!(ray::Matrix{T}, pixel::Krang.AbstractPixel, res::Int) where T
    actual_res = unsafe_trunc(Int, sum((Krang._isreal2.(Krang.roots(pixel))))) == 4 ? res+1 : res 
    τf = total_mino_time(pixel)

    Δτ = τf/actual_res
    for I in range(1, res)

        _, curr_rad, curr_inc, curr_az, _, _ = emission_coordinates(pixel, Δτ*I)
        curr_az = curr_az % T(2π)
        ray[1,I] = curr_rad * sin(curr_inc)*cos(curr_az)
        ray[2,I] = curr_rad * sin(curr_inc)*sin(curr_az)
        ray[3,I] = curr_rad * cos(curr_inc)

    end
end

function generate_ray(pixel::AbstractPixel{T}, res::Int) where T
    ray = zeros(T, 3, res)
    generate_ray!(ray, pixel, res)
    return ray
end

@kernel function _generate_rays!(rays::AbstractArray{T}, pixels::AbstractMatrix{P}, res::Int) where {T, P<:AbstractPixel}
    (I,J,K,L) = @index(Global, NTuple)
    l = Int(L)
    pixel = pixels[I,J]
    actual_res = unsafe_trunc(Int, sum((Krang._isreal2.(Krang.roots(pixel))))) == 4 ? res+1 : res 
    τf =  total_mino_time(pixel)

    Δτ = τf/actual_res

    _, curr_rad, curr_inc, curr_az, _, _ = emission_coordinates(pixel, Δτ*l)
    curr_az = curr_az % T(2π)
    k = Int(K)
    if k == 1
        rays[I,J,K,L] = curr_rad * sin(curr_inc)*cos(curr_az)
    elseif k == 2
        rays[I,J,K,L] = curr_rad * sin(curr_inc)*sin(curr_az)
    else
        rays[I,J,K,L] = curr_rad * cos(curr_inc)
    end
end

function generate_rays(pixels::AbstractMatrix{P}, res::Int; A=Array)  where {P<: AbstractPixel{T}} where T
    dims = (size(pixels)..., 3, res)
    rays = A{T}(undef, dims...)
    backend = get_backend(rays)
    _generate_rays!(backend)(rays, pixels, res, ndrange = dims)
    return rays
end

function line_intersection(origin::AbstractVector{T}, line_point_2, t_a, t_b, t_c) where T
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

function raytrace(camera, mesh_geometry; res=100) 
    faces = begin
        temp = getfield(getfield(mesh_geometry, :simplices), :faces)
        len = length(temp)
        reshape([i[j] for i in temp for j in 1:3], 3, len)
    end
    vertices = begin
        temp = mesh_geometry.position
        len = length(temp)
       reshape([(i[j]) for i in temp for j in 1:3], 3, len)
    end
    intersections = zeros(Int, size(camera.screen.pixels))
    Threads.@threads for I in CartesianIndices(camera.screen.pixels)
        intersections[I] = raytrace(camera.screen.pixels[I], faces, vertices;res)
    end
    return intersections
end
function raytrace(pixel::AbstractPixel{T}, faces::Matrix{OffsetInteger{-1, UInt32}}, vertices::Matrix{T}; res=100) where T
    intersections = 0
    ray = zeros(T, 3, res)
    generate_ray!(ray, pixel, res)
    origin = @view ray[:,1]
    for i in 2:res
        line_point_2 = @view ray[:,i]
        for j in 1:(size(faces)[2])
            f1, f2, f3 = @view faces[:,j]
            v1 = @view vertices[:,f1]; 
            v2 = @view vertices[:,f2]; 
            v3 = @view vertices[:,f3];
            didintersect, point = line_intersection(origin, line_point_2, v1, v2, v3)
            intersections += didintersect ? 1 : 0
        end
        origin = line_point_2
    end
    return intersections
end