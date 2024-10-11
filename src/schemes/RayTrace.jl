export generate_ray, raytrace

struct RayTrace <: AbstractScheme end

function generate_ray!(ray::Matrix{T}, pixel::Krang.AbstractPixel, res::Int) where T
    numreals = Int(sum((Krang._isreal2.(Krang.roots(pixel)))))
    τf = zero(T)
    actual_res = res
    if numreals == 4 
        τf = 2Krang.I0_inf(pixel)
        actual_res += 1 # The last point on scattering rays is at infinity
    else
        rh = Krang.horizon(pixel.metric)
        radial_roots = roots(pixel)
        r1, r2, _, r4 = radial_roots
        r21, r31, r32, r41, r42, _ = _get_root_diffs(radial_roots...)
        r1, r2, r21 = real.((r1, r2, r21))

        if numreals == 2
            A = √abs(r32 * r42)
            B = √abs(r31 * r41)
            k = (((A + B)^2 - r21^2) / (4 * A * B))
            temprat = B * (rh - r2) / (A * (rh - r1))
            x3_s = clamp(((one(T) - temprat) / (one(T) + temprat)), -one(T), one(T))
            coef = one(T) * √inv(A * B)
            τf = I0_inf(pixel) - coef * JacobiElliptic.F((acos(x3_s)), k)
        else
            C = √abs(r31 * r42)
            D = √abs(r32 * r41)
            k4 = 4C * D / (C + D)^2
            a2 = abs(imag(r1))
            b1 = real(r4)
        
            k4 = T(4) * C * D / (C + D)^2
        
            go = √max((T(4)a2^2 - (C - D)^2) / ((C + D)^2 - T(4)a2^2), zero(T))
            x4_s = (rh + b1) / a2
            coef = 2 / (C + D)
            τf = I0_inf(pixel) - coef*JacobiElliptic.F(atan(x4_s) + atan(go), k4)
        end
    end

    Δτ = T(τf/(actual_res))
    for I in range(1, res)

        _, curr_rad, curr_inc, curr_az, _, _ = emission_coordinates(pixel, Δτ*I)
        curr_az = mod2pi(curr_az)
        ray[1,I] = curr_rad * sin(curr_inc)*cos(curr_az)
        ray[2,I] = curr_rad * sin(curr_inc)*sin(curr_az)
        ray[3,I] = curr_rad * cos(curr_inc)

    end

end

function generate_ray(pixel::Krang.AbstractPixel{T}, res::Int) where T
    ray = zeros(T, 3, res)
    generate_ray!(ray, pixel, res)
    return ray
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