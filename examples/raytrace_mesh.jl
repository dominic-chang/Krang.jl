using Krang, GeometryBasics, FileIO
import CairoMakie as CMk
using Rotations

spin = 0.9f0
metric = Krang.Kerr(spin)
θo = Float32(20/180*π)
ρmax = 20f0
sze = 500

# Magnetic Field Lines
bunny = load("/Users/dominicchang/Software/Krang.jl/examples/stanfordbunny.obj") 
#maglines = load((@__DIR__) * "/Magnetic_Field_OBJ/maglines1.obj")
rot = Rotations.AngleAxis(π / 2, 0.0, 0.0, 1.0)
points = Ref(rot) .* (bunny.position .* 150)
points .+=   Ref([15.0, 3.0, 0.])
faces = getfield(getfield(bunny, :simplices), :faces)
bunny_mesh = GeometryBasics.Mesh(map(x->Point(Float32.(x)),points), faces)

#GLMk.mesh!(ax, mag_mesh, color=red_cb, linewidth=0.01, shading=false, overdraw=false)
CMk.mesh(
    bunny_mesh,
    color=[parse(CMk.Colorant, "rgba(84%, 11%, 38%,1.0)") for tri in bunny_mesh.position],
    colormap = :blues, transparency=true,
    shading=false
)

# Create screen
camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
ray_segment = Krang.raytrace.(Ref(camera.screen.pixels[1]), [0.1f0, 10f9])

line = Line(Point(0,0,-1), Point(0,0,1))
simplex_points1 = GeometryBasics.@SVector [Point(1,0,0.0), Point(0,1,0.0), Point(-1,-1,0.0)]
simplex1 = Triangle(simplex_points1...)
simplex_points2 = GeometryBasics.@SVector [Point(1,0,5), Point(0,1,5), Point(-1,-1,5)]
simplex2 = Triangle(simplex_points2...)

function calc_normal(t::Triangle)
    v1, v2, _ = t.points
    normal = Point(v1[2]*v2[3] - v1[3]*v2[2], v1[3]*v2[1] - v1[1]*v2[3], v1[1]*v2[2] - v1[2]*v2[1])
    return normal# ./ √(sum(normal .^ 2))
end

function cross(a::Point{T}, b::Point{T}) where T
    #return [zero(T) -a[3] a[2]; a[3] zero(T) -a[1]; -a[2] a[1] zero(T)] * b
    return Point(a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1])
end

#Möller–Trumbore_intersection_algorithm
#https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
function line_intersection(line::Line, triangle::Triangle)
    origin, lp2 = line.points
    direction = lp2 - origin
    dir_len = sqrt(direction' * direction)
    direction /= dir_len
    t_a, t_b, t_c = triangle.points
	e1 = t_b - t_a
	e2 = t_c - t_a

	ray_cross_e2 = cross(direction, e2)
	det = e1' * ray_cross_e2

	( det ≈ 0.0) && return false #This ray is parallel to this triangle.

	inv_det = 1.0 / det
	s = origin - t_a
	u = inv_det * (s' * ray_cross_e2)
	!(1.0 > u > 0.0) && return false

	s_cross_e1 = cross(s, e1)
	v = inv_det * (direction' * s_cross_e1)
	if v < 0.0 || u + v > 1.0 
		return false
    end
	# At this stage we can compute t to find out where the intersection point is on the line.
	t = inv_det * (e2' * s_cross_e1)

	if t > 0.0 # ray intersection
		#intersection_point = origin + (direction * t)
		#return Some(intersection_point);
        return t < dir_len
	else  # This means that there is a line intersection but not a ray intersection.
		return false
    end
end

#function line_intersection(line::Line, simplex::Triangle)
#    lx1, ly1, lz1 = line.points[2] - line.points[1]
#    
#    (sx1, sy1, sz1) = simplex.points[2] - simplex.points[1]
#    (sx2, sy2, sz2) = simplex.points[3] - simplex.points[1]
#
#    lx0, ly0, lz0 = line.points[1]
#    sx0, sy0, sz0 = simplex.points[1]
#    a = (-ly0*lz1*sx2 + lz1*sx2*sy0 - lx1*lz0*sy2 + lx0*lz1*sy2 - 
#      lz1*sx0*sy2 + lx1*sy2*sz0 + lx1*ly0*sz2 - lx1*sy0*sz2 + 
#      ly1*(lz0*sx2 - sx2*sz0 - lx0*sz2 + sx0*sz2))/(-lz1*sx2*sy1 + 
#      lz1*sx1*sy2 + ly1*sx2*sz1 - lx1*sy2*sz1 - ly1*sx1*sz2 + 
#      lx1*sy1*sz2)
#    b =  (ly0*lz1*sx1 - lz1*sx1*sy0 + lx1*lz0*sy1 - lx0*lz1*sy1 + 
#        lz1*sx0*sy1 - lx1*sy1*sz0 - lx1*ly0*sz1 + lx1*sy0*sz1 + 
#        ly1*(-lz0*sx1 + sx1*sz0 + lx0*sz1 - sx0*sz1))/(-lz1*sx2*sy1 + 
#        lz1*sx1*sy2 + ly1*sx2*sz1 - lx1*sy2*sz1 - ly1*sx1*sz2 + 
#        lx1*sy1*sz2)
#    t =  (lz0*(sx2*sy1 - sx1*sy2) - (sx2*sy1 + sx1*sy2)*sz0 - 
#        ly0*sx2*sz1 + sx2*sy0*sz1 + lx0*sy2*sz1 - sx0*sy2*sz1 + 
#        ly0*sx1*sz2 - sx1*sy0*sz2 - lx0*sy1*sz2 + 
#        sx0*sy1*sz2)/(-lz1*sx2*sy1 + lz1*sx1*sy2 + ly1*sx2*sz1 - 
#        lx1*sy2*sz1 - ly1*sx1*sz2 + lx1*sy1*sz2)
#    return (0 ≤ a ≤ 1) && (0 ≤ b ≤ 1) && (0 ≤ t ≤ 1) 
#end

function straight_line(pix, τ)
    α, β = Krang.screen_coordinate(pix)
    return -β, -α, 20. -τ
end

function test(line, bunny_mesh)
    for _ in 1:100
        for I in LinearIndices(bunny_mesh) 
            line_intersection(line, bunny_mesh[I])
        end
    end
end
#@profview test(line, bunny_mesh)

using Polyester
function trace!(obs_screen::Matrix{T}, pixels::Matrix{Krang.IntensityPixel{T}}, bunny_mesh::GeometryBasics.Mesh, res, n) where T
    for I in CartesianIndices(pixels)
        curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(camera.screen.pixels[I], zero(T), true, n)
        prev_x = curr_rad * sin(curr_inc)*cos(curr_az)
        prev_y = curr_rad * sin(curr_inc)*sin(curr_az)
        prev_z = curr_rad *  cos(curr_inc)
        #prev_x, prev_y, prev_z = straight_line(pixels[I], 0)
        for θs in T(π/res):T(π/res):T(π)
        #for τ in 0:1:40
            curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(camera.screen.pixels[I], θs, true, n)
            if isnan(curr_rad) continue end
            curr_x = curr_rad * sin(curr_inc)*cos(curr_az)
            curr_y = curr_rad * sin(curr_inc)*sin(curr_az)
            curr_z = curr_rad * cos(curr_inc)
            #curr_x, curr_y, curr_z = straight_line(pixels[I], τ)
            curr_line = Line(Point(prev_x::T, prev_y::T, prev_z::T), Point(curr_x::T, curr_y::T, curr_z::T))
            tot = zero(T)
            @batch for J in LinearIndices(bunny_mesh)
                if line_intersection(curr_line, bunny_mesh[J])
                    obs_screen[I] += one(T)
                end
            end
            prev_x, prev_y, prev_z = curr_x, curr_y, curr_z

        end
        curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(camera.screen.pixels[I], zero(T), false, n)
        prev_x = curr_rad * sin(curr_inc)*cos(curr_az)
        prev_y = curr_rad * sin(curr_inc)*sin(curr_az)
        prev_z = curr_rad *  cos(curr_inc)
        #prev_x, prev_y, prev_z = straight_line(pixels[I], 0)
        for θs in T(π/res):T(π/res):T(π)
        #for τ in 0:1:40
            curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(camera.screen.pixels[I], θs, false, n)
            if isnan(curr_rad) continue end
            curr_x = curr_rad * sin(curr_inc)*cos(curr_az)
            curr_y = curr_rad * sin(curr_inc)*sin(curr_az)
            curr_z = curr_rad * cos(curr_inc)
            #curr_x, curr_y, curr_z = straight_line(pixels[I], τ)
            curr_line = Line(Point(prev_x, prev_y, prev_z), Point(curr_x, curr_y, curr_z))
            tot = zero(T)
            @batch for J in LinearIndices(bunny_mesh)
                if line_intersection(curr_line, bunny_mesh[J])
                    obs_screen[I] += one(T)
                end
            end

            prev_x, prev_y, prev_z = curr_x, curr_y, curr_z

        end
    end
end

function trace_pixel(pixel::Krang.IntensityPixel{T}, bunny_mesh::GeometryBasics.Mesh, res, n) where T
    curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(pixel, eps(T), true, n)
    prev_x = curr_rad * sin(curr_inc)*cos(curr_az)
    prev_y = curr_rad * sin(curr_inc)*sin(curr_az)
    prev_z = curr_rad *  cos(curr_inc)
    tot = zero(T)
    θo = Krang.inclination(pixel)
    for θs in 0:T(π/res):T(π)
        if θs == θo continue end
        curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(pixel, θs, true, n)
        if isnan(curr_rad) continue end
        curr_x = (curr_rad * sin(curr_inc)*cos(curr_az))
        curr_y = (curr_rad * sin(curr_inc)*sin(curr_az))
        curr_z = (curr_rad * cos(curr_inc))
        curr_line = Line(Point(prev_x, prev_y, prev_z), Point(curr_x, curr_y, curr_z))
        
        for J in LinearIndices(bunny_mesh)
            if line_intersection(curr_line, bunny_mesh[J])
                tot += one(T)
            end
        end
        prev_x, prev_y, prev_z = curr_x, curr_y, curr_z
    end

    curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(pixel, eps(T), false, n)
    prev_x = curr_rad * sin(curr_inc)*cos(curr_az)
    prev_y = curr_rad * sin(curr_inc)*sin(curr_az)
    prev_z = curr_rad *  cos(curr_inc)
    for θs in 0:T(π/res):T(π)
        if θs == θo continue end

        curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(pixel, θs, false, n)
        if isnan(curr_rad) continue end
        curr_x = (curr_rad * sin(curr_inc)*cos(curr_az))
        curr_y = (curr_rad * sin(curr_inc)*sin(curr_az))
        curr_z = (curr_rad * cos(curr_inc))
        curr_line = Line(Point(prev_x, prev_y, prev_z), Point(curr_x, curr_y, curr_z))
        for J in LinearIndices(bunny_mesh)
            if line_intersection(curr_line, bunny_mesh[J])
                tot += one(T)
            end
        end
        prev_x, prev_y, prev_z = curr_x, curr_y, curr_z
    end
    return tot
end

fig = CMk.Figure();
for n in 0:2
    obs_screen = zeros(Float32, size(camera.screen.pixels))
    #@time trace!(obs_screen, camera.screen.pixels, bunny_mesh,180, 10)
    using Metal
    #f(pix) = trace_pixel(pix, (bunny_mesh), 36,0)
    #f.(camera.screen.pixels)
    obs_screen = zeros(size(camera.screen.pixels))
    @time Threads.@threads for I in CartesianIndices(camera.screen.pixels)
        #obs_screen[I] = trace_pixel(camera.screen.pixels[I], bunny_mesh, 72, 0  )
        obs_screen[I] = trace_pixel(camera.screen.pixels[I], bunny_mesh, 36, n  )
    end
    max(obs_screen...)

    ax = CMk.Axis(fig[1,n+1], aspect=1)
    CMk.hidedecorations!(ax)
    CMk.heatmap!(ax,obs_screen)
end
display(fig)