using Revise
using Krang, GeometryBasics, FileIO
import CairoMakie as CMk
using Rotations
using Metal
using Adapt
using StaticArrays

spin = 0.9f0
metric = Krang.Kerr(spin)
θo = 20f0/180*π
ρmax = 20f0
sze = 50

# Import external mesh
using FileIO
bunny_mesh = load(joinpath((@__DIR__),"stanfordbunny.obj") )
rot = Rotations.AngleAxis(π / 2, 0.0, 0.0, 1.0)
points = (Ref(rot) .* (bunny.position .* 150)) .+   Ref([15, 3, 0])
faces = getfield(getfield(bunny, :simplices), :faces)
bunny_mesh = GeometryBasics.Mesh([Point(Float32.(x)) for x in points], faces)
bunny_mesh = (bunny_mesh...)

#GLMk.mesh!(ax, mag_mesh, color=red_cb, linewidth=0.01, shading=false, overdraw=false)
CMk.mesh(
    bunny_mesh,
    color=[parse(CMk.Colorant, "rgba(84%, 11%, 38%,1.0)") for tri in bunny_mesh.position],
    colormap = :blues, transparency=true,
    shading=false
)

# Create screen
camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
ray_segment = Krang.emission_coordinates.(Ref(camera.screen.pixels[1]), [0.1f0, 10f9])

line = Line(Point(0,0,-1), Point(0,0,1))
simplex_points1 = (Point(1,0,0.0), Point(0,1,0.0), Point(-1,-1,0.0))
simplex1 = Triangle(simplex_points1...)
simplex_points2 = (Point(1,0,5), Point(0,1,5), Point(-1,-1,5))
simplex2 = Triangle(simplex_points2...)


function cross(a::Point{T}, b::Point{T}) where T
    #return [zero(T) -a[3] a[2]; a[3] zero(T) -a[1]; -a[2] a[1] zero(T)] * b
    return Point(a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1])
end

#Möller–Trumbore_intersection_algorithm
#https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
function line_intersection(line::Line, triangle::Triangle{3,T}) where T
    origin, lp2 = line.points
    direction = lp2 - origin
    dir_len = sqrt(direction' * direction)
    direction /= dir_len
    t_a, t_b, t_c = triangle.points
	e1 = t_b - t_a
	e2 = t_c - t_a

	ray_cross_e2 = cross(direction, e2)
	det = e1' * ray_cross_e2

	if ( -eps(T) < det < eps(T)) return false end#This ray is parallel to this triangle.

	inv_det = inv(det)
	s = origin - t_a
	u = inv_det * (s' * ray_cross_e2)
	if !(one(T) > u > zero(T)) return false end

	s_cross_e1 = cross(s, e1)
	v = inv_det * (direction' * s_cross_e1)
	if ((v < zero(T)) || (u + v > one(T) )) return false end
	# At this stage we can compute t to find out where the intersection point is on the line.
	t = inv_det * (e2' * s_cross_e1)

		#intersection_point = origin + (direction * t)
		#return Some(intersection_point);

	if t > zero(T) return t < dir_len end# ray intersection
        
	# This means that there is a line intersection but not a ray intersection.
	return false
end

using Polyester
function get_trajectory(pixel::Krang.IntensityPixel{T}, res, n) where T
    curr_rad, curr_inc, curr_az = T(Inf), Krang.inclination(pixel), zero(T)#Krang.emission_coordinates_fast_light(pixel, eps(T), true, n)
    prev_x = curr_rad * sin(curr_inc)*cos(curr_az)
    prev_y = curr_rad * sin(curr_inc)*sin(curr_az)
    prev_z = curr_rad *  cos(curr_inc)
    θo = Krang.inclination(pixel)
    lines = []
    for θs in 0:T(π/res):T(π)
        if (θs == θo || (π-θs) == θo)  continue end
        curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(pixel, θs, Krang.screen_coordinate(pixel)[2] ≥ 0, n)# .+ Krang.emission_coordinates_fast_light(pixel, θs, false, n)
        if iszero(curr_rad) continue end
        curr_x = (curr_rad * sin(curr_inc)*cos(curr_az))
        curr_y = (curr_rad * sin(curr_inc)*sin(curr_az))
        curr_z = (curr_rad * cos(curr_inc))
        curr_line = Line(Point(prev_x, prev_y, prev_z), Point(curr_x, curr_y, curr_z))
        
        append!(lines, curr_line)
        prev_x, prev_y, prev_z = curr_x, curr_y, curr_z
    end

    return lines
end

#function trace_pixel!(pixel::Krang.IntensityPixel{T}, bunny_mesh::NTuple{N,Triangle{3, T}}, res::Int, n::Int) where {N,T}
#function trace_pixel!(pixel::Krang.IntensityPixel{T}, bunny_mesh::NTuple{N,Triangle{3, T}}, res::Int, n::Int) where {N,T}
function trace_pixel!(pixel::Krang.IntensityPixel{T}, bunny_mesh::AbstractArray{Triangle{3,T}}, res::Int, n::Int) where {T}
#function trace_pixel!(pixel::Krang.IntensityPixel{T}, bunny_mesh::Triangle, res::Int, n::Int) where {T}
    curr_rad = T(Inf)
    curr_inc = Krang.inclination(pixel) 
    curr_az = zero(T)
    prev_x = curr_rad * sin(curr_inc)*cos(curr_az)
    prev_y = curr_rad * sin(curr_inc)*sin(curr_az)
    prev_z = curr_rad *  cos(curr_inc)
    tot = zero(T)
    θo = Krang.inclination(pixel)
    θs = zero(T)
    #for θs in 0:T(π/T(res)):T(π)
    while θs < T(π)
        θs = θs + T(π/T(res))

        if (θs == θo || (T(π)-θs) == θo)  continue end

        curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(pixel, θs, Krang.screen_coordinate(pixel)[2]>0, n)

        if iszero(curr_rad) continue end

        curr_x = (curr_rad * sin(curr_inc)*cos(curr_az))
        curr_y = (curr_rad * sin(curr_inc)*sin(curr_az))
        curr_z = (curr_rad * cos(curr_inc))
        curr_line = Line(Point(prev_x, prev_y, prev_z), Point(curr_x, curr_y, curr_z))
        
        tot += mapreduce(x::Triangle{3,T}-> line_intersection(curr_line, x) ? one(T) : zero(T), +, bunny_mesh)
        #tot  += curr_az
        #tot += line_intersection(curr_line, bunny_mesh)

        prev_x, prev_y, prev_z = curr_x, curr_y, curr_z
    end

    return tot
end

using LaTeXStrings
fig = CMk.Figure(resolution=(900,300));
obs_screen = zeros(Float32, size(camera.screen.pixels))
#obs_screen = zeros(size(camera.screen.pixels))

using Metal
n = 0
line = Line(Point(-10f0,0f0,0f0), Point(10f0,0f0,0f0))
mtl_bunny_mesh = adapt(MtlArray, bunny_mesh)
buny_mesh = ones(Float32, 100)
foo(a, b) = b 
a = foo.(MtlArray(ones(Float32, 100)), Ref(bunny_mesh))

a = Array(trace_pixel!.(MtlArray(camera.screen.pixels), Ref(Vector(bunny_mesh)), 72, n))
b = collect(map(pix->trace_pixel!(pix, bunny_mesh, 72, n), (camera.screen.pixels)))
#b = trace_pixel!.((camera.screen.pixels), Ref(Tuple(bunny_mesh[begin:begin+1])), 72, n)
maximum(abs.(a .- b))
CMk.heatmap(a)
CMk.heatmap(b)
#mapreduce(x->trace_pixel!.(MtlArray(camera.screen.pixels), Ref(x), 72, n), +, (bunny_mesh))

for n in 0:2
    #@time Threads.@threads for I in CartesianIndices(camera.screen.pixels)
        #obs_screen[I] = trace_pixel(camera.screen.pixels[I], bunny_mesh, 72, 0  )
    #obs_screen = Array(trace_pixel!.((camera.screen.pixels), Ref(bunny_mesh), 72, n  ))
    obs_screen = collect(mapreduce(x->trace_pixel!.((camera.screen.pixels), Ref(Tuple(x)), 72, n), +, (bunny_mesh[begin:begin+2000], bunny_mesh[begin+2001:begin+4000], bunny_mesh[begin+4001:end])))


    #max(obs_screen...)

    ax = CMk.Axis(fig[1,n+1], aspect=1, title="n=$n")
    CMk.hidedecorations!(ax)
    CMk.heatmap!(ax,obs_screen)
end
display(fig)

met = Krang.Kerr(spin)
pix = Krang.IntensityPixel(met, 5.2f0, 2f0, θo)

begin
    fig = CMk.Figure()
    ax = CMk.Axis3(fig[1,1])
    CMk.xlims!(ax, (-2_0, 2_0)) 
    CMk.ylims!(ax, (-2_0, 2_0)) 
    CMk.zlims!(ax, (-2_0, 2_0)) 
    lines_to_plot = []
    lines_to_plot = get_trajectory(pix, 10_000,0)
    #println(lines_to_plot[end-10])
    append!(lines_to_plot, get_trajectory(pix, 10_000,1)[end:-1:begin])
    append!(lines_to_plot, get_trajectory(pix, 10_000,2))
    filter!(x->x[1] < Inf, lines_to_plot)
    lines_to_plot = lines_to_plot#[begin:end-150]

    CMk.lines!(ax, lines_to_plot, color=range(1, length(lines_to_plot)))
    display(fig)
end
