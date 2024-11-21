"""
    Level Sets should be a functor
"""
abstract type AbstractLevelSetGeometry{T} <: AbstractGeometry end

function raytrace(camera::AbstractCamera, mesh::Mesh{<:AbstractLevelSetGeometry, <:AbstractMaterial}; res=100) 
    intersections = Array{typeof(yield(mesh.material))}(undef, size(camera.screen.pixels)...)#zeros(yield(mesh.material), size(camera.screen.pixels))
    Threads.@threads for I in CartesianIndices(camera.screen.pixels)
        intersections[I] = raytrace(camera.screen.pixels[I], mesh;res)
    end
    return intersections
end

function raytrace(pixel::AbstractPixel{T}, mesh::Mesh{<:AbstractLevelSetGeometry, <:AbstractMaterial}; res=100) where T
    geometry = mesh.geometry
    material = mesh.material
    intersections = yield(material)
    ray = Vector{Intersection{T}}(undef, res)
    generate_ray!(ray, pixel, res)
    (;rs, θs, ϕs, νr, νθ) =  ray[1]
    origin = (rs * sin(θs)*cos(ϕs), rs * sin(θs)*sin(ϕs), rs * cos(θs))
    z = yield(material)
    for i in 1:res
        (;ts, rs, θs, ϕs, νr, νθ) =  ray[i]
        line_point_2 = (rs * sin(θs)*cos(ϕs), rs * sin(θs)*sin(ϕs), rs * cos(θs))
        if isinf(line_point_2[1]) || isinf(line_point_2[2]) || isinf(line_point_2[3])
            continue
        end
        didintersect, point = line_intersection(origin, line_point_2, geometry)
        rs = sqrt(sum(point .^ 2))
        θs = acos(point[3]/rs)
        ϕs = atan(point[2], point[1])
        νr = true
        νθ = true

        if rs < Krang.horizon(pixel.metric)
            continue
        end

        intersection = Intersection(ts, rs, θs, ϕs, νr, νθ)
        intersections += !didintersect ? z : material(pixel, intersection)#didintersect ? 1 : 0
        origin = line_point_2
    end
    return intersections
end

function line_intersection(origin::NTuple{3,T}, line_point_2, geometry::AbstractLevelSetGeometry{T}) where T
    didintersect = geometry(origin...)*geometry(line_point_2...) < zero(T)
    if didintersect
        direction = (line_point_2[1] - origin[1], line_point_2[2] - origin[2], line_point_2[3] - origin[3])

        t = find_zero((x)->geometry((origin .+ (direction .* x))...),  (zero(T), one(T)))

        return true, origin .+ (direction .* t)
    end
    return false, (zero(T), zero(T), zero(T))

end