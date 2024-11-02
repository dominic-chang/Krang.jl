"""
    Level Sets should be a functor
"""
abstract type AbstractLevelSetGeometry{T} <: AbstractGeometry end

function raytrace(camera, mesh::Mesh{<:AbstractLevelSetGeometry, <:AbstractMaterial}; res=100) 
    intersections = zeros(Float64, size(camera.screen.pixels))
    Threads.@threads for I in CartesianIndices(camera.screen.pixels)
        intersections[I] = raytrace(camera.screen.pixels[I], mesh;res)
    end
    return intersections
end

function raytrace(pixel::AbstractPixel{T}, mesh::Mesh{<:AbstractLevelSetGeometry, <:AbstractMaterial}; res=100) where T
    geometry = mesh.geometry
    intersections = 0
    ray = zeros(T, 3, res)
    generate_ray!(ray, pixel, res)
    origin = @view ray[:,1]
    for i in 1:res
        line_point_2 = @view ray[:,i]
        if isinf(line_point_2[1]) || isinf(line_point_2[2]) || isinf(line_point_2[3])
            continue
        end
        didintersect, point = line_intersection(origin, line_point_2, geometry)
        intersections += didintersect ? 0 : mesh.material(point...)#didintersect ? 1 : 0
        origin = line_point_2
    end
    return intersections
end

function line_intersection(origin::AbstractVector{T}, line_point_2, geometry::AbstractLevelSetGeometry{T}) where T
    didintersect = geometry(origin...)*geometry(line_point_2...) < zero(T)
    if didintersect
        direction = (line_point_2[1] - origin[1], line_point_2[2] - origin[2], line_point_2[3] - origin[3])

        t = find_zero((x)->geometry((origin .+ (direction .* x))...),  (zero(T), one(T)))

        return true, origin .+ (direction .* t)
    end
    return false, (zero(T), zero(T), zero(T))

end