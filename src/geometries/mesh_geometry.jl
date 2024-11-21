export load, translate, scale, rotate, MeshGeometry

@doc """
    Load a mesh from a file.
    
    # Arguments
    - `filename::String`: The path to the file containing the mesh.
    
    # Returns
    - A `Mesh` object representing the mesh.
"""
const MeshGeometry = GeometryBasics.Mesh

function translate(mesh::MeshGeometry, x, y, z)
    points = (Ref([x, y, z]) .+ mesh.position)
    faces = getfield(getfield(mesh, :simplices), :faces)
    return GeometryBasics.Mesh([GeometryBasics.Point(x) for x in points], faces)
end

function scale(mesh::MeshGeometry, multiple)
    points = (multiple .* mesh.position)
    faces = getfield(getfield(mesh, :simplices), :faces)
    return GeometryBasics.Mesh([GeometryBasics.Point(x) for x in points], faces)
end

function rotate(mesh::MeshGeometry, angle, x, y, z)
    angleaxis = Rotations.AngleAxis(angle, x, y, z)
    points = Ref(angleaxis) .* (mesh.position)
    faces = getfield(getfield(mesh, :simplices), :faces)
    return GeometryBasics.Mesh([GeometryBasics.Point(x) for x in points], faces)
end

function raytrace(camera::AbstractCamera, mesh_geometry::MeshGeometry; res=100) 
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

function raytrace(pixel::AbstractPixel{T}, faces::Matrix{GeometryBasics.OffsetInteger{-1, UInt32}}, vertices::Matrix{T}; res=100) where T
    intersections = 0
    ray = Vector{Intersection{T}}(undef,res)
    generate_ray!(ray, pixel, res)
    (;rs, θs, ϕs) = ray[1]
    origin = (rs * sin(θs)*cos(ϕs), rs * sin(θs)*sin(ϕs), rs * cos(θs))

    for i in 2:res
        (;rs, θs, ϕs) = ray[i]
        line_point_2 = (rs * sin(θs)*cos(ϕs), rs * sin(θs)*sin(ϕs), rs * cos(θs))
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