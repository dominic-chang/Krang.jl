export load, translate, scale, rotate, MeshGeometry

@doc """
    Load a mesh from a file.
    
    # Arguments
    - `filename::String`: The path to the file containing the mesh.
    
    # Returns
    - A `Mesh` object representing the mesh.
"""
struct MeshGeometry <: AbstractGeometry
    geometryBasicsMesh::GeometryBasics.Mesh
end

function translate(mesh::MeshGeometry, x, y, z)
    points = (Ref([x, y, z]) .+ mesh.geometryBasicsMesh.position)
    faces = getfield(mesh.geometryBasicsMesh, :faces)
    return MeshGeometry(GeometryBasics.Mesh([GeometryBasics.Point(x) for x in points], faces))
end

function scale(mesh::MeshGeometry, multiple)
    points = (multiple .* mesh.geometryBasicsMesh.position)
    faces = getfield(mesh.geometryBasicsMesh, :faces)
    return MeshGeometry(GeometryBasics.Mesh([GeometryBasics.Point(x) for x in points], faces))
end

function rotate(mesh::MeshGeometry, angle, x, y, z)
    angleaxis = Rotations.AngleAxis(angle, x, y, z)
    points = Ref(angleaxis) .* (mesh.geometryBasicsMesh.position)
    faces = getfield(mesh.geometryBasicsMesh, :faces)
    return MeshGeometry(GeometryBasics.Mesh([GeometryBasics.Point(x) for x in points], faces))
end

function raytrace(camera::AbstractCamera, mesh::Krang.Mesh{<:MeshGeometry, <:AbstractMaterial}; res = 100)
    mesh_geometry = mesh.geometry.geometryBasicsMesh
    material = mesh.material
    faces = begin
        temp = getfield(mesh_geometry, :faces)
        len = length(temp)
        reshape([i[j] for i in temp for j = 1:3], 3, len)
    end
    vertices = begin
        temp = mesh_geometry.position
        len = length(temp)
        reshape([(i[j]) for i in temp for j = 1:3], 3, len)
    end
    intersections = zeros(Float64, size(camera.screen.pixels))
    Threads.@threads for I in CartesianIndices(camera.screen.pixels)
        intersections[I] = _raytrace(camera.screen.pixels[I], faces, vertices, material; res)
    end
    return intersections
end

function _raytrace(
    pixel::AbstractPixel{T},
    faces::Matrix{GeometryBasics.OffsetInteger{-1,UInt32}},
    vertices::Matrix{T},
    material::AbstractMaterial;
    res = 100,
) where {T}
    metric = pixel.metric
    intersections = 0
    ray = Vector{Intersection{T}}(undef, res)
    generate_ray!(ray, pixel, res)
    (; rs, θs, ϕs, νr, νθ) = ray[1]
    origin = boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(metric, rs, θs, ϕs)#(rs * sin(θs)*cos(ϕs), rs * sin(θs)*sin(ϕs), rs * cos(θs))

    for i = 2:res
        (; rs, θs, ϕs) = ray[i]
        if rs <= Krang.horizon(pixel.metric) || iszero(rs)
            continue
        end
        line_point_2 =
            boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(metric, rs, θs, ϕs)
        for j = 1:(size(faces)[2])
            f1, f2, f3 = @view faces[:, j]
            v1 = @view vertices[:, f1]
            v2 = @view vertices[:, f2]
            v3 = @view vertices[:, f3]
            didintersect, point = line_intersection(origin, line_point_2, v1, v2, v3)
            # Take νr and νθ from the first intersection 
            # TODO: Come up with a better way to handle this
            rnew = sqrt(sum(point .^ 2))
            θnew = acos(point[3] / rnew)
            ϕnew = atan(point[2], point[1])
            intersections += didintersect ? material(pixel, Intersection(0.0, rnew, θnew, ϕnew, νr, νθ)) : 0
        end
        origin = line_point_2
    end
    return intersections
end
