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
    return GeometryBasics.Mesh([Point(x) for x in points], faces)
end

function scale(mesh::MeshGeometry, multiple)
    points = (multiple .* mesh.position)
    faces = getfield(getfield(mesh, :simplices), :faces)
    return GeometryBasics.Mesh([Point(x) for x in points], faces)
end

function rotate(mesh::MeshGeometry, angle, x, y, z)
    angleaxis = Rotations.AngleAxis(angle, x, y, z)
    points = Ref(angleaxis) .* (mesh.position)
    faces = getfield(getfield(mesh, :simplices), :faces)
    return GeometryBasics.Mesh([Point(x) for x in points], faces)
end