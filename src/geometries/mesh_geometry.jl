export load, translate, scale, rotate

"""
    $TYPEDEF

Mesh Geometry Type. Embeds a mesh in 3D space assuming Boyer-Lindquist coordinates.
"""
const MeshGeometry = GeometryBasics.Mesh

function translate(mesh, x, y, z)
    points = (Ref([x, y, z]) .+ mesh.position)
    faces = getfield(getfield(mesh, :simplices), :faces)
    return GeometryBasics.Mesh([Point(x) for x in points], faces)
end

function scale(mesh, multiple)
    points = (multiple .* mesh.position)
    faces = getfield(getfield(mesh, :simplices), :faces)
    return GeometryBasics.Mesh([Point(x) for x in points], faces)
end

function rotate(mesh, angle, x, y, z)
    angleaxis = Rotations.AngleAxis(angle, x, y, z)
    points = Ref(angleaxis) .* (mesh.position)
    faces = getfield(getfield(mesh, :simplices), :faces)
    return GeometryBasics.Mesh([Point(x) for x in points], faces)
end