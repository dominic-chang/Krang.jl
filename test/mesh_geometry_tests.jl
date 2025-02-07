@testset "Mesh Geometry Manipulation" begin
    bunny = FileIO.load(
        download(
            "https://graphics.stanford.edu/~mdfisher/Data/Meshes/bunny.obj",
            "bunny.obj",
        ),
    )
    rot = Rotations.AngleAxis(π / 2, 0.0, 0.0, 1.0)
    points = (Ref(rot) .* (bunny.position .* 150)) .+ Ref([15, 3, 0])
    faces = getfield(bunny, :faces)
    bunny_mesh1 = GeometryBasics.Mesh([Point(x) for x in points], faces)

    bunny_mesh = translate(rotate(scale(bunny, 150), π / 2, 0.0, 0.0, 1.0), 15, 3, 0)
    @test all(map(x -> x[1] == x[2], zip(bunny_mesh.position, bunny_mesh1.position)))
end
