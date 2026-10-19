@testset "Mesh Geometry Manipulation" begin
    struct GeometryTestMaterial <: Krang.AbstractMaterial end

    disk = Krang.Disk()
    attributed_disk = Krang.Disk(attributes = (name = :disk,))
    mesh = Krang.Mesh(disk, GeometryTestMaterial())
    append_scene = getfield(Krang, Symbol("⊕"))

    @test disk.attributes === nothing
    @test attributed_disk.attributes == (name = :disk,)
    @test append_scene(Krang.Scene(), mesh) == Krang.add(Krang.Scene(), mesh)
    @test :load ∉ names(Krang)

    bunny = FileIO.load(
        Downloads.download(
            "https://graphics.stanford.edu/~mdfisher/Data/Meshes/bunny.obj",
            "bunny.obj",
        ),
    )
    rot = Rotations.AngleAxis(π / 2, 0.0, 0.0, 1.0)
    points = (Ref(rot) .* (bunny.position .* 150)) .+ Ref([15, 3, 0])
    faces = getfield(bunny, :faces)
    bunny_mesh1 = GeometryBasics.Mesh([Point(x) for x in points], faces)
    bunny_mesh_geometry = MeshGeometry(bunny)
    bunny_mesh =
        translate(rotate(scale(bunny_mesh_geometry, 150), π / 2, 0.0, 0.0, 1.0), 15, 3, 0)
    @test all(
        map(
            x -> x[1] == x[2],
            zip(bunny_mesh.geometryBasicsMesh.position, bunny_mesh1.position),
        ),
    )

    @testset "Triangulated cone matches ConeGeometry" begin
        "Constructs cone mesh with apex at the origin, out to r = rmax."
        function cone_mesh(half_angle, rmax; nϕ = 64, nr = 8)
            sθ, cθ = sincos(half_angle)
            pts = GeometryBasics.Point{3,Float64}[GeometryBasics.Point(0.0, 0.0, 0.0)]
            for i = 1:nr, j = 1:nϕ
                r = rmax * i / nr
                ϕ = 2π * (j - 1) / nϕ
                push!(pts, GeometryBasics.Point(r * sθ * cos(ϕ), r * sθ * sin(ϕ), r * cθ))
            end
            idx(i, j) = 1 + (i - 1) * nϕ + mod1(j, nϕ)
            faces = GeometryBasics.GLTriangleFace[]
            for j = 1:nϕ # fan of triangles around the apex
                push!(faces, GeometryBasics.GLTriangleFace(1, idx(1, j), idx(1, j + 1)))
            end
            for i = 1:(nr-1), j = 1:nϕ # quad bands, each split into two triangles
                push!(
                    faces,
                    GeometryBasics.GLTriangleFace(
                        idx(i, j),
                        idx(i + 1, j),
                        idx(i + 1, j + 1),
                    ),
                )
                push!(
                    faces,
                    GeometryBasics.GLTriangleFace(
                        idx(i, j),
                        idx(i + 1, j + 1),
                        idx(i, j + 1),
                    ),
                )
            end
            return GeometryBasics.Mesh(pts, faces)
        end

        struct ConeHitMaterial{N} <: Krang.AbstractMaterial
            subimgs::NTuple{N,Int}
        end
        function (mat::ConeHitMaterial)(pix, intersection)
            (; rs,) = intersection
            return Krang.horizon(pix.metric) < rs < 10.0
        end

        metric = Krang.Kerr(0.999)
        θo = 89 / 180 * π
        ρmax = 10.0
        sze = 30
        camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze)

        # The mesh runs out to r = 12 for a little extra buffer
        cone = MeshGeometry(cone_mesh(π / 4, 12.0))
        from_mesh =
            raytrace(camera, Krang.Mesh(cone, ConeHitMaterial((0, 1, 2, 3))), res = 400)
        from_cone = raytrace(
            camera,
            Krang.Mesh(Krang.ConeGeometry(π / 4), ConeHitMaterial((0, 1, 2, 3))),
        )

        # `ConeGeometry` only computes the sub-images listed in `subimgs`, while the mesh
        # marcher finds every crossing along the ray, so a clamp is used. 
        nsub = 4
        delta = min.(from_mesh, nsub) .- min.(from_cone, nsub)
        @test sum(from_cone) > 0 # the cone is actually being hit
        # Counting mismatches rather than summing them matters: `sum(delta)` lets a ray that
        # gains a crossing cancel one that loses it. The bound still catches a wrong surface:
        # opening the cone by 2% puts 4 pixels out, stopping the mesh at r = 6 puts 164 out.
        @test count(!iszero, delta) == 0
        # Above the clamp flat facets may legitimately cost a crossing, but not more than two.
        @test maximum(abs, from_mesh .- from_cone) <= 2
    end
end
