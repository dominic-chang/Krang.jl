@testset "Raytracer Functions" begin
    @testset "Generate ray" begin
        for T in (Float32, Float64)
            pixel = Krang.IntensityPixel(Krang.Kerr(T(0.5)), T(4), T(4), T(0.5))
            ray = Krang.generate_ray(pixel, 2)

            @test ray isa Vector{Krang.Intersection{T}}
            @test length(ray) == 2
        end
    end

    @testset "Threaded rendering" begin
        struct ThreadedTestMaterial <: Krang.AbstractMaterial
            weight::Float64
        end
        (material::ThreadedTestMaterial)(pixel, geometry::Krang.ConeGeometry) =
            material.weight

        camera = Krang.IntensityCamera(Krang.Kerr(0.5), 0.5, -1.0, 1.0, -1.0, 1.0, 2)
        scene = (
            Krang.Mesh(Krang.ConeGeometry(0.2), ThreadedTestMaterial(1.0)),
            Krang.Mesh(Krang.ConeGeometry(0.4), ThreadedTestMaterial(2.0)),
        )
        store = zeros(2, 2)

        rendered = Krang.render_cpu_threaded!(store, camera, scene)

        @test rendered === store
        @test store == fill(3.0, 2, 2)
    end

    @testset "Boyer-Lindquist to Kerr-Schild transformations" begin
        met = Krang.Kerr(0.99)
        rs = 1e10
        θs = π / 4.0
        ϕs = π / 4.0

        x = rs * sin(θs) * cos(ϕs)
        y = rs * sin(θs) * sin(ϕs)
        z = rs * cos(θs)

        # Boyer-Lindquist and Kerr-Schild coordinates are the same far from the black hole
        trat, xrat, yrat, zrat =
            Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild(met, 1e10, rs, θs, ϕs) ./
            (1e10, x, y, z)
        @test trat ≈ 1.0 atol = 1e-5
        @test xrat ≈ 1.0 atol = 1e-5
        @test yrat ≈ 1.0 atol = 1e-5
        @test zrat ≈ 1.0 atol = 1e-5

        # Boyer-Lindquist and Kerr-Schild coordinates are the same far from the black hole
        xrat, yrat, zrat =
            Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
                met,
                rs,
                θs,
                ϕs,
            ) ./ (x, y, z)
        @test trat ≈ 1.0 atol = 1e-5
        @test xrat ≈ 1.0 atol = 1e-5
        @test yrat ≈ 1.0 atol = 1e-5
        @test zrat ≈ 1.0 atol = 1e-5

        # ϕ_kerr_schild / ϕ_BL must return T (not Union{T, Float64}) on Float32.
        let metf = Krang.Kerr(0.5f0)
            @test (@inferred Float32 Krang.ϕ_kerr_schild(metf, 5.0f0, 0.0f0)) isa Float32
            @test (@inferred Float32 Krang.ϕ_BL(metf, 5.0f0, 0.0f0)) isa Float32
        end
    end

    @testset "Level Set" begin
        metric = Krang.Kerr(0.999)
        θo = 89 / 180 * π
        ρmax = 10.0
        sze = 200

        camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze)

        struct Cone{T} <: Krang.AbstractLevelSetGeometry{T}
            θo::T
        end

        function (cone::Cone)(x, y, z)
            r = hypot(x, y, z)
            return z / r - cos(cone.θo)
        end
        struct XMaterial{N} <: Krang.AbstractMaterial
            subimgs::NTuple{N,Int}
        end
        Krang.horizon(metric)
        function (mat::XMaterial)(pix, intersection)
            (; rs,) = intersection
            return horizon(pix.metric) < rs < 10.0
        end

        mesh1 = Krang.Mesh(Cone(π / 4), XMaterial((0, 1, 2, 3)))
        intersections1 = raytrace(camera, mesh1, res = 1_000)

        mesh2 = Krang.Mesh(Krang.ConeGeometry(π / 4), XMaterial((0, 1, 2, 3)))
        intersections2 = raytrace(camera, mesh2)

        sum(intersections1 .- intersections2) .≈ 0.0
    end
end
