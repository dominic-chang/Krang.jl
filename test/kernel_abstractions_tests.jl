@static if Sys.isapple()
    try
        @eval using Metal
        metal_loaded = true
    catch _
        metal_loaded = false
    end
end

@testset "Kernel Abstractions" begin
    Base.retry_load_extensions()
    @test !isnothing(Base.get_extension(Krang, :KrangKernelAbstractionsExt))

    met = Krang.Kerr(0.3)
    θo = π / 3
    αmin, αmax = -2.0, 2.0
    βmin, βmax = -1.5, 1.5
    screen_res = 4
    Mat = Sys.isapple() ? Metal.MtlMatrix : Matrix

    function max_root_delta(lhs, rhs)
        return maximum(abs.(collect(lhs) .- collect(rhs)))
    end

    @testset "IntensityScreen constructor" begin
        cpu_screen = Krang.IntensityScreen(met, αmin, αmax, βmin, βmax, θo, screen_res)
        ka_screen = Krang.IntensityScreen(met, αmin, αmax, βmin, βmax, θo, screen_res, Mat)

        @test size(ka_screen.pixels) == size(cpu_screen.pixels)
        @test ka_screen.αrange == cpu_screen.αrange
        @test ka_screen.βrange == cpu_screen.βrange

        for idx in CartesianIndices(cpu_screen.pixels)
            cpu_pixel = cpu_screen.pixels[idx]
            ka_pixel = ka_screen.pixels[idx]
            @test Krang.screen_coordinate(ka_pixel)[1] ≈ Krang.screen_coordinate(cpu_pixel)[1]
            @test Krang.screen_coordinate(ka_pixel)[2] ≈ Krang.screen_coordinate(cpu_pixel)[2]
            @test Krang.inclination(ka_pixel) == Krang.inclination(cpu_pixel)
            @test Krang.η(ka_pixel) ≈ Krang.η(cpu_pixel)
            @test Krang.λ(ka_pixel) ≈ Krang.λ(cpu_pixel)
            @test max_root_delta(Krang.roots(ka_pixel), Krang.roots(cpu_pixel)) ≤ 1e-10
            @test Krang.I0_inf(ka_pixel) ≈ Krang.I0_inf(cpu_pixel)
            @test Krang.total_mino_time(ka_pixel) ≈ Krang.total_mino_time(cpu_pixel)
        end
    end

    @testset "generate_rays" begin
        screen = Krang.IntensityScreen(met, αmin, αmax, βmin, βmax, θo, 2)
        ray_res = 5
        rays = Krang.generate_rays(screen.pixels, ray_res)

        @test size(rays) == (2, 2, ray_res)

        for idx in CartesianIndices(screen.pixels)
            cpu_ray = Krang.generate_ray(screen.pixels[idx], ray_res)
            for k = 1:ray_res
                ka_intersection = rays[idx[1], idx[2], k]
                cpu_intersection = cpu_ray[k]
                @test ka_intersection.ts ≈ cpu_intersection.ts
                @test ka_intersection.rs ≈ cpu_intersection.rs
                @test ka_intersection.θs ≈ cpu_intersection.θs
                @test ka_intersection.ϕs ≈ cpu_intersection.ϕs
                @test ka_intersection.νr == cpu_intersection.νr
                @test ka_intersection.νθ == cpu_intersection.νθ
            end
        end
    end

    @testset "render!" begin
        struct TestKernelMaterial <: Krang.AbstractMaterial
            weight::Float64
        end

        function (mat::TestKernelMaterial)(pix, geometry::Krang.ConeGeometry)
            α, β = Krang.screen_coordinate(pix)
            return mat.weight * (α - 2β + geometry.opening_angle)
        end

        camera = Krang.IntensityCamera(met, θo, αmin, αmax, βmin, βmax, 3)
        scene = (
            Krang.Mesh(Krang.ConeGeometry(0.2), TestKernelMaterial(1.0)),
            Krang.Mesh(Krang.ConeGeometry(0.4), TestKernelMaterial(-0.5)),
        )
        store = zeros(Float64, size(camera.screen.pixels))
        rendered = Krang.render!(store, camera, scene)

        expected = similar(store)
        for idx in CartesianIndices(expected)
            pix = camera.screen.pixels[idx]
            expected[idx] = sum(mesh.material(pix, mesh.geometry) for mesh in scene)
        end

        @test rendered == expected
        @test store == expected
    end
end
