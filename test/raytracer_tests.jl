@testset "Raytracer Functions" begin
    @testset "Boyer-Lindquist to Kerr-Schild transformations" begin
        met = Krang.Kerr(0.99)
        rs = 1e10
        θs = π/4.0
        ϕs = π/4.0

        x = rs * sin(θs) * cos(ϕs)
        y = rs * sin(θs) * sin(ϕs)
        z = rs * cos(θs)

        # Boyer-Lindquist and Kerr-Schild coordinates are the same far from the black hole
        trat, xrat, yrat, zrat =  Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild(met, 1e10, rs, θs, ϕs) ./ (1e10, x, y, z) 
        @test trat ≈ 1.0 atol = 1e-5 
        @test xrat ≈ 1.0 atol = 1e-5 
        @test yrat ≈ 1.0 atol = 1e-5 
        @test zrat ≈ 1.0 atol = 1e-5 

        # Boyer-Lindquist and Kerr-Schild coordinates are the same far from the black hole
        xrat, yrat, zrat =  Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(met, rs, θs, ϕs) ./ (x, y, z) 
        @test trat ≈ 1.0 atol = 1e-5 
        @test xrat ≈ 1.0 atol = 1e-5 
        @test yrat ≈ 1.0 atol = 1e-5 
        @test zrat ≈ 1.0 atol = 1e-5 

    end

    @testset "Level Set" begin
        metric = Krang.Kerr(0.999)
        θo = 89/180*π 
        ρmax = 10.0
        sze = 200 
            
        camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze)

        struct Cone{T} <: Krang.AbstractLevelSetGeometry{T}
            θo::T
        end
        
        function (cone::Cone)(x,y,z)
            r = hypot(x,y,z)
            return z/r - cos(cone.θo)
        end
        struct XMaterial{N} <: Krang.AbstractMaterial 
            subimgs::NTuple{N,Int}
        end
        Krang.horizon(metric)
        function (mat::XMaterial)(pix, intersection) where T
            (;rs,) = intersection
            return horizon(pix.metric)<rs <10.0
        end

        mesh1 = Krang.Mesh(Cone(π/4), XMaterial((0,1,2,3)))
        intersections1 = raytrace(camera, mesh1, res=1_000)

        mesh2 = Krang.Mesh(Krang.ConeGeometry(π/4), XMaterial((0,1,2,3)))
        intersections2 = raytrace(camera, mesh2)

        sum(intersections1 .- intersections2) .≈ 0.0
    end
end