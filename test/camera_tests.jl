@testset "Camera" begin
    met = Krang.Kerr(0.0)
    α = √27*cos(π/4)
    β = √27*sin(π/4)
    θo = π/4
    crit_roots = (-6.0 + 0im, 0.0 + 0im, 3.0 + 0im, 3.0 + 0im)
    @testset "Abstract Camera" begin

        @testset "Abstract Pixel" begin
            struct TestPixel{T} <: Krang.AbstractPixel{T}
                metric::Kerr{T}
                screen_coordinate::NTuple{2, T}
                θo::T
                function TestPixel(met::Kerr{T}, α, β, θo) where {T}
                    new{T}(
                        met,
                        (α, β), 
                        θo
                    )
                end
            end
            abspix = TestPixel(met, α, β, θo)

            @test Krang.screen_coordinate(abspix) == (α, β)
            @test Krang.inclination(abspix) == θo
            @test Krang.η(abspix) == Krang.η(met, α, β, θo)
            @test Krang.λ(abspix) == Krang.λ(met, α, θo)
            @test max((abs.(Krang.roots(abspix) .- crit_roots))...) ≈ 0.0 atol = 1e-5
            @test Krang.I0_inf(abspix) ≈ Krang.Ir_inf(met, crit_roots)
            @test Krang.Ir_inf(abspix) ≈ Krang.Ir_inf(met, crit_roots)
        end
    end

    @testset "Intensity Camera" begin
        intcamera = Krang.IntensityCamera(met, -√27, √27, -√27, √27, θo, 10)
        intscreen = intcamera.screen

        @testset "Intensity Pixel" begin
            intpix = Krang.IntensityPixel(met, α, β, θo)

            @test Krang.screen_coordinate(intpix) == (α, β)
            @test Krang.inclination(intpix) == θo
            @test Krang.η(intpix) == Krang.η(met, α, β, θo)
            @test Krang.λ(intpix) == Krang.λ(met, α, θo)
            @test max((abs.(Krang.roots(intpix) .- crit_roots))...) ≈ 0.0 atol = 1e-5
            @test Krang.I0_inf(intpix) ≈ Krang.Ir_inf(met, crit_roots)
            @test Krang.Ir_inf(intpix) ≈ Krang.Ir_inf(met, crit_roots)
        end

        @testset "Intensity Screen" begin
            intscreen = Krang.IntensityScreen(met, -√27, √27, -√27, √27, θo, 10)
        end

    end
end