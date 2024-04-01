@testset "Camera" begin
    met = Krang.Kerr(0.0)
    α = √27*cos(π/4)
    β = √27*sin(π/4)
    θo = π/4
    crit_roots = (-6.0 + 0im, 0.0 + 0im, 3.0 + 0im, 3.0 + 0im)
    @testset "Abstract Pixels API" begin
        struct TestPixel{T} <: Krang.AbstractPixel 
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
        intpix = Krang.IntensityPixel(met, α, β, θo)

        @test Krang.screen_coordinate(intpix) == (α, β)
        @test Krang.inclination(intpix) == θo
        @test Krang.η(intpix) == Krang.η(met, α, β, θo)
        @test Krang.λ(intpix) == Krang.λ(met, α, θo)
        @test max((abs.(Krang.roots(intpix) .- crit_roots))...) ≈ 0.0 = 1e-5
        @test isnan(Krang.I0_inf(intpix)) ≈ isnan(Krang.Ir_inf(met, crit_roots))
        @test isnan(Krang.Ir_inf(intpix)) ≈ isnan(Krang.Ir_inf(met, crit_roots))
    end
end