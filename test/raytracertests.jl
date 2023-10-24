@testset "Raytracer Functions" begin

    @testset "Radial Emission" begin
        a = 0.99
        met = Kang.Kerr(a)

        @test isnan(emission_radius(met, 10, 1.0, π/2, π/2, true, 2)[1])

        @testset "Case 2" begin
            α = 10.0
            β = 10.0
            θo = π / 4

            ηcase1 = η(met, α, β, θo)
            λcase1 = λ(met, α, θo)
            roots = get_radial_roots(met, ηcase1, λcase1)
            _, _, _, root = roots
            @test sum(Kang._isreal2.(roots)) == 4
            rs = 1.1 * real(root)
            τ1 = Ir(met, true, rs, ηcase1, λcase1)[1]
            @test Kang.emission_radius(met, α, β, τ1, θo)[1] / rs ≈ 1 atol = 1e-5
        end

        @testset "Case 3" begin
            α = 1.0
            β = 1.0
            θo = π / 4

            ηcase3 = η(met, α, β, θo)
            λcase3 = λ(met, α, θo)
            roots = get_radial_roots(met, ηcase3, λcase3)
            @test sum(Kang._isreal2.(roots)) == 2
            rs = 1.1horizon(met)
            τ3 = Ir(met, true, rs, ηcase3, λcase3)[1]
            @test Kang.emission_radius(met, α, β, τ3, θo)[1] / rs ≈ 1 atol = 1e-5
        end
        @testset "Case 4" begin
            α = 0.1
            β = 0.1
            θo = π / 4

            ηcase4 = η(met, α, β, θo)
            λcase4 = λ(met, α, θo)
            roots = get_radial_roots(met, ηcase4, λcase4)
            @test sum(Kang._isreal2.(roots)) == 0
            rs = 1.1horizon(met)
            τ4 = Ir(met, true, rs, ηcase4, λcase4)[1]
            @test Kang.emission_radius(met, α, β, τ4, θo)[1] / rs ≈ 1 atol = 1e-5
        end
    end
end