const metal_loaded = if Sys.isapple()
    try
        @eval using Metal
        true
    catch _
        false
    end
else
    false
end

@testset "Metal" begin
    if !metal_loaded
        @test true
    else
        # Mostly to check if Metal compilation happens.
        Mat = Metal.MtlArray

        met = Krang.Kerr(0.5f0)

        α_s = Float32[i for i in -10:0.01:10]
        β_s = Float32[i for i in -10:0.01:10]
        θo_s = Float32[i * π / 180 for i in range(1, 40, length(α_s))]
        θs_s = Float32[i * π / 180 for i in range(80, 100, length(α_s))]
        α_mtl = Mat(α_s)
        β_mtl = Mat(β_s)
        θo_mtl = Mat(θo_s)
        θs_mtl = Mat(θs_s)

        λ_s = Krang.λ.(Ref(met), α_s, θo_s)
        λ_mtl = Krang.λ.(Ref(met), α_mtl, θo_mtl)

        η_s = Krang.η.(Ref(met), α_s, β_s, θo_s)
        η_mtl = Krang.η.(Ref(met), α_mtl, β_mtl, θo_mtl)

        radial_roots_s = Krang.get_radial_roots.(Ref(met), η_s, λ_s)
        radial_roots_mtl = Krang.get_radial_roots.(Ref(met), η_mtl, λ_mtl)
        @test maximum(abs.(sum.(collect.(Array(radial_roots_mtl)) .- collect.(radial_roots_s)))) ≈ 0f0 atol = 1e-5

        numreals = sum.(Ref(Krang._isreal2), radial_roots_s)
        numreals_mtl = sum.(Ref(Krang._isreal2), radial_roots_mtl)

        Ir_inf_s = Krang.Ir_inf.(Ref(met), map(x -> real.(x), radial_roots_s))
        Ir_inf_mtl = Krang.Ir_inf.(Ref(met), map(x -> real.(x), radial_roots_mtl))

        τ_total_s = Krang.total_mino_time.(Ref(met), radial_roots_s)
        τ_total_mtl = Krang.total_mino_time.(Ref(met), radial_roots_mtl)

        Gθo_Gθhat_s = Krang._absGθo_Gθhat.(Ref(met), θo_s, η_s, λ_s)
        Gθo_Gθhat_mtl = Krang._absGθo_Gθhat.(Ref(met), θo_mtl, η_mtl, λ_mtl)

        pix_s = Krang.IntensityPixel.(Ref(met), α_s, β_s, θo_s)
        pix_mtl = Krang.IntensityPixel.(Ref(met), α_mtl, β_mtl, θo_mtl)

        Ivals_s = Krang.radial_inf_integrals.(Ref(met), radial_roots_s)
        Ivals_mtl = Krang.radial_inf_integrals.(Ref(met), radial_roots_mtl)

        Iϕ_s = Krang.Iϕ_inf.(Ref(met), radial_roots_s, λ_s)
        Iϕ_mtl = Krang.Iϕ_inf.(Ref(met), radial_roots_mtl, λ_mtl)

        It_s = Krang.It_inf.(Ref(met), radial_roots_s, λ_s)
        It_mtl = Krang.It_inf.(Ref(met), radial_roots_mtl, λ_mtl)

        Gϕo_Gϕhat_s = Krang._absGϕo_Gϕhat.(Ref(met), θo_s, η_s, λ_s)
        Gϕo_Gϕhat_mtl = Krang._absGϕo_Gϕhat.(Ref(met), θo_mtl, η_mtl, λ_mtl)

        Gto_Gthat_s = Krang._absGto_Gthat.(Ref(met), θo_s, η_s, λ_s)
        Gto_Gthat_mtl = Krang._absGto_Gthat.(Ref(met), θo_mtl, η_mtl, λ_mtl)

        pix_s = Krang.IntensityPixel.(Ref(met), α_s, β_s, θo_s)
        pix_mtl = Krang.IntensityPixel.(Ref(met), α_mtl, β_mtl, θo_mtl)

        struct MetalConeMaterial{N,T} <: Krang.AbstractMaterial
            subimgs::NTuple{N,Int}
            weight::T
        end

        Krang.isFastLight(::MetalConeMaterial) = true
        Krang.isAxisymmetric(::MetalConeMaterial) = true

        function (mat::MetalConeMaterial)(pix, intersection)
            α, β = Krang.screen_coordinate(pix)
            return mat.weight * (intersection.rs + intersection.θs + α - β)
        end

        mesh = Krang.Mesh(Krang.ConeGeometry(π / 4f0), MetalConeMaterial((0, 1), 0.125f0))
        cone_obs_s = raytrace.(pix_s, Ref(mesh))
        cone_obs_mtl = raytrace.(pix_mtl, Ref(mesh))
        cone_delta = abs.(Array(cone_obs_mtl) .- cone_obs_s)

        @test maximum(cone_delta) < 0.04f0
        @test sum(cone_delta) / length(cone_delta) < 0.003f0
        @test sum(abs, Array(cone_obs_mtl)) > 0f0
    end
end
