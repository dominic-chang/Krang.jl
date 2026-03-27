@static if Sys.isapple()
    try
        @eval using Metal
        metal_loaded = true
    catch _
        metal_loaded = false
    end
end
Mat = Sys.isapple() ? Metal.MtlArray : Array

met = Krang.Kerr(0.5f0)

α_s = Float32[i for i in -10:0.01:10] 
β_s = Float32[i for i in -10:0.01:10]
θo_s = Float32[i*π/180 for i in range(1, 40, length(α_s))]
θs_s = Float32[i*π/180 for i in range(80, 100, length(α_s))]
α_mtl = Mat(α_s)
β_mtl = Mat(β_s)
θo_mtl = Mat(θo_s)
θs_mtl = Mat(θs_s)

λ_s = Krang.λ.(Ref(met), α_s, θo_s)
λ_mtl = Krang.λ.(Ref(met), α_mtl, θo_mtl)
@test Array(λ_mtl) ≈ λ_s

η_s = Krang.η.(Ref(met), α_s, β_s, θo_s)
η_mtl = Krang.η.(Ref(met), α_mtl, β_mtl, θo_mtl)
@test Array(η_mtl) ≈ η_s

radial_roots_s = Krang.get_radial_roots.(Ref(met), η_s, λ_s)
radial_roots_mtl = Krang.get_radial_roots.(Ref(met), η_mtl, λ_mtl)
@test maximum(abs.(sum.(collect.(Array(radial_roots_mtl)) .- collect.(radial_roots_s)))) ≈ 0f0 atol=1e-5

numreals = sum.(Ref(Krang._isreal2), radial_roots_s)
numreals_mtl = sum.(Ref(Krang._isreal2), radial_roots_mtl)
@test Array(numreals_mtl) == numreals

Ir_inf_s = Krang.Ir_inf.(Ref(met), map(x->real.(x),radial_roots_s))
Ir_inf_mtl = Krang.Ir_inf.(Ref(met), map(x->real.(x),radial_roots_mtl))
@test Array(Ir_inf_mtl) ≈ Ir_inf_s

τ_total_s = Krang.total_mino_time.(Ref(met), radial_roots_s)
τ_total_mtl = Krang.total_mino_time.(Ref(met), radial_roots_mtl)
@test maximum(abs.((Array(τ_total_mtl) ./ τ_total_s) .- 1f0)) ≈ 0f0 atol=1e-2

Gθo_Gθhat_s = Krang._absGθo_Gθhat.(Ref(met), θo_s, η_s, λ_s)
Gθo_Gθhat_mtl = Krang._absGθo_Gθhat.(Ref(met), θo_mtl, η_mtl, λ_mtl)
@test maximum(abs.(sum.(collect.(Array(Gθo_Gθhat_mtl)) .- collect.(Gθo_Gθhat_s)))) ≈ 0f0 atol=1e-4

pix_s = Krang.IntensityPixel.(Ref(met), α_s, β_s, θo_s)
pix_mtl = Krang.IntensityPixel.(Ref(met), α_mtl, β_mtl, θo_mtl)

pix_s = Krang.SlowLightIntensityPixel.(Ref(met), α_s, β_s, θo_s)
pix_mtl = Krang.SlowLightIntensityPixel.(Ref(met), α_mtl, β_mtl, θo_mtl)

@test pix_s ≈ pix_mtl


