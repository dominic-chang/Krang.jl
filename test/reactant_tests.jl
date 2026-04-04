push!(LOAD_PATH, dirname(@__DIR__))
using Krang
using Test
using Reactant 
using StaticArrays
import Reactant: to_rarray

same_or_approx(x, y; atol = 1e-12) = isequal(x, y) || (isnan(x) && isnan(y)) || isapprox(x, y; atol)

struct ReactantRaytraceMaterial{N} <: Krang.AbstractMaterial
    subimgs::NTuple{N,Int}
end

function (mat::ReactantRaytraceMaterial)(pix, intersection)
    rs = intersection.rs
    return rs / (one(rs) + rs)
end

@testset "Reactant Extension" begin
    met = Krang.Kerr(0.5)
    reactant_ext = Base.get_extension(Krang, :KrangReactantExt)
    αs = collect(range(-10.0,10.0,10))
    βs = collect(range(-10.0,10.0,10))
    αr = to_rarray(αs)
    βr = to_rarray(βs)
#
    θo = π / 4
#
    base_η = Krang.η.(Ref(met), αs, βs, θo)
    reactant_η = @jit Krang.η.(Ref(met), αr, βr, θo)
    @test Array(reactant_η) ≈ base_η
#
    base_λ = Krang.λ.(Ref(met), αs, θo)
    reactant_λ = @jit Krang.λ.(Ref(met), αr, θo)
    @test Array(reactant_λ) ≈ base_λ
#
    function temp_get_radial_roots_num_reals(met, η, λ)
        return sum(Krang._isreal2, Krang.get_radial_roots(met, η,λ))
    end
    base_numreals = temp_get_radial_roots_num_reals.(Ref(met), base_η, base_λ)#Krang.get_radial_roots.(Ref(met), base_η, base_λ)
    reactant_numreals = @jit temp_get_radial_roots_num_reals.(Ref(met), reactant_η, reactant_λ)#reactant_ext.get_radial_roots(met, reactant_η, reactant_λ)
    @test base_numreals == Array(reactant_numreals)
#
    function temp_get_radial_roots_radial_integrals1(met, η, λ)
        return Krang.radial_inf_integrals(met, Krang.get_radial_roots(met, η,λ))[1]
    end
    function temp_get_radial_roots_radial_integrals2(met, η, λ)
        return Krang.radial_inf_integrals(met, Krang.get_radial_roots(met, η,λ))[2]
    end
    function temp_get_radial_roots_radial_integralsp(met, η, λ)
        return Krang.radial_inf_integrals(met, Krang.get_radial_roots(met, η,λ))[3]
    end
    function temp_get_radial_roots_radial_integralsm(met, η, λ)
        return Krang.radial_inf_integrals(met, Krang.get_radial_roots(met, η,λ))[4]
    end
#
    base_I1 = temp_get_radial_roots_radial_integrals1.(Ref(met), base_η, base_λ)
    base_I2 = temp_get_radial_roots_radial_integrals2.(Ref(met), base_η, base_λ)
    base_Ip = temp_get_radial_roots_radial_integralsp.(Ref(met), base_η, base_λ)
    base_Im = temp_get_radial_roots_radial_integralsm.(Ref(met), base_η, base_λ)
    reactant_I1 = @jit temp_get_radial_roots_radial_integrals1.(Ref(met), reactant_η, reactant_λ)
    reactant_I2 = @jit temp_get_radial_roots_radial_integrals2.(Ref(met), reactant_η, reactant_λ)
    reactant_Ip = @jit temp_get_radial_roots_radial_integralsp.(Ref(met), reactant_η, reactant_λ)
    reactant_Im = @jit temp_get_radial_roots_radial_integralsm.(Ref(met), reactant_η, reactant_λ)
#
    @test reactant_I1 ≈ base_I1
    @test reactant_I2 ≈ base_I2
    @test reactant_Ip ≈ base_Ip
    @test reactant_Im ≈ base_Im
#
    function temp_get_radial_roots_I0_inf(met, η, λ)
        return Krang.Ir_inf(met, Krang.get_radial_roots(met, η,λ))
    end
    base_I0_inf = temp_get_radial_roots_I0_inf.(Ref(met), base_η, base_λ)
    reactant_I0_inf = @jit temp_get_radial_roots_I0_inf.(Ref(met), reactant_η, reactant_λ)
    @test reactant_I0_inf ≈ base_I0_inf
#
    function temp_total_mino_time(met, η, λ)
        return Krang.total_mino_time(met, Krang.get_radial_roots(met, η,λ))
    end
    base_τ_total = temp_total_mino_time.(Ref(met), base_η, base_λ)
    reactant_τ_total = @jit temp_total_mino_time.(Ref(met), reactant_η, base_λ)
    @test reactant_τ_total ≈ base_τ_total
#
    function temp_Iϕ_inf(met, η, λ)
        return Krang.Iϕ_inf(met, Krang.get_radial_roots(met, η, λ), λ)
    end
    base_Iϕ_inf = temp_Iϕ_inf.(Ref(met), base_η, base_λ)
    reactant_Iϕ_inf = @jit temp_Iϕ_inf.(Ref(met), reactant_η, base_λ)
    @test reactant_Iϕ_inf ≈ base_Iϕ_inf
#
    function temp_It_inf(met, η, λ)
        return Krang.It_inf(met, Krang.get_radial_roots(met, η, λ), λ)
    end
    base_It_inf = temp_It_inf.(Ref(met), base_η, base_λ)
    reactant_It_inf = @jit temp_It_inf.(Ref(met), reactant_η, base_λ)
    @test reactant_It_inf ≈ base_It_inf
#
    function temp_absGθo(met, θo, η, λ)
        return Krang._absGθo_Gθhat(met, θo, η, λ)[1]
    end
    base_Gθ = temp_absGθo.(Ref(met), θo, base_η, base_λ)
    reactant_Gθ = @jit temp_absGθo.(Ref(met), θo, reactant_η, reactant_λ)
    @test reactant_Gθ ≈ base_Gθ
#
    function temp_Gθhat(met, θo, η, λ)
        return Krang._absGθo_Gθhat(met, θo, η, λ)[2]
    end
    base_Gθ = temp_Gθhat.(Ref(met), θo, base_η, base_λ)
    reactant_Gθ = @jit temp_Gθhat.(Ref(met), θo, reactant_η, reactant_λ)
    @test reactant_Gθ ≈ base_Gθ
#
    function temp_absGϕo(met, θo, η, λ)
        return Krang._absGϕo_Gϕhat(met, θo, η, λ)[1]
    end
    base_Gϕ = temp_absGϕo.(Ref(met), θo, base_η, base_λ)
    reactant_Gϕ = @jit temp_absGϕo.(Ref(met), θo, reactant_η, reactant_λ)
    @test reactant_Gϕ ≈ base_Gϕ
#
    function temp_Gϕhat(met, θo, η, λ)
        return Krang._absGϕo_Gϕhat(met, θo, η, λ)[2]
    end
    base_Gϕ = temp_Gϕhat.(Ref(met), θo, base_η, base_λ)
    reactant_Gϕ = @jit temp_Gϕhat.(Ref(met), θo, reactant_η, reactant_λ)
    @test reactant_Gϕ ≈ base_Gϕ
#
    function temp_absGto(met, θo, η, λ)
        return Krang._absGto_Gthat(met, θo, η, λ)[1]
    end
    base_Gt = temp_absGto.(Ref(met), θo, base_η, base_λ)
    reactant_Gt = @jit temp_absGto.(Ref(met), θo, reactant_η, reactant_λ)
    @test reactant_Gt ≈ base_Gt
#
    function temp_Gthat(met, θo, η, λ)
        return Krang._absGto_Gthat(met, θo, η, λ)[2]
    end
    base_Gt = temp_Gthat.(Ref(met), θo, base_η, base_λ)
    reactant_Gt = @jit temp_Gthat.(Ref(met), θo, reactant_η, reactant_λ)
    @test reactant_Gt ≈ base_Gt
end
@testset "Reactant Emission Radius" begin
    metric = Krang.Kerr(0.01)
    θs = π / 2
    θo = 0.01 / 180 * π
    αvals = [4.5, 5.0]
    βvals = [4.5, 5.0]
#
    emission_radius_rs(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_radius(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[1]
    emission_radius_νr(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_radius(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[2]
    emission_radius_νθ(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_radius(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[3]
    emission_radius_numreals(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_radius(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[4]
    emission_radius_success(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_radius(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[5]
#
    base_rs = emission_radius_rs.(Ref(metric), αvals, βvals, θo, θs, true, 0)
    base_νr = emission_radius_νr.(Ref(metric), αvals, βvals, θo, θs, true, 0)
    base_νθ = emission_radius_νθ.(Ref(metric), αvals, βvals, θo, θs, true, 0)
    base_numreals = emission_radius_numreals.(Ref(metric), αvals, βvals, θo, θs, true, 0)
    base_success = emission_radius_success.(Ref(metric), αvals, βvals, θo, θs, true, 0)
#
    αr = to_rarray(αvals)
    βr = to_rarray(βvals)
    reactant_rs = @jit emission_radius_rs.(Ref(metric), αr, βr, θo, θs, true, 0)
    reactant_νr = @jit emission_radius_νr.(Ref(metric), αr, βr, θo, θs, true, 0)
    reactant_νθ = @jit emission_radius_νθ.(Ref(metric), αr, βr, θo, θs, true, 0)
    reactant_numreals = @jit emission_radius_numreals.(Ref(metric), αr, βr, θo, θs, true, 0)
    reactant_success = @jit emission_radius_success.(Ref(metric), αr, βr, θo, θs, true, 0)
#
    @test Array(reactant_rs) ≈ base_rs
    @test Array(reactant_νr) == base_νr
    @test Array(reactant_νθ) == base_νθ
    @test Array(reactant_numreals) == base_numreals
    @test Array(reactant_success) == base_success
end

@testset "Reactant Emission Coordinates" begin
    metric = Krang.Kerr(0.01)
    θo = 0.01 / 180 * π
    θs = π / 2
    αvals = [4.5, 5.0]
    βvals = [4.5, 5.0]

    emission_coordinates_t(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[1]
    emission_coordinates_r(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[2]
    emission_coordinates_ϕ(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[3]
    emission_coordinates_νr(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[4]
    emission_coordinates_νθ(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[5]
    emission_coordinates_success(met, α, β, θobs, θsrc, isindir, n) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), θsrc, isindir, n)[6]

    base_t = emission_coordinates_t.(Ref(metric), αvals, βvals, θo, θs, true, 0)
    base_r = emission_coordinates_r.(Ref(metric), αvals, βvals, θo, θs, true, 0)
    base_ϕ = emission_coordinates_ϕ.(Ref(metric), αvals, βvals, θo, θs, true, 0)
    base_νr = emission_coordinates_νr.(Ref(metric), αvals, βvals, θo, θs, true, 0)
    base_νθ = emission_coordinates_νθ.(Ref(metric), αvals, βvals, θo, θs, true, 0)
    base_success = emission_coordinates_success.(Ref(metric), αvals, βvals, θo, θs, true, 0)

    αr = to_rarray(αvals)
    βr = to_rarray(βvals)
    reactant_t = @jit emission_coordinates_t.(Ref(metric), αr, βr, θo, θs, true, 0)
    reactant_r = @jit emission_coordinates_r.(Ref(metric), αr, βr, θo, θs, true, 0)
    reactant_ϕ = @jit emission_coordinates_ϕ.(Ref(metric), αr, βr, θo, θs, true, 0)
    reactant_νr = @jit emission_coordinates_νr.(Ref(metric), αr, βr, θo, θs, true, 0)
    reactant_νθ = @jit emission_coordinates_νθ.(Ref(metric), αr, βr, θo, θs, true, 0)
    reactant_success = @jit emission_coordinates_success.(Ref(metric), αr, βr, θo, θs, true, 0)

    @test Array(reactant_t) ≈ base_t
    @test Array(reactant_r) ≈ base_r
    @test Array(reactant_ϕ) ≈ base_ϕ
    @test Array(reactant_νr) == base_νr
    @test Array(reactant_νθ) == base_νθ
    @test Array(reactant_success) == base_success
end

@testset "Reactant Emission Coordinates From Tau" begin
    metric = Krang.Kerr(0.5)
    θo = 0.01 / 180 * π
    θs = π / 2
    αvals = [4.5, 5.0]
    βvals = [4.5, 5.0]

    τvals = map(αvals, βvals) do α, β
        pix = Krang.SlowLightIntensityPixel(metric, α, β, θo)
        Krang.Gθ(pix, θs, true, 0)[1]
    end

    emission_coordinates_τ_t(met, α, β, θobs, τ) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), τ)[1]
    emission_coordinates_τ_r(met, α, β, θobs, τ) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), τ)[2]
    emission_coordinates_τ_θ(met, α, β, θobs, τ) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), τ)[3]
    emission_coordinates_τ_ϕ(met, α, β, θobs, τ) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), τ)[4]
    emission_coordinates_τ_νr(met, α, β, θobs, τ) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), τ)[5]
    emission_coordinates_τ_νθ(met, α, β, θobs, τ) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), τ)[6]
    emission_coordinates_τ_success(met, α, β, θobs, τ) =
        Krang.emission_coordinates(Krang.SlowLightIntensityPixel(met, α, β, θobs), τ)[7]

    base_t = emission_coordinates_τ_t.(Ref(metric), αvals, βvals, θo, τvals)
    base_r = emission_coordinates_τ_r.(Ref(metric), αvals, βvals, θo, τvals)
    base_θ = emission_coordinates_τ_θ.(Ref(metric), αvals, βvals, θo, τvals)
    base_ϕ = emission_coordinates_τ_ϕ.(Ref(metric), αvals, βvals, θo, τvals)
    base_νr = emission_coordinates_τ_νr.(Ref(metric), αvals, βvals, θo, τvals)
    base_νθ = emission_coordinates_τ_νθ.(Ref(metric), αvals, βvals, θo, τvals)
    base_success = emission_coordinates_τ_success.(Ref(metric), αvals, βvals, θo, τvals)

    αr = to_rarray(αvals)
    βr = to_rarray(βvals)
    τr = to_rarray(τvals)
    reactant_t = @jit emission_coordinates_τ_t.(Ref(metric), αr, βr, θo, τr)
    reactant_r = @jit emission_coordinates_τ_r.(Ref(metric), αr, βr, θo, τr)
    reactant_θ = @jit emission_coordinates_τ_θ.(Ref(metric), αr, βr, θo, τr)
    reactant_ϕ = @jit emission_coordinates_τ_ϕ.(Ref(metric), αr, βr, θo, τr)
    reactant_νr = @jit emission_coordinates_τ_νr.(Ref(metric), αr, βr, θo, τr)
    reactant_νθ = @jit emission_coordinates_τ_νθ.(Ref(metric), αr, βr, θo, τr)
    reactant_success = @jit emission_coordinates_τ_success.(Ref(metric), αr, βr, θo, τr)

    @test Array(reactant_t) ≈ base_t
    @test Array(reactant_r) ≈ base_r
    @test Array(reactant_θ) ≈ base_θ
    @test Array(reactant_ϕ) ≈ base_ϕ
    @test Array(reactant_νr) == base_νr
    @test Array(reactant_νθ) == base_νθ
    @test Array(reactant_success) == base_success
end

#@testset "Reactant Synchrotron Intensity" begin
#    metric = Krang.Kerr(0.01)
#    θs = π / 2
#    θo = 0.01 / 180 * π
#    αvals = [4.5, 5.0]
#    βvals = [4.5, 5.0]
#    magfield = @SVector [0.0, 0.0, 1.0]
#    βfluid = @SVector [0.0, 0.0, 0.0]
#
#    rsvals = map(αvals, βvals) do α, β
#        pix = Krang.SlowLightIntensityPixel(metric, α, β, θo)
#        rs, _, _, _, issuccess = Krang.emission_radius(pix, θs, true, 0)
#        @test issuccess
#        rs
#    end
#
#    base_intensity = Krang.synchrotronIntensity.(
#        Ref(metric),
#        αvals,
#        βvals,
#        rsvals,
#        θs,
#        θo,
#        Ref(magfield),
#        Ref(βfluid),
#        true,
#        false,
#    )
#
#    αr = to_rarray(αvals)
#    βr = to_rarray(βvals)
#    rsr = to_rarray(rsvals)
#    reactant_norm, reactant_redshift, reactant_lp = @jit Krang.synchrotronIntensity.(
#        Ref(metric),
#        αr,
#        βr,
#        rsr,
#        θs,
#        θo,
#        Ref(magfield),
#        Ref(βfluid),
#        true,
#        false,
#    )
#
#    @test Array(reactant_norm) ≈ first.(base_intensity)
#    @test Array(reactant_redshift) ≈ getindex.(base_intensity, 2)
#    @test Array(reactant_lp) ≈ getindex.(base_intensity, 3)
#end
#
#
#
#@testset "Broadcasted Intensity Pixel Raytrace" begin
#    met = Krang.Kerr(0.5)
#    αvals = collect(range(8.0, 10.0, 4))
#    βvals = collect(range(8.0, 10.0, 4))
#    αr = to_rarray(αvals)
#    βr = to_rarray(βvals)
#    
#    θray = π / 4
#    mesh = Krang.Mesh(Krang.ConeGeometry(π / 4), ReactantRaytraceMaterial((0, 1, 2)))
#
#    function test_pixel(met, α, β, θo)
#        Krang.SlowLightIntensityPixel(met, α, β, θo).I0_inf
#    end
#    @jit test_pixel.(Ref(met), αr, βr, θray)
#    
#    raytrace_intensity(met, α, β, θo, mesh) =
#        Krang.raytrace(Krang.SlowLightIntensityPixel(met, α, β, θo), mesh)
#    
#    base_obs = raytrace_intensity.(Ref(met), αvals, βvals, θray, Ref(mesh))
#    #reactant_obs = @jit raytrace_intensity.(Ref(met), αr, βr, θray, Ref(mesh))
#    
#    #@test reactant_obs ≈ base_obs
#end
