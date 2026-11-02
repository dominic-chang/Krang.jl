using KernelAbstractions

"Cheap axisymmetric cone material, so the render test exercises the kernel and not a material."
struct KAConeMaterial{N,T} <: Krang.AbstractMaterial
    subimgs::NTuple{N,Int}
    weight::T
end
Krang.lightTrait(::KAConeMaterial) = Krang.FastLightTrait()
Krang.symmetryTrait(::KAConeMaterial) = Krang.AxisymmetricTrait()
function (mat::KAConeMaterial)(pix, intersection)
    α, β = Krang.screen_coordinate(pix)
    return mat.weight * (intersection.rs + intersection.θs + α - β)
end

# Metal is optional: on a machine without it the suite still exercises the CPU backend, which
# is the one that has to keep working everywhere.
# Metal kernels here take minutes to compile, so they are opt-in: set
# KRANG_KA_METAL_TESTS=true to include them.
const ka_metal_loaded =
    if Sys.isapple() && get(ENV, "KRANG_KA_METAL_TESTS", "false") == "true"
        try
            @eval using Metal
            Metal.functional()
        catch _
            false
        end
    else
        false
    end

@testset "Kernel Abstractions" begin
    Base.retry_load_extensions()
    @test !isnothing(Base.get_extension(Krang, :KrangKernelAbstractionsExt))

    met = Krang.Kerr(0.3)
    θo = π / 3
    αmin, αmax = -2.0, 2.0
    βmin, βmax = -1.5, 1.5
    res = 4

    "Backends to sweep, with the array type and the tolerance its precision supports."
    backends = Any[("CPU", Matrix, Float64, 0.0)]
    ka_metal_loaded && push!(backends, ("Metal", Metal.MtlMatrix, Float32, 1e-4))

    for (name, Mat, T, atol) in backends
        metT = Krang.Kerr(T(0.3))
        cmp(a, b) = atol == 0 ? a == b : isapprox(a, b; atol = atol, rtol = atol)

        @testset "$name IntensityScreen constructor" begin
            serial = Krang.IntensityScreen(metT, αmin, αmax, βmin, βmax, θo, res)
            kern = Krang.IntensityScreen(metT, αmin, αmax, βmin, βmax, θo, res, Mat)

            @test size(kern.pixels) == size(serial.pixels)
            @test kern.αrange == serial.αrange
            @test kern.βrange == serial.βrange

            # Materialise on the host once; indexing a device-backed view per element would
            # scalar-index the backing arrays, which GPU array types forbid.
            host = Array(kern.pixels)
            for idx in CartesianIndices(serial.pixels)
                want, got = serial.pixels[idx], host[idx]
                @test cmp(Krang.screen_coordinate(got)[1], Krang.screen_coordinate(want)[1])
                @test cmp(Krang.screen_coordinate(got)[2], Krang.screen_coordinate(want)[2])
                @test Krang.inclination(got) == Krang.inclination(want)
                @test cmp(Krang.η(got), Krang.η(want))
                @test cmp(Krang.λ(got), Krang.λ(want))
                @test maximum(
                    abs.(collect(Krang.roots(got)) .- collect(Krang.roots(want))),
                ) <= max(atol, 1e-10)
                @test cmp(Krang.I0_inf(got), Krang.I0_inf(want))
                @test cmp(Krang.total_mino_time(got), Krang.total_mino_time(want))
            end
        end

        # The slow-light kernel does not finish compiling on Metal; see the extension's
        # docstring. Exercised on CPU only so the gap stays visible rather than hanging CI.
        # No `return` to skip: a `@testset` body is not a function, so a `return` escapes the
        # *enclosing* testset and silently drops everything after it.
        @testset "$name SlowLightIntensityScreen constructor" begin
            if name != "CPU"
                @test_broken false  # slow-light screen kernel hangs the Metal compiler
            else
                serial =
                    Krang.SlowLightIntensityScreen(metT, αmin, αmax, βmin, βmax, θo, res)
                kern = Krang.SlowLightIntensityScreen(
                    metT,
                    αmin,
                    αmax,
                    βmin,
                    βmax,
                    θo,
                    res,
                    Mat,
                )

                @test size(kern.pixels) == size(serial.pixels)
                host = Array(kern.pixels)
                for idx in CartesianIndices(serial.pixels)
                    want, got = serial.pixels[idx], host[idx]
                    @test cmp(Krang.I0_inf(got), Krang.I0_inf(want))
                    @test cmp(Krang.total_mino_time(got), Krang.total_mino_time(want))
                    @test cmp(got.Iϕ_inf, want.Iϕ_inf)
                    @test cmp(got.It_inf, want.It_inf)
                    @test cmp(got.absGθo_Gθhat[1], want.absGθo_Gθhat[1])
                    @test cmp(got.absGϕo_Gϕhat[2], want.absGϕo_Gϕhat[2])
                    @test cmp(got.absGto_Gthat[2], want.absGto_Gthat[2])
                end
            end
        end
    end

    @testset "render! matches render" begin
        camera = Krang.IntensityCamera(met, θo, αmin, αmax, βmin, βmax, res)
        scene = Krang.Scene((
            Krang.Mesh(Krang.ConeGeometry(π / 4), KAConeMaterial((0, 1), 0.125)),
        ))

        want = render(camera, scene)
        store = similar(want)
        got = render!(store, camera, scene)

        @test got === store
        @test got == want

        # The size check must fire rather than launching a mismatched kernel.
        @test_throws DimensionMismatch render!(similar(want, res + 1, res), camera, scene)
    end
end
