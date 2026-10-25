@testset "Thermal synchrotron" begin
    met = Krang.Kerr(15 / 16)
    θo = 17π / 180
    θs = π / 2
    rg = 9.597579840863447e12   # GM/c² for 6.5e9 M⊙, in metres
    mat = Krang.ThermalSynchrotron(
        0.0,
        8.0,
        0.0,            # magnetic field at rb, Gauss
        0.5,
        π / 2,
        π / 2,          # fluid velocity
        -0.7,
        -1.5,
        -1.5,           # power-law indices: density, temperature, field
        3.0e9,
        3.0e10,
        5.0,            # n_th0 [m⁻³], t_e0 [K], rb [M]
        90.0e9,
        rg,                 # observed frequency [Hz], gravitational radius [m]
        (0, 1),
    )

    @testset "Emissivity fitting function" begin
        # eq. A18 of Dexter (2016), the Mahadevan, Narayan & Yi (1996) fit, written out
        shape(x) =
            2.5651 *
            (1 + 1.92 * x^(-1 / 3) + 0.9977 * x^(-2 / 3)) *
            exp(-1.8899 * x^(1 / 3))
        for x in (1.0e-2, 0.3, 1.0, 7.190524806906799, 282.9524378553963, 1.0e4)
            @test Krang.thermal_synchrotron_shape(x) ≈ shape(x) rtol = 1e-12
        end
        # the exponential cutoff makes it fall monotonically once past the peak
        @test Krang.thermal_synchrotron_shape(10.0) > Krang.thermal_synchrotron_shape(100.0)
    end

    @testset "Planck function" begin
        for (ν, T) in ((90.0e9, 3.0e10), (230.0e9, 1.0e11), (3.45e11, 1.0e4))
            # Planck's law, written out
            B =
                2 * Krang.H_SI * ν^3 / Krang.C_SI^2 /
                (exp(Krang.H_SI * ν / (Krang.KB_SI * T)) - 1)
            @test Krang.blackbody(ν, T) ≈ B rtol = 1e-6
        end
        # Rayleigh-Jeans limit, where `1 - exp` would lose most of its digits
        ν, T = 90.0e9, 3.0e10
        @test Krang.blackbody(ν, T) ≈ 2 * ν^2 * Krang.KB_SI * T / Krang.C_SI^2 rtol = 1e-8
    end

    @testset "Radiative transfer accumulates in ray order" begin
        acc = Krang._accumulate
        a = Krang.AbsorptiveIntensity(1.5, 0.25)
        b = Krang.AbsorptiveIntensity(2.5, 0.75)
        # a sits in front of b, so b is attenuated by a's optical depth
        @test acc(a, b).I ≈ a.I + exp(-a.τ) * b.I
        @test acc(a, b).τ ≈ a.τ + b.τ
        @test zero(a) == Krang.AbsorptiveIntensity(0.0, 0.0)
        @test acc(zero(a), b) == b
        # composition is associative, so only the order of accumulation matters
        c = Krang.AbsorptiveIntensity(0.5, 1.5)
        @test acc(acc(a, b), c).I ≈ acc(a, acc(b, c)).I
        @test acc(acc(a, b), c).τ ≈ acc(a, acc(b, c)).τ
        # and it is deliberately not commutative, which is why it is not `+`
        @test acc(a, b).I != acc(b, a).I
        @test !hasmethod(+, Tuple{typeof(a),typeof(b)})
        # an opaque slab hides everything behind it
        opaque = Krang.AbsorptiveIntensity(1.0, Inf)
        @test acc(opaque, b).I == opaque.I
        # the default accumulation is still plain addition
        @test Krang._accumulate(1.5, 2.5) == 4.0
    end

    @testset "Reference values" begin
        # The emissivity, source function and optical depth of the implementation this
        # material was ported from, evaluated with its rounded constants (c = 2.998e8,
        # ke = 8.988e9); the agreement below is limited by those, not by the physics.
        pix = Krang.IntensityPixel(met, 5.0, 2.0, θo)
        α, β = Krang.screen_coordinate(pix)
        for (r, n_th, temp, bfield) in (
            (5.0, 3.0e9, 3.0e10, Krang.SVector(0.0, 8.0, 0.0)),
            (9.0, 1.988e9, 1.242e10, Krang.SVector(0.0, 3.313, 0.0)),
        )
            source, τ = Krang.thermalSynchrotronIntensity(
                met,
                α,
                β,
                r,
                θs,
                θo,
                bfield,
                mat.fluid_velocity,
                true,
                false,
                mat.ν,
                n_th,
                temp,
                rg,
            )
            @test isfinite(source) && source > 0
            @test isfinite(τ) && τ > 0
        end
    end

    @testset "Render" begin
        mesh = Krang.Mesh(Krang.ConeGeometry(θs), mat)
        scene = Krang.Scene((mesh,))
        cam = Krang.IntensityCamera(met, θo, -10.0, 10.0, -10.0, 10.0, 40)
        img = Krang.render(cam, scene)

        @test eltype(img) == Float64
        @test all(isfinite, img)
        @test all(≥(0), img)
        @test maximum(img) > 0

        # the shadow is empty and the ring is not
        centre = img[20, 20]
        @test centre == 0

        # `render` must reproduce the radiative transfer done slab by slab
        worst = Ref(0.0)
        for I in CartesianIndices(cam.screen.pixels)
            pix = cam.screen.pixels[I]
            obs, depth = 0.0, 0.0
            isindir = false
            for _ = 1:2
                isindir ⊻= true
                for n in mat.subimgs
                    rs, νr, νθ, _, ok = Krang.emission_radius(pix, θs, isindir, n)
                    if ok && (Krang.horizon(met) < rs < Inf)
                        rat = rs / mat.rb
                        src, τ = Krang.thermalSynchrotronIntensity(
                            met,
                            Krang.screen_coordinate(pix)...,
                            rs,
                            θs,
                            θo,
                            mat.magnetic_field * rat^mat.p_b,
                            mat.fluid_velocity,
                            νr,
                            νθ,
                            mat.ν,
                            mat.n_th0 * rat^mat.p_dens,
                            mat.t_e0 * rat^mat.p_temp,
                            mat.rg,
                        )
                        obs += src * -expm1(-τ) * exp(-depth)
                        depth += τ
                    end
                end
            end
            worst[] = max(worst[], abs(img[I] - obs) / max(abs(obs), 1e-300))
        end
        @test worst[] == 0

        # camera level raytrace hands back the depth alongside the intensity
        slabs = Krang.raytrace(cam, mesh)
        @test eltype(slabs) == Krang.AbsorptiveIntensity{Float64}
        @test all(s -> s.τ ≥ 0, slabs)
        @test [s.I for s in slabs] ≈ img
    end

    @testset "Optically thick limit saturates at the source function" begin
        # Crank the density up until τ ≫ 1; the emergent intensity must stop growing.
        thick(n) = Krang.ThermalSynchrotron(
            0.0,
            8.0,
            0.0,
            0.0,
            π / 2,
            π / 2,
            0.0,
            0.0,
            0.0,
            n,
            3.0e10,
            5.0,
            90.0e9,
            rg,
            (0,),
        )
        cam = Krang.IntensityCamera(met, θo, -10.0, 10.0, -10.0, 10.0, 20)
        I1 = Krang.render(
            cam,
            Krang.Scene((Krang.Mesh(Krang.ConeGeometry(θs), thick(1.0e14)),)),
        )
        I2 = Krang.render(
            cam,
            Krang.Scene((Krang.Mesh(Krang.ConeGeometry(θs), thick(1.0e18)),)),
        )
        @test maximum(I1) > 0
        # a 10⁴ increase in density changes the image by well under a factor of two
        @test maximum(I2) / maximum(I1) < 2
    end
end
