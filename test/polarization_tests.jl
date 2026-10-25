@testset "Polarization" begin
    @testset "Null momentum" begin
        metric = Kerr(0.99)
        α = 5.0
        β = 5.0
        θo = π / 4
        pix = Krang.SlowLightIntensityPixel(metric, α, β, θo)
        λtemp = λ(metric, 5.0, π / 4)
        ηtemp = η(metric, 5.0, 5.0, π / 4)
        emrs, _, _, issuccess = emission_radius(pix, π / 5, true, 0)
        p_d = p_bl_d(metric, emrs, π / 4, ηtemp, λtemp, true, true)
        met_uu = metric_uu(metric, emrs, π / 4)

        # Null momentum magnitude should be 0.
        @test (p_d' * met_uu * p_d) ≈ 0.0 atol = 1e-10
    end

    @testset "Tetrad transformation" begin
        metric = Kerr(0.99)
        rs = 4.0
        θs = π / 3
        jacbz_d_u = jac_bl_d_zamo_u(metric, rs, θs)
        jacbz_u_d = jac_bl_u_zamo_d(metric, rs, θs)
        jaczb_d_u = jac_zamo_d_bl_u(metric, rs, θs)
        jaczb_u_d = jac_zamo_u_bl_d(metric, rs, θs)
        jacfz_u_d = jac_fluid_u_zamo_d(metric, 0.5, π / 4, π / 4)
        jaczf_u_d = jac_fluid_u_zamo_d(metric, -0.5, π / 4, π / 4)
        met_uu = metric_uu(metric, rs, θs)
        met_dd = metric_dd(metric, rs, θs)
        minkowski = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
        identity = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

        @testset "Polarization routine" begin
            metric = Kerr(0.01)
            θs = π / 2
            θo = 0.01 / 180 * π
            α = 5.0
            β = 5.0
            pix = Krang.SlowLightIntensityPixel(metric, α, β, θo)
            λtemp = λ(metric, α, θo)
            ηtemp = η(metric, α, β, θo)
            rs, _, _, issuccess = emission_radius(pix, θs, true, 0)

            magfield = @SVector[0.0, 0.0, 1.0]
            βfluid = @SVector[0.0, 0.0, 0.0]

            eα, eβ, redshift1, lp1 = Krang.synchrotronPolarization(
                metric,
                α,
                β,
                rs,
                θs,
                θo,
                magfield,
                βfluid,
                true,
                false,
            )
            i, redshift2, lp2 = Krang.synchrotronIntensity(
                metric,
                α,
                β,
                rs,
                θs,
                θo,
                magfield,
                βfluid,
                true,
                false,
            )
            @test evpa(eα, eβ) ≈ -3π / 4 atol = 1e-3
            @test hypot(eα, eβ) / i ≈ 1 atol = 1e-3
            @test redshift1 / redshift2 ≈ 1.0 atol = 1e-3
            @test lp1 / lp2 ≈ 1.0 atol = 1e-3
        end

        @testset "Polarization pixel interface" begin
            metric = Kerr(0.01)
            θs = π / 2
            θo = 0.01 / 180 * π
            α = 5.0
            β = 5.0
            pix = Krang.SlowLightIntensityPixel(metric, α, β, θo)
            λtemp = λ(metric, α, θo)
            ηtemp = η(metric, α, β, θo)
            rs, _, _, issuccess = emission_radius(pix, θs, true, 0)

            magfield = @SVector[0.0, 0.0, 1.0]
            βfluid = @SVector[0.0, 0.0, 0.0]

            @testset "Intensity" begin
                subimgs = (0,)
                spec = 1.0
                rpeak = 5.0
                p1 = 1.0
                p2 = 1.0

                mat = Krang.ElectronSynchrotronPowerLawIntensity(
                    magfield...,
                    βfluid...,
                    spec,
                    rpeak,
                    p1,
                    p2,
                    subimgs,
                )
                geometry = Krang.ConeGeometry((θs))

                int = render(pix, (Krang.Mesh(geometry, mat),))

                profile(r) = (r / rpeak)^p1 / (1 + (r / rpeak)^(p1 + p2))
                norm, redshift, lp = Krang.synchrotronIntensity(
                    metric,
                    α,
                    β,
                    rs,
                    θs,
                    θo,
                    magfield,
                    βfluid,
                    true,
                    false,
                )
                int2 = norm^(1 + spec) * lp * profile(rs) * redshift^(3 + spec)
                @test int ≈ int2 atol = 1e-3
            end

            @testset "Polarization" begin
                @inline profile(r) = 1
                subimgs = (0,)
                spec = 1.0
                rpeak = 5.0
                p1 = 1.0
                p2 = 1.0

                mat = Krang.ElectronSynchrotronPowerLawPolarization(
                    magfield...,
                    βfluid...,
                    spec,
                    rpeak,
                    p1,
                    p2,
                    subimgs,
                )
                geometry = Krang.ConeGeometry((θs))

                stokes = render(pix, (Krang.Mesh(geometry, mat),))

                eα, eβ, redshift, lp = Krang.synchrotronPolarization(
                    metric,
                    α,
                    β,
                    rs,
                    θs,
                    θo,
                    magfield,
                    βfluid,
                    true,
                    false,
                )
                q_temp = (-(eα^2 - eβ^2))
                u_temp = (-2 * eα * eβ)

                profile(r) = (r / rpeak)^p1 / (1 + (r / rpeak)^(p1 + p2))
                mag = hypot(q_temp, u_temp)
                int2 = mag^(1 + spec) * lp * profile(rs) * redshift^(3 + spec)
                q = int2 * q_temp / mag
                u = int2 * u_temp / mag
                @test stokes ≈ Krang.StokesParams(int2, q, u, 0.0)
            end
        end
    end
    @testset "Expanded frame transformations" begin
        # Each frame transformation has a method that returns the matrix and one that
        # applies it to a vector without forming it; they must agree.
        @testset "a=$a" for a in (0.01, 0.5, 0.94, -0.8)
            metric = Kerr(a)
            @testset "r=$r, θ=$θ" for r in (2.5, 4.0, 8.0), θ in (0.4, π / 3, π / 2, 2π / 3)
                v_d = @SVector [1.7, -0.6, 0.9, -1.3]
                v_u = @SVector [1.7, -0.6, 0.9, -1.3]
                @test jac_zamo_d_bl_u(metric, r, θ, v_d) ≈
                      jac_zamo_d_bl_u(metric, r, θ) * v_d
                @test jac_bl_d_zamo_u(metric, r, θ, v_d) ≈
                      jac_bl_d_zamo_u(metric, r, θ) * v_d
                @test jac_zamo_u_bl_d(metric, r, θ, v_u) ≈
                      jac_zamo_u_bl_d(metric, r, θ) * v_u
                @test jac_bl_u_zamo_d(metric, r, θ, v_u) ≈
                      jac_bl_u_zamo_d(metric, r, θ) * v_u

                @testset "boost β=$βv, θz=$θz, ϕz=$ϕz" for βv in (0.0, 0.4, -0.8),
                    θz in (0.0, π / 4, 1.1),
                    ϕz in (0.0, π / 3, -0.5)

                    @test jac_fluid_u_zamo_d(metric, βv, θz, ϕz, v_u) ≈
                          jac_fluid_u_zamo_d(metric, βv, θz, ϕz) * v_u
                end

                # `p_zamo_u` collapses metric_uu → ZAMO into one step
                met_uu = metric_uu(metric, r, θ)
                jaczb_u_d = jac_zamo_u_bl_d(metric, r, θ)
                for ηtemp in (1.0, 20.0),
                    λtemp in (-6.0, 0.0, 3.5),
                    νr in (true, false),
                    νθ in (true, false)

                    p_d = p_bl_d(metric, r, θ, ηtemp, λtemp, νr, νθ)
                    @test Krang.p_zamo_u(metric, r, θ, ηtemp, λtemp, νr, νθ) ≈
                          jaczb_u_d * (met_uu * p_d)
                end
            end
        end

        @testset "κ is the Penrose-Walker constant" begin
            # `synchrotronPolarization` used to build κ as `(A - iB)(r - ia cosθ)`
            # contracted with the momentum; that product is exactly `penrose_walker`.
            @testset "a=$a" for a in (0.01, 0.5, 0.94, -0.8)
                metric = Kerr(a)
                for r in (2.5, 4.0, 8.0), θ in (0.4, π / 3, π / 2, 2π / 3)
                    # BL-ordered {t, r, θ, ϕ} vectors, as the matrices expect
                    p_u = @SVector [1.7, -0.6, 0.9, -1.3]
                    f_u = @SVector [0.4, 1.1, -0.2, 0.8]
                    z, o = 0.0, 1.0
                    A = @SMatrix [
                        z o z z
                        -o z z a*sin(θ)^2
                        z z z z
                        z -a*sin(θ)^2 z z
                    ]
                    B = @SMatrix [
                        z z -a*sin(θ) z
                        z z z z
                        a*sin(θ) z z -(r^2 + a^2)*sin(θ)
                        z z (r^2+a^2)*sin(θ) z
                    ]
                    κ_matrix = sum(p_u .* (((A - B * im) * (r - a * cos(θ) * im)) * f_u))
                    # `penrose_walker` orders its arguments {t, r, ϕ, θ}
                    κ1, κ2 = penrose_walker(
                        metric,
                        r,
                        θ,
                        (@SVector [p_u[1], p_u[2], p_u[4], p_u[3]]),
                        (@SVector [f_u[1], f_u[2], f_u[4], f_u[3]]),
                    )
                    @test κ1 ≈ real(κ_matrix)
                    @test κ2 ≈ imag(κ_matrix)
                end
            end
        end
    end
    @testset "Penrose Walker Constant" begin
        @test all(
            penrose_walker(
                Krang.Kerr(0.0),
                1.0,
                0.0,
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
            ) .== (1.0, 0.0),
        )
    end
end
