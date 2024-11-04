@testset "Polarization" begin
    @testset "Null momentum" begin
        metric = Kerr(0.99)
        α = 5.0
        β = 5.0
        θo = π/4
        pix = Krang.SlowLightIntensityPixel(metric, α, β, θo)
        λtemp = λ(metric, 5.0, π/4)
        ηtemp = η(metric, 5.0, 5.0, π/4)
        emrs,_,_,issuccess = emission_radius(pix, π/5, true, 0)
        p_d = p_bl_d(metric, emrs, π/4, ηtemp, λtemp, true, true)
        met_uu = metric_uu(metric, emrs, π/4)

        # Null momentum magnitude should be 0.
        @test (p_d' * met_uu * p_d) ≈ 0.0 atol=1e-10 
    end

    @testset "Tetrad transformation" begin
        metric = Kerr(0.99)
        rs = 4.0
        θs = π/3
        jacbz_d_u = jac_bl_d_zamo_u(metric, rs, θs)
        jacbz_u_d = jac_bl_u_zamo_d(metric, rs, θs)
        jaczb_d_u = jac_zamo_d_bl_u(metric, rs, θs)
        jaczb_u_d = jac_zamo_u_bl_d(metric, rs, θs)
        jacfz_u_d = jac_fluid_u_zamo_d(metric, 0.5, π/4, π/4)
        jaczf_u_d = jac_fluid_u_zamo_d(metric, -0.5, π/4, π/4)
        trivjacfz_u_d = jac_fluid_u_zamo_d(metric, 0.0, 0.0, 0.0)
        met_uu = metric_uu(metric, rs, θs)
        met_dd = metric_dd(metric, rs, θs)
        minkowski = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
        identity = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

        @test maximum(abs, trivjacfz_u_d .- identity) ≈ 0 atol=1e-10
        @testset "Inversion" begin
            @test maximum(abs, (jacbz_u_d * jaczb_u_d) .- identity) ≈ 0 atol=1e-10
            @test maximum(abs, (jacbz_d_u * jaczb_d_u) .- identity) ≈ 0 atol=1e-10
            @test maximum(abs, (jacfz_u_d * jaczf_u_d) .- identity) ≈ 0 atol=1e-10
        end

        # Check that the local metric is Minkowski in the ZAMO basis.
        @testset "Metric transformation" begin
            @test maximum(abs, (jaczb_u_d * met_uu * jaczb_u_d') .- minkowski) ≈ 0 atol=1e-10
            @test maximum(abs, (jaczb_d_u * met_dd * jaczb_d_u') .- minkowski) ≈ 0 atol=1e-10
        end

        @testset "Polarization routine" begin
            metric = Kerr(0.01)
            θs = π/2
            θo = 0.01/180*π
            α = 5.0
            β = 5.0
            pix = Krang.SlowLightIntensityPixel(metric, α, β, θo)
            λtemp = λ(metric, α, θo)
            ηtemp = η(metric, α, β, θo)
            rs,_,_,issuccess = emission_radius(pix, θs, true, 0)

            magfield = @SVector[0.,0.,1.0]
            βfluid = @SVector[0.,0.,.0]

            eα, eβ, redshift1, lp1 = Krang.synchrotronPolarization(metric, α, β, rs, θs, θo, magfield, βfluid, true, false)
            i, redshift2, lp2 = Krang.synchrotronIntensity(metric, α, β, rs, θs, θo, magfield, βfluid, true, false)
            @test evpa(eα, eβ) ≈ -3π/4 atol = 1e-3
            @test hypot(eα, eβ)/i ≈ 1 atol = 1e-3
            @test redshift1/redshift2 ≈ 1.0 atol = 1e-3
            @test lp1/lp2 ≈ 1.0 atol = 1e-3
        end

        @testset "Polarization pixel interface"  begin
            metric = Kerr(0.01)
            θs = π/2
            θo = 0.01/180*π
            α = 5.0
            β = 5.0
            pix = Krang.SlowLightIntensityPixel(metric, α, β, θo)
            λtemp = λ(metric, α, θo)
            ηtemp = η(metric, α, β, θo)
            rs,_,_,issuccess = emission_radius(pix, θs, true, 0)

            magfield = @SVector[0.,0.,1.0]
            βfluid = @SVector[0.,0.,0.]

            @testset "Intensity" begin
                subimgs = (0, )
                spec = 1.0
                rpeak = 5.0
                p1 = 1.0
                p2 = 1.0

                mat = Krang.ElectronSynchrotronPowerLawIntensity(magfield..., βfluid..., spec, rpeak, p1, p2, subimgs)
                geometry = Krang.ConeGeometry((θs))

                int = render(pix, (Krang.Mesh(geometry, mat),))

                profile(r) = (r/rpeak)^p1/(1+(r/rpeak)^(p1+p2))
                norm, redshift, lp = Krang.synchrotronIntensity(metric, α, β, rs, θs, θo, magfield, βfluid, true, false)
                int2 = norm^(1 + spec) * lp * profile(rs)*redshift^(3 + spec)
                @test int ≈ int2 atol = 1e-3
            end

            @testset "Polarization" begin
                @inline profile(r) = 1               
                subimgs = (0, )
                spec = 1.0
                rpeak = 5.0
                p1 = 1.0
                p2 = 1.0

                mat = Krang.ElectronSynchrotronPowerLawPolarization(magfield..., βfluid..., spec, rpeak, p1, p2, subimgs)
                geometry = Krang.ConeGeometry((θs))

                stokes = render(pix, (Krang.Mesh(geometry, mat),))

                eα, eβ, redshift, lp = Krang.synchrotronPolarization(metric, α, β, rs, θs, θo, magfield, βfluid, true, false)
                q = (-(eα^2 - eβ^2) )
                u = (-2 * eα * eβ )

                profile(r) = (r/rpeak)^p1/(1+(r/rpeak)^(p1+p2))
                int2 = hypot(q,u)^(1 + spec) * lp * profile(rs)*redshift^(3 + spec)
                @test stokes ≈ Krang.StokesParams(int2, q, u, 0.0)
            end
        end
    end
    @testset "Penrose Walker Constant" begin
        @test all(penrose_walker(Krang.Kerr(0.0), 1.0, 0.0, [1.0,0.0,0.0,0.0], [0.0,1.0,0.0,0.0]) .== (1.0, 0.0))
    end
end