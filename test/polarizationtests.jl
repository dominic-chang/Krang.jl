@testset "Polarization" begin
    @testset "Null momentum" begin
        metric = Kerr(0.99)
        λtemp = λ(metric, 5.0, π/4)
        ηtemp = η(metric, 5.0, 5.0, π/4)
        emrs,_,_ = emission_radius(metric, 5.0, 5.0, π/4, π/5, true, 0)
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
        met_uu = metric_uu(metric, rs, θs)
        met_dd = metric_dd(metric, rs, θs)
        minkowski = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
        identity = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]


        @testset "Inversion" begin
            @test maximum(abs, (jacbz_u_d * jaczb_u_d) .- identity) ≈ 0 atol=1e-10
            @test maximum(abs, (jacbz_d_u * jaczb_d_u) .- identity) ≈ 0 atol=1e-10
        end

        # Check that the local metric is Minkowski in the ZAMO basis.
        @testset "Metric transformation" begin
            @test maximum(abs, (jaczb_u_d * met_uu * jaczb_u_d') .- minkowski) ≈ 0 atol=1e-10
            @test maximum(abs, (jaczb_d_u * met_dd * jaczb_d_u') .- minkowski) ≈ 0 atol=1e-10
        end
    end
end