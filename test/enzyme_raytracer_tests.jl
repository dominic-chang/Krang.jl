using EnzymeTestUtils: test_forward, test_reverse

@testset "Enzyme Raytracer Functions" begin
    met = Krang.Kerr(0.99)
    tBL = 10.0
    rBL = 8.0
    θBL = π / 4
    ϕBL = π / 5

    @testset "Boyer-Lindquist to Kerr-Schild" begin
        f_r(r) = begin
            t, x, y, z = Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild(
                met,
                tBL,
                r,
                θBL,
                ϕBL,
            )
            return t + 2x - 3y + 4z
        end

        f_θ(θ) = begin
            t, x, y, z = Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild(
                met,
                tBL,
                rBL,
                θ,
                ϕBL,
            )
            return t + 2x - 3y + 4z
        end

        f_ϕ(ϕ) = begin
            t, x, y, z = Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild(
                met,
                tBL,
                rBL,
                θBL,
                ϕ,
            )
            return t + 2x - 3y + 4z
        end

        test_forward(
            f_r,
            Enzyme.DuplicatedNoNeed,
            (rBL, Enzyme.Duplicated);
            rtol = 1e-5,
            atol = 1e-6,
        )
        test_forward(
            f_θ,
            Enzyme.DuplicatedNoNeed,
            (θBL, Enzyme.Duplicated);
            rtol = 1e-5,
            atol = 1e-6,
        )
        test_forward(
            f_ϕ,
            Enzyme.DuplicatedNoNeed,
            (ϕBL, Enzyme.Duplicated);
            rtol = 1e-5,
            atol = 1e-6,
        )
        test_reverse(f_r, Enzyme.Active, (rBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_θ, Enzyme.Active, (θBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_ϕ, Enzyme.Active, (ϕBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
    end

    @testset "Horizon-colliding root stays differentiable" begin
        # On the screen line λ = (r₋² + a²)/a a radial root lands on the inner horizon and
        # Iϕ_inf used to evaluate to NaN. Reverse mode then poisoned the whole gradient
        # through 0 × NaN partials, even where the value itself is multiplied by zero.
        aspin = 0.1488920854854475
        θline = deg2rad(165.1348867591876)
        αline, βline = -0.5836283673902987, 0.26401133520389547

        Iϕline(q) = begin
            m = Krang.Kerr(q[3])
            ηq, λq = Krang.η(m, q[1], q[2], q[4]), Krang.λ(m, q[1], q[4])
            return Krang.Iϕ_inf(m, Krang.get_radial_roots(m, ηq, λq), λq)
        end
        unused(q) =
            Krang.total_mino_time(
                Krang.Kerr(q[3]),
                Krang.get_radial_roots(
                    Krang.Kerr(q[3]),
                    Krang.η(Krang.Kerr(q[3]), q[1], q[2], q[4]),
                    Krang.λ(Krang.Kerr(q[3]), q[1], q[4]),
                ),
            ) + 0.0 * Iϕline(q)

        q0 = [αline, βline, aspin, θline]
        for g in (Iϕline, unused)
            dq = zeros(4)
            Enzyme.autodiff(Enzyme.Reverse, g, Enzyme.Active, Enzyme.Duplicated(q0, dq))
            @test all(isfinite, dq)
        end

        # and the surviving gradient is the right one
        dq = zeros(4)
        Enzyme.autodiff(Enzyme.Reverse, Iϕline, Enzyme.Active, Enzyme.Duplicated(q0, dq))
        h = 1e-6
        fd = (Iϕline(q0 .+ [h, 0, 0, 0]) - Iϕline(q0 .- [h, 0, 0, 0])) / 2h
        @test dq[1] ≈ fd rtol = 1e-5
    end

    @testset "Cardano cancellation stays differentiable" begin
        # Where the resolvent cubic's cube-root term cancels, every radial root came back
        # NaN, which took the gradient of anything built from the pixel with it.
        rootsum(q) = begin
            m = Krang.Kerr(q[3])
            ηq, λq = Krang.η(m, q[1], q[2], q[4]), Krang.λ(m, q[1], q[4])
            return sum(x -> real(x) + imag(x)^2, Krang.get_radial_roots(m, ηq, λq))
        end
        q0 = [
            -1.172784384687551,
            -0.1443508333390526,
            0.9387019683755741,
            deg2rad(94.91559403834371),
        ]
        dq = zeros(4)
        Enzyme.autodiff(Enzyme.Reverse, rootsum, Enzyme.Active, Enzyme.Duplicated(q0, dq))
        @test all(isfinite, dq)

        h = 1e-7
        fd = (rootsum(q0 .+ [h, 0, 0, 0]) - rootsum(q0 .- [h, 0, 0, 0])) / 2h
        @test dq[1] ≈ fd rtol = 1e-5
    end

    @testset "Fast-Light Kerr-Schild" begin
        f_r(r) = begin
            x, y, z = Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
                met,
                r,
                θBL,
                ϕBL,
            )
            return x + 2y - 3z
        end

        f_θ(θ) = begin
            x, y, z = Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
                met,
                rBL,
                θ,
                ϕBL,
            )
            return x + 2y - 3z
        end

        f_ϕ(ϕ) = begin
            x, y, z = Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
                met,
                rBL,
                θBL,
                ϕ,
            )
            return x + 2y - 3z
        end

        test_forward(
            f_r,
            Enzyme.DuplicatedNoNeed,
            (rBL, Enzyme.Duplicated);
            rtol = 1e-5,
            atol = 1e-6,
        )
        test_forward(
            f_θ,
            Enzyme.DuplicatedNoNeed,
            (θBL, Enzyme.Duplicated);
            rtol = 1e-5,
            atol = 1e-6,
        )
        test_forward(
            f_ϕ,
            Enzyme.DuplicatedNoNeed,
            (ϕBL, Enzyme.Duplicated);
            rtol = 1e-5,
            atol = 1e-6,
        )
        test_reverse(f_r, Enzyme.Active, (rBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_θ, Enzyme.Active, (θBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_ϕ, Enzyme.Active, (ϕBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
    end
end
