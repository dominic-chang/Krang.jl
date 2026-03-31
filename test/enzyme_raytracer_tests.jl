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

        test_forward(f_r, Enzyme.DuplicatedNoNeed, (rBL, Enzyme.Duplicated); rtol = 1e-5, atol = 1e-6)
        test_forward(f_θ, Enzyme.DuplicatedNoNeed, (θBL, Enzyme.Duplicated); rtol = 1e-5, atol = 1e-6)
        test_forward(f_ϕ, Enzyme.DuplicatedNoNeed, (ϕBL, Enzyme.Duplicated); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_r, Enzyme.Active, (rBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_θ, Enzyme.Active, (θBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_ϕ, Enzyme.Active, (ϕBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
    end

    @testset "Fast-Light Kerr-Schild" begin
        f_r(r) = begin
            x, y, z =
                Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
                    met,
                    r,
                    θBL,
                    ϕBL,
                )
            return x + 2y - 3z
        end

        f_θ(θ) = begin
            x, y, z =
                Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
                    met,
                    rBL,
                    θ,
                    ϕBL,
                )
            return x + 2y - 3z
        end

        f_ϕ(ϕ) = begin
            x, y, z =
                Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
                    met,
                    rBL,
                    θBL,
                    ϕ,
                )
            return x + 2y - 3z
        end

        test_forward(f_r, Enzyme.DuplicatedNoNeed, (rBL, Enzyme.Duplicated); rtol = 1e-5, atol = 1e-6)
        test_forward(f_θ, Enzyme.DuplicatedNoNeed, (θBL, Enzyme.Duplicated); rtol = 1e-5, atol = 1e-6)
        test_forward(f_ϕ, Enzyme.DuplicatedNoNeed, (ϕBL, Enzyme.Duplicated); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_r, Enzyme.Active, (rBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_θ, Enzyme.Active, (θBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
        test_reverse(f_ϕ, Enzyme.Active, (ϕBL, Enzyme.Active); rtol = 1e-5, atol = 1e-6)
    end
end
