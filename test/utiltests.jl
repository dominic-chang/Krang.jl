@testset "Utility Functions" begin
    met = Kang.Kerr(0.9)
    ηtemp = η(met, 5.0, 5.0, π / 4)
    λtemp = λ(met, 5.0, π / 4)
    @test Kang._pow(-4.0 + 5.0im, 1 / 5) ≈ (-4.0 + 5.0im)^(1 / 5)
    @test abs(Kang._pow(-4.0, 1 / 5)) ≈ abs((-4.0 + 0.0im)^(1 / 5))
    @test Kang._isreal2(1.0 + 0.0im) == true
    @test λtemp == -5.0 * sin(π / 4)
    @test ηtemp == (25.0 - 0.9^2) * cos(π / 4)^2 + 25.0
    @test θ_potential(met, ηtemp, λtemp, π / 4) ≈ β(met, λtemp, ηtemp, π / 4)^2
    @test α(met, λtemp, π / 4) ≈ -λtemp / sin(π / 4)

    ssmet = Kerr(0.0)
    @test r_potential(ssmet, 27.0, 0.0, 5.0) ≈ 5.0 * (5.0 - 3.0) * (5.0 - 3.0) * (5.0 + 6.0)
    @test maximum(abs, Kang.get_radial_roots(ssmet, 27.0, 0.0) .- (-6.0 + 0im, 0.0 + 0im, 3.0 + 0im, 3.0 + 0im)) ≈ 0.0

    @testset "Radial Integrals 1" begin
        a = 0.99
        met = Kang.Kerr(a)
        @testset "Case 2" begin
            α = 10.0
            β = 10.0
            θo = π / 4
            ηcase1 = η(met, α, β, θo)
            λcase1 = λ(met, α, θo)
            roots = get_radial_roots(met, ηcase1, λcase1)
            _, _, _, root = roots
            @test sum(Kang._isreal2.(roots)) == 4

            rs = 1.1 * real(root)
            τ1 = Ir(met, true, θo, rs, α, β)[1]
            I0, I1, I2, Ip, Im = Kang.radial_integrals_case2(met, rs, real.(roots), τ1, true)
            @testset "I0/Ir" begin
                f(r, p) = inv(√(r_potential(met, ηcase1, λcase1, r)))
                prob = IntegralProblem(f, rs, Inf; nout=1)
                sol = solve(prob, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                @test τ1 / sol.u ≈ 1.0 atol = 1e-5
            end
            @testset "I1" begin
                f1(r, p) = r * inv(√(r_potential(met, ηcase1, λcase1, r)))
                prob1 = IntegralProblem(f1, rs, 1e6; nout=1)
                sol1 = solve(prob1, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                @test I1 / (sol1.u - log(1e6)) ≈ 1.0 atol = 1e-5
            end

            @testset "I2" begin
                #Regularized Time            
                f2(r, p) = r^2 * inv(√(r_potential(met, ηcase1, λcase1, r)))
                prob2 = IntegralProblem(f2, rs, 1e6; nout=1)
                sol2 = solve(prob2, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                @test I2 / (sol2.u - 1e6) ≈ 1.0 atol = 1e-3
            end
        end

        @testset "Case 3" begin
            α = 1.0
            β = 3.0
            θo = π / 4
            ηcase3 = η(met, α, β, θo)
            λcase3 = λ(met, α, θo)
            roots = get_radial_roots(met, ηcase3, λcase3)
            @test sum(Kang._isreal2.(roots)) == 2
            rs = 1.1horizon(met)
            τ3 = Ir(met, true, θo, rs, α, β)[1]
            I0, I1, I2, Ip, Im = Kang.radial_integrals_case3(met, rs, roots, τ3)
            @testset "I0/Ir" begin
                f(r, p) = inv(√(r_potential(met, ηcase3, λcase3, r)))
                prob = IntegralProblem(f, rs, Inf; nout=1)
                sol = solve(prob, HCubatureJL(); reltol=1e-5, abstol=1e-5)
                @test τ3 / sol.u ≈ 1.0 atol = 1e-5
            end
            @testset "I1" begin
                f1(r, p) = r * inv(√(r_potential(met, ηcase3, λcase3, r)))
                prob1 = IntegralProblem(f1, rs, 1e6; nout=1)
                sol1 = solve(prob1, HCubatureJL(); reltol=1e-8, abstol=1e-8)
                @test I1 / (sol1.u - log(1e6)) ≈ 1.0 atol = 1e-5
            end
            @testset "I2" begin
                #Regularized Time            
                f2(r, p) = r^2 * inv(√(r_potential(met, ηcase3, λcase3, r)))
                prob2 = IntegralProblem(f2, rs, 1e6; nout=1)
                sol2 = solve(prob2, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                @test I2 / (sol2.u - 1e6) ≈ 1.0 atol = 1e-5
            end
        end
        @testset "Case 4" begin
            α = 0.1
            β = 0.1
            θo = π / 4
            ηcase4 = η(met, α, β, θo)
            λcase4 = λ(met, α, θo)
            roots = get_radial_roots(met, ηcase4, λcase4)
            @test sum(Kang._isreal2.(roots)) == 0
            rs = 1.1horizon(met)
            τ4 = Ir(met, true, θo, rs, α, β)[1]
            I0, I1, I2, Ip, Im = Kang.radial_integrals_case4(met, rs, roots, τ4)
            @testset "I0/Ir" begin
                f(r, p) = inv(√(r_potential(met, ηcase4, λcase4, r)))
                prob = IntegralProblem(f, rs, Inf; nout=1)
                sol = solve(prob, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                @test τ4 / sol.u ≈ 1.0 atol = 1e-5
            end
            @testset "I1" begin
                f1(r, p) = r * inv(√(r_potential(met, ηcase4, λcase4, r)))
                prob1 = IntegralProblem(f1, rs, 1e10; nout=1)
                sol1 = solve(prob1, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                @test I1 / (sol1.u - log(1e10)) ≈ 1.0 atol = 1e-5
            end
            @testset "I2" begin
                #Regularized Time            
                f2(r, p) = r^2 * inv(√(r_potential(met, ηcase4, λcase4, r)))
                prob2 = IntegralProblem(f2, rs, 1e6; nout=1)
                sol2 = solve(prob2, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                @test I2 / (sol2.u - 1e6) ≈ 1.0 atol = 1e-5
            end
        end
    end
    @testset "Radial Integrals 2" begin
        a = 0.99
        met = Kang.Kerr(a)
        @testset "Case 2" begin
            α = 10.0
            β = 10.0
            θo = π / 4
            ηcase1 = η(met, α, β, θo)
            λcase1 = λ(met, α, θo)
            roots = get_radial_roots(met, ηcase1, λcase1)
            _, _, _, root = roots
            @test sum(Kang._isreal2.(roots)) == 4

            rs = 1.1 * real(root)
            τ1 = Ir(met, true, θo, rs, α, β)[1]
            @testset "Ir" begin
                f(r, p) = inv(√(r_potential(met, ηcase1, λcase1, r)))
                prob = IntegralProblem(f, rs, Inf; nout=1)
                sol = solve(prob, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                @test τ1 / sol.u ≈ 1.0 atol = 1e-5
            end

            @testset "Iϕ" begin
                fϕ(r, p) = a * (2r - a * λcase1) * inv((r^2 - 2r + a^2) * √(r_potential(met, ηcase1, λcase1, r)))
                probϕ = IntegralProblem(fϕ, rs, Inf; nout=1)
                solϕ = solve(probϕ, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                Iϕ = Kang.Iϕ_case2(met, real.(roots), real.(Kang._get_root_diffs(roots...)), rs, τ1, true, λcase1)
                @test Iϕ / solϕ.u ≈ 1.0 atol = 1e-5
            end

            @testset "It" begin
                #Regularized Time            
                ft(r, p) = -((r^2 * (r^2 - 2r + a^2) + 2r * (r^2 + a^2 - a * λcase1)) * inv((r^2 - 2r + a^2) * √(r_potential(met, ηcase1, λcase1, r))))

                probt = IntegralProblem(ft, rs, 1e6; nout=1)
                solt = solve(probt, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                It = Kang.It_case2(met, real.(roots), real.(Kang._get_root_diffs(roots...)), rs, τ1, true, λcase1)
                @test It / (solt.u + 1e6 + 2log(1e6)) ≈ 1.0 atol = 1e-3
            end
        end

        @testset "Case 3" begin
            α = 1.0
            β = 3.0
            θo = π / 4
            ηcase3 = η(met, α, β, θo)
            λcase3 = λ(met, α, θo)
            roots = get_radial_roots(met, ηcase3, λcase3)
            @test sum(Kang._isreal2.(roots)) == 2
            rs = 1.1horizon(met)
            τ3 = Ir(met, true, θo, rs, α, β)[1]
            @testset "Ir" begin
                f(r, p) = inv(√(r_potential(met, ηcase3, λcase3, r)))
                prob = IntegralProblem(f, rs, Inf; nout=1)
                sol = solve(prob, HCubatureJL(); reltol=1e-5, abstol=1e-5)
                @test τ3 / sol.u ≈ 1.0 atol = 1e-5
            end
            @testset "Iϕ" begin
                fϕ(r, p) = a * (2r - a * λcase3) * inv((r^2 - 2r + a^2) * √(r_potential(met, ηcase3, λcase3, r)))
                probϕ = IntegralProblem(fϕ, rs, Inf; nout=1)
                solϕ = solve(probϕ, HCubatureJL(); reltol=1e-8, abstol=1e-8)
                Iϕ = Kang.Iϕ_case3(met, roots, Kang._get_root_diffs(roots...), rs, τ3, λcase3)
                @test Iϕ / (solϕ.u) ≈ 1.0 atol = 1e-5
            end
            @testset "It" begin
                #Regularized Time            
                ft(r, p) = -((r^2 * (r^2 - 2r + a^2) + 2r * (r^2 + a^2 - a * λcase3)) * inv((r^2 - 2r + a^2) * √(r_potential(met, ηcase3, λcase3, r))))
                probt = IntegralProblem(ft, rs, 1e6; nout=1)
                solt = solve(probt, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                It = Kang.It_case3(met, roots, Kang._get_root_diffs(roots...), rs, τ3, λcase3)
                #@test It ≈ solt.u + 1e6 + 2log(1e6) atol = 1e-3
                @test It / (solt.u + 1e6 + 2log(1e6)) ≈ 1.0 atol = 1e-3

            end
        end
        @testset "Case 4" begin
            α = 0.05
            β = 0.1
            θo = π / 4
            ηcase4 = η(met, α, β, θo)
            λcase4 = λ(met, α, θo)
            roots = get_radial_roots(met, ηcase4, λcase4)
            @test sum(Kang._isreal2.(roots)) == 0
            rs = 1.1horizon(met)
            τ4 = Ir(met, true, θo, rs, α, β)[1]
            @testset "Ir" begin
                f(r, p) = inv(√(r_potential(met, ηcase4, λcase4, r)))
                prob = IntegralProblem(f, rs, Inf; nout=1)
                sol = solve(prob, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                @test τ4 / sol.u ≈ 1.0 atol = 1e-5
            end
            @testset "Iϕ" begin
                fϕ(r, p) = a * (2r - a * λcase4) * inv((r^2 - 2r + a^2) * √(r_potential(met, ηcase4, λcase4, r)))
                probϕ = IntegralProblem(fϕ, rs, Inf; nout=1)
                solϕ = solve(probϕ, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                Iϕ = Kang.Iϕ_case4(met, roots, Kang._get_root_diffs(roots...), rs, τ4, λcase4)
                @test Iϕ / solϕ.u ≈ 1.0 atol = 2e-2
            end
            @testset "It" begin
                #Regularized Time            
                ft(r, p) = -((r^2 * (r^2 - 2r + a^2) + 2r * (r^2 + a^2 - a * λcase4)) * inv((r^2 - 2r + a^2) * √(r_potential(met, ηcase4, λcase4, r))))
                probt = IntegralProblem(ft, rs, 1e6; nout=1)
                solt = solve(probt, HCubatureJL(); reltol=1e-10, abstol=1e-10)
                It = Kang.It_case4(met, roots, Kang._get_root_diffs(roots...), rs, τ4, λcase4)

                @test It / (solt.u + 1e6 + 2log(1e6)) ≈ 1.0 atol = 1e-2
            end
        end
    end
    @testset "Inclination Integrals" begin
        a = 0.99
        met = Kang.Kerr(a)
        @testset "θo:$θo" for θo in [π / 4, 3π / 4]
            @testset "Ordinary Geodesics" begin
                α = 10.0
                @testset "isindir: $isindir" for isindir in [true, false]
                    β = isindir ? 10.0 : -10.0
                    β *= sign(cos(θo))
                    tempη = η(met, α, β, θo)
                    tempλ = λ(met, α, θo)
                    a2 = met.spin^2
                    Δθ = (1.0 - (tempη + tempλ^2) / a2) / 2
                    Δθ2 = Δθ^2
                    desc = √(Δθ2 + tempη / a2)
                    up = Δθ + desc

                    θturning = acos(√up) * (1 + 1e-10)
                    @testset "Gθ" begin
                        @testset "θs:$θs" for θs in clamp.((1, 3) .* (acos(-√up) - acos(√up)) ./ 4 .+ acos(√up), 0.0, π)
                            fθ(θ, p) = inv(√(Kang.θ_potential(met, tempη, tempλ, θ)))
                            probts = IntegralProblem(fθ, θs, π / 2; nout=1)
                            solθs = solve(probts, HCubatureJL(); reltol=1e-12, abstol=1e-12)
                            probto = IntegralProblem(fθ, θo, π / 2; nout=1)
                            solθo = solve(probto, HCubatureJL(); reltol=1e-12, abstol=1e-12)
                            probt = IntegralProblem(fθ, θturning, π / 2; nout=1)
                            solθt = solve(probt, HCubatureJL(); reltol=1e-12, abstol=1e-12)

                            τ = 0
                            if isindir
                                if sign(cos(θs) * cos(θo)) > 0
                                    τ = abs(2solθt.u - (solθo.u + solθs.u))
                                else
                                    τ = abs(2solθt.u + solθs.u - solθo.u)
                                end
                            else
                                if sign(cos(θs) * cos(θo)) > 0
                                    τ = abs(solθo.u - solθs.u)
                                else
                                    τ = abs(solθo.u + solθs.u)
                                end
                            end
                            tempGθ, _, _, _, _ = Kang.Gθ(met, α, β, θs, θo, isindir, 0)

                            @test tempGθ / τ ≈ 1.0 atol = 1e-3
                            @test Kang.Gs(met, α, β, θo, τ) / ((θs > π / 2 ? -1 : 1) * (solθs.u)) ≈ 1.0 atol = 1e-3
                        end
                    end
                    @testset "Gϕ" begin
                        @testset "θs:$θs" for θs in clamp.((1, 3) .* (acos(-√up) - acos(√up)) ./ 4 .+ acos(√up), 0.0, π)
                            fθ(θ, p) = csc(θ)^2 * inv(√(Kang.θ_potential(met, tempη, tempλ, θ)))
                            probts = IntegralProblem(fθ, θs, π / 2; nout=1)
                            solθs = solve(probts, HCubatureJL(); reltol=1e-12, abstol=1e-12)
                            probto = IntegralProblem(fθ, θo, π / 2; nout=1)
                            solθo = solve(probto, HCubatureJL(); reltol=1e-12, abstol=1e-12)
                            probt = IntegralProblem(fθ, θturning, π / 2; nout=1)
                            solθt = solve(probt, HCubatureJL(); reltol=1e-12, abstol=1e-12)

                            τ = 0
                            if isindir
                                if sign(cos(θs) * cos(θo)) > 0
                                    τ = abs(2solθt.u - (solθo.u + solθs.u))
                                else
                                    τ = abs(2solθt.u + solθs.u - solθo.u)
                                end
                            else
                                if sign(cos(θs) * cos(θo)) > 0
                                    τ = abs(solθo.u - solθs.u)
                                else
                                    τ = abs(solθo.u + solθs.u)
                                end
                            end
                            tempGϕ, _, _, _ = Kang.Gϕ(met, α, β, θs, θo, isindir, 0)

                            @test tempGϕ / τ ≈ 1.0 atol = 1e-3
                        end
                    end
                    @testset "Gt" begin
                        @testset "θs:$θs" for θs in clamp.((1, 3) .* (acos(-√up) - acos(√up)) ./ 4 .+ acos(√up), 0.0, π)
                            ft(θ, p) = cos(θ)^2 * inv(√(Kang.θ_potential(met, tempη, tempλ, θ)))
                            probts = IntegralProblem(ft, θs, π / 2; nout=1)
                            solθs = solve(probts, HCubatureJL(); reltol=1e-12, abstol=1e-12)
                            probto = IntegralProblem(ft, θo, π / 2; nout=1)
                            solθo = solve(probto, HCubatureJL(); reltol=1e-12, abstol=1e-12)
                            probt = IntegralProblem(ft, θturning, π / 2; nout=1)
                            solθt = solve(probt, HCubatureJL(); reltol=1e-12, abstol=1e-12)

                            τ = 0
                            if isindir
                                if sign(cos(θs) * cos(θo)) > 0
                                    τ = abs(2solθt.u - (solθo.u + solθs.u))
                                else
                                    τ = abs(2solθt.u + solθs.u - solθo.u)
                                end
                            else
                                if sign(cos(θs) * cos(θo)) > 0
                                    τ = abs(solθo.u - solθs.u)
                                else
                                    τ = abs(solθo.u + solθs.u)
                                end
                            end
                            tempGt, tempGs, tempGo, _, _ = Kang.Gt(met, α, β, θs, θo, isindir, 0)

                            @test tempGt / τ ≈ 1.0 atol = 1e-3
                        end
                    end
                end
            end
        end
    end
end