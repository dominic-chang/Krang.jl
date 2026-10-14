# Miscellaneous functions
##----------------------------------------------------------------------------------------------------------------------
export λ,
    η,
    α,
    β,
    αboundary,
    βboundary,
    r_potential,
    θ_potential,
    get_radial_roots,
    Ir,
    Gθ,
    total_mino_time


"""
Checks if a complex number is real to within √eps of its own magnitude
(or √eps absolute, whichever is larger).
"""
function _isreal2(num)
    T = typeof(real(num))
    return (abs2(imag(num)) <= eps(T) * max(abs2(num), one(T)))
end

"""
Regularized elliptic integral of the third kind
    
# Arguments

- `n`: Parameter
- `ϕ`: Arguments
- `k`: Parameter
"""
regularized_Pi(n, ϕ, k) = regularized_Pi(n, ϕ, k, JacobiElliptic.F(ϕ, k))

"""
`regularized_Pi` with `F(ϕ|k)` supplied by the caller.

Inside one radial-integral case every elliptic integral is evaluated at the same amplitude
and modulus, so `F(ϕ, k)` comes up several times with identical arguments — five times in
case 4, twice in case 3. Threading the value through computes it once. The `_*_at_F`
variants below exist for the same reason.
"""
@inline function regularized_Pi(n, ϕ, k, Fϕ)
    ω = k / n
    return Fϕ - JacobiElliptic.Pi(ω, ϕ, k)
end

regularized_R1(α, φ, j) = regularized_R1(α, φ, j, JacobiElliptic.F(φ, j))

"`regularized_R1` with `F(φ|j)` supplied by the caller; see [`regularized_Pi`](@ref)."
@inline function regularized_R1(α, φ, j, Fφ)
    n = α^2 / (α^2 - 1)
    return 1 / (1 - α^2) * (regularized_Pi(n, φ, j, Fφ))
end

regularized_R2(α, φ, j) = regularized_R2(α, φ, j, JacobiElliptic.CarlsonAlg._F_E(φ, j)...)

"`regularized_R2` with `(F, E)` at `(φ, j)` from the caller; see [`regularized_Pi`](@ref)."
@inline function regularized_R2(α, φ, j, Fφ, Eφ)
    return 1 / (1 - α^2) * (Fφ - α^2 / (j + (1 - j) * α^2) * (
        Eφ
        #-α*sin(φ)*sqrt(1-j*sin(φ)^2)/(1+α*cos(φ)) #Linear Divergent term
    )) #+ 
    #inv(j+(1-j)*α^2)*(2*j-α^2/(α^2-1))*regularized_R1(α, φ, j)
end

regularizedS1(α, φ, j) = regularizedS1(α, φ, j, JacobiElliptic.F(φ, j))

"`regularizedS1` with `F(φ|j)` supplied by the caller; see [`regularized_Pi`](@ref)."
@inline function regularizedS1(α, φ, j, Fφ)
    α2 = α * α
    return inv(1 + α2) * (Fφ + α2 * regularized_Pi(1 + α2, φ, j, Fφ))# - α*f2(α, sin(φ), j)) # logarithmic divergence removed
end

function regularizedS2(α, φ, j)
    Fφ, Eφ = JacobiElliptic.CarlsonAlg._F_E(φ, j)
    return regularizedS2(α, φ, j, regularizedS1(α, φ, j, Fφ), Fφ, Eφ)
end

"""
`regularizedS2` with `regularizedS1(α, φ, j)` supplied by the caller.

`regularizedS1` costs three elliptic evaluations, and every caller of `regularizedS2`
already evaluates it at the same `(α, φ, j)`, so passing it in avoids repeating them.
"""
@inline function regularizedS2(α, φ, j, S1reg, Fφ, Eφ)
    return -inv((1 + α^2) * (1 - j + α^2)) * ((1 - j) * Fφ + α^2 * Eφ) +# + α^2*√(1-j*sin(φ)^2)-α^3) +
           (inv(1 + α^2) + (1 - j) / (1 - j + α^2)) * S1reg
end

function p1(α, j)
    √abs((α^2 - 1) / (j + (1 - j) * α^2))
end

function f1(α, sinφ, j)
    p1temp = p1(α, j)
    tempsinφ = √(1 - j * sinφ^2)
    return p1temp / 2 * log(abs((p1temp * tempsinφ + sinφ) / (p1temp * tempsinφ - sinφ)))
end

function f2(α, sinφ, j)
    p2 = √((1 + α^2) / (1 - j + α^2))
    return p2 / 2 * log(
        abs(
            (1 - p2) / (1 + p2) * (1 + p2 * √(1 - j * sinφ^2)) /
            (1 - p2 * √(1 - j * sinφ^2)),
        ),
    )
end

function R1(α, φ, j)
    #FIXME: This function is undefined when n=1 in Pi(n, ϕ, m) and when α^2 =1
    epsT = eps(α)
    denom = ((1 - α^2) + epsT)
    return (JacobiElliptic.Pi(-α^2 / denom + epsT, φ, j) - α * f1(α, sin(φ), j)) / denom
end

# Combined product `(2·r21·√AB / D) · R1(α, φ, j)` for α = S/D (with
# S = B·r2 + A·r1, D = B·r2 − A·r1, R = 4·A·B·r1·r2 = S² − D²).
# Rewriting via R1's expansion + the n = S²/R, p1 = √(R/(j·D² + (1−j)·S²))
# identities gives the equivalent form `(2·r21·√AB) · (S·f1 − D·Π) / R` which
# is smooth at D → 0, avoiding the catastrophic Float32 cancellation in
# `B·r2 − A·r1` near the photon sphere. R is computed directly as
# `4·A·B·r1·r2` (not via `S² − D²`). `√abs` mirrors `p1`'s defensive form.
#
# `R` is floored at `S²·√eps`. Callers pass `(r1, r2) = (rh - r1, rh - r2)` for a horizon
# `rh`, and `R(rh) = (rh² + a² - aλ)²` is a perfect square, so `rh - r2 → 0⁺` exactly on
# the screen line `λ = (rh² + a²)/a`, quadratically in `λ - λ*`. There `R → 0` while
# `n_α = S²/R → ∞`: `Pi` underflows to zero past `n_α ≈ 1e14` and `R = 0` turns the whole
# expression into `NaN`, which the `(rh - aλ/2)` prefactor at the call site cannot kill
# because it vanishes on the same line (`0 × ∞`, and `0 × NaN` partials under reverse AD).
# The limit itself is finite — the pole at `r = rh` is a principal value whose divergence
# cancels — so bounding `n_α` by `1/√eps` evaluates it at the smallest `R` that `Pi` still
# resolves. That reproduces the true `R → 0` limit to ~7 digits and leaves every
# non-degenerate `R` (where `n_α = O(1)`) untouched. The floor also absorbs the
# `r2 = rh + 1ulp` rounding that would otherwise make `R` negative and the result wrong.
@inline function _R1_alpha_combined(A::T, B::T, r1::T, r2::T, r21::T, φ, j) where {T}
    AB = A * B
    sqrtAB = √AB
    S = B * r2 + A * r1
    D = B * r2 - A * r1
    R = max(4 * AB * r1 * r2, S * S * √eps(T))
    # `n_α = S²/(S² - D²) ≥ 1` exactly, but `R` is computed as `4·A·B·r1·r2`, not as `S² - D²`, so once
    # `D² < eps·S²` (a thin curve on the screen, e.g. `D = B·(r₊ - r2) - A·(r₊ - r1) → 0` in the
    # outer-horizon term of `Iϕ_inf_case3`) `S²/R` can round to just below 1 and `S²/R + eps` lands on
    # exactly 1. With `φ > π/2` the pole of `Π(n_α, φ, j)` at θ = π/2 is then inside the range and `Pi`
    # returns `Inf`, although `D·Π` (and the result) has a finite limit; the `Inf` also turns reverse-mode
    # partials of everything upstream into `NaN`. Clamp to the exact bound before adding `eps`.
    n_α = max(S * S / R, one(T)) + eps(T)
    sinφ = sin(φ)
    tempsinφ = √(one(T) - j * sinφ * sinφ)
    p1_α = √abs(R / (j * D * D + (one(T) - j) * S * S))
    f1_α = (p1_α / 2) * log(abs((p1_α * tempsinφ + sinφ) / (p1_α * tempsinφ - sinφ)))
    Π = JacobiElliptic.Pi(n_α, φ, j)
    return 2 * r21 * sqrtAB * (S * f1_α - D * Π) / R
end

# Horizon term `coef·Π(n, φ, k)` of the four-real-root case, with
# `coef = 2·r43/(√(r31·r42)·rh3·rh4)` and `n = rh3·r41/(rh4·r31)` for `rhi = rh - ri`.
#
# `R(rh) = (rh² + a² - aλ)²` is a perfect square, so `R(rh) = rh1·rh2·rh3·rh4` vanishes —
# and one root collides with the horizon — exactly on the screen line `λ = (rh² + a²)/a`.
# `R'(rh) = ∓(r₊ - r₋)·(η + (λ - a)²)` fixes which root: at `r₋` it is `r4` (or the unused
# `r2`), at `r₊` it is `r3` (or the unused `r1`). So `rh4 → 0` is the collision this term
# has to survive, and there `coef → ∞` while `n → ∞`; the product stays finite because
# `n·Π(n, φ, k) → Δφ·cot φ - F(φ, k) + E(φ, k)`. Numerically both factors overflow
# (`Pi` underflows to zero past `n ≈ 1e14`, and `rh4 = 0` gives `Inf`/`NaN`), which the
# vanishing `(rh - aλ/2)` prefactor at the call site cannot kill — it is zero on the same
# line, so the result is `0 × ∞`, and `0 × NaN` partials under reverse-mode AD.
# Flooring `|rh4|` so that `|n| ≤ 1/√eps` keeps both factors in range and reproduces the
# limit. Near the collision `rh1`, `rh2`, `rh3 > 0`, so `R(rh) ≥ 0` forces `rh4 ≥ 0`: the
# floor pushes onto that side, which also absorbs a `rh4` rounded one ulp past zero. The
# window is `√eps`-relative, so an `rh4` that is legitimately negative and `O(1)` — the
# horizon sitting between `r2` and `r3` — is left untouched.
@inline function _Pi_horizon_combined(r31, r42, r43, r41, rh3, rh4::T, φ, k) where {T}
    floor4 = abs(rh3 * r41 / r31) * √eps(T)
    rh4 = ifelse(abs(rh4) < floor4, floor4, rh4)
    n = rh3 * r41 / (rh4 * r31)
    n += ifelse(n == k, eps(T), zero(T)) # Fix for ∂_m Pi blowing up when n = m
    return 2 / √(r31 * r42) * r43 / (rh3 * rh4) * JacobiElliptic.Pi(n, φ, k)
end

function R2(α, φ, j)
    #FIXME: This function is undefined when α*cos(φ) = 1
    epsT = eps(α)
    denom = ((1 - α^2) + epsT)
    return (
        JacobiElliptic.F(φ, j) -
        α^2 / (j + (1 - j) * α^2) * (
            JacobiElliptic.E(φ, j) -
            α * sin(φ) * √(1 - j * sin(φ)^2) / ((1 + α * cos(φ)) + epsT)
        )
    )/denom + inv(j + (1 - j) * α^2) * (2 * j + α^2 / ((1 - α^2) + epsT)) * R1(α, φ, j)
end

S1(α, φ, j) = S1(α, φ, j, JacobiElliptic.F(φ, j))

"`S1` with `F(φ|j)` supplied by the caller; see [`regularized_Pi`](@ref)."
@inline function S1(α, φ, j, Fφ)
    α2 = α * α
    return inv(1 + α2) * (Fφ + α2 * JacobiElliptic.Pi(1 + α2, φ, j) - α * f2(α, sin(φ), j))
end

S2(α, φ, j) = S2(α, φ, j, JacobiElliptic.CarlsonAlg._F_E(φ, j)...)

"`S2` with `(F, E)` at `(φ, j)` from the caller; see [`regularized_Pi`](@ref)."
@inline function S2(α, φ, j, Fφ, Eφ)
    return -inv((1 + α^2) * (1 - j + α^2)) * (
        (1 - j) * Fφ +
        α^2 * Eφ +
        α^2 * √(1 - j * sin(φ)^2) * (α - tan(φ)) / (1 + α * tan(φ)) - α^3
    ) + (inv(1 + α^2) + (1 - j) / (1 - j + α^2)) * S1(α, φ, j, Fφ)
end

"""
Energy reduced azimuthal angular momentum

# Arguments

- `metric`: Kerr
- `α`: Horizontal Bardeen screen coordinate
- `θo`: Observer inclination
"""
function λ(::Kerr, α, θo)
    return -α * sin(θo)
end

"""
Energy reduced Carter integral

# Arguments

- `metric`: Kerr
- `α`: Horizontal Bardeen screen coordinate
- `β`: Bardeen vertical coordinate
- `θo`: Observer inclination
"""
function η(metric::Kerr, α, β, θo)
    return (α^2 - metric.spin^2) * cos(θo)^2 + β^2
end

"""
Horizontal Bardeen Screen Coordinate

# Arguments

- `metric`: Kerr
- `α`: Horizontal Bardeen screen coordinate
- `θo`: Observer inclination
"""
function α(::Kerr, λ, θo)
    return -λ / sin(θo)
end

"""
Vertical Bardeen Screen Coordinate

# Arguments

- `metric`: Kerr
- `λ`: Energy-reduced azimuthal angular momentum
- `η`: Energy reduced Carter integral 
- `θo`: Observer inclination
"""
function β(metric::Kerr, λ, η, θo)
    temp = η - (α(metric, λ, θo)^2 - metric.spin^2) * cos(θo)^2
    return sqrt(max(temp, zero(temp)))
end

"""
Defines a horizontal boundary on the assymptotic observer's screen that emission that from θs must fall within.

# Arguments

- `metric`: Kerr metric
- `θs`  : Emission Inclination
"""
@inline function αboundary(metric::Kerr, θs)
    return metric.spin * sin(θs)
end

"""
Defines a vertical boundary on the assymptotic observer's screen that emission that from θs must fall within.

# Arguments

- `metric`: Kerr{T} metric
- `α`   : Horizontal Bardeen screen coordinate
- `θo`  : Observer inclination
- `θs`  : Emission Inclination
"""
@inline function βboundary(metric::Kerr, α, θo, θs)
    a = metric.spin
    cosθs2 = cos(θs)^2
    temp = (cos(θo)^2 - cosθs2) * (α^2 - a^2 * (1 - cosθs2)) / (cosθs2 - 1)
    return √max(temp, zero(temp)) #eq 15 DOI 10.3847/1538-4357/acafe3 
end

"""
Radial potential of spacetime

# Arguments

- `metric`: Kerr{T} metric
- `η`  : Reduced Carter constant
- `λ`  : Reduced azimuthal angular momentum
- `r`  : Boyer Lindquist radius
"""
function r_potential(metric::Kerr, η, λ, r)
    a = metric.spin
    λ2 = λ^2
    return a * (a * (r * (r + 2) - η) - 4 * λ * r) +
           r * ((η + η) + (λ2 + λ2) + r * (-η - λ2 + r^2)) # Eq 7 PhysRevD.101.044032
end

"""
Theta potential of a Kerr black hole

# Arguments

- `metric`: Kerr{T} metric
- `η`  : Reduced Carter constant
- `λ`  : Reduced azimuthal angular momentum
- `θ`  : Boyer Lindquist inclination
"""
function θ_potential(metric::Kerr, η, λ, θ)
    a = metric.spin
    return η + a^2 * cos(θ)^2 - λ^2 * cot(θ)^2
end

##----------------------------------------------------------------------------------------------------------------------
# Radial functions
##----------------------------------------------------------------------------------------------------------------------
"""
Returns roots of \$r^4 + (a^2-η-λ^2)r^2 + 2(η+(a-λ)^2)r - a^2η\$

# Arguments

- `metric`: Kerr{T} metric
- `η`  : Reduced Carter constant
- `λ`  : Reduced azimuthal angular momentum
"""
function get_radial_roots(metric::Kerr{T}, η, λ) where {T}
    a = metric.spin

    a2 = a^2
    A = a2 - η - λ * λ
    A2 = A + A
    B = T(2) * (η + (λ - a)^2)
    C = -a2 * η

    P = -A^2 / T(12) - C
    Q = -A / T(3) * (A * A / T(36) + zero(T)im - C) - B^2 / T(8)

    negΔ3 = T(4) * P * P * P + T(27) * (Q^2)
    sq = sqrt(negΔ3 / T(108))
    # Cardano's formula holds with either sign of `sq`: the two choices swap `ωp` and
    # `v = -P/(3ωp)` (their product is `-P/3`), so the three resolvent roots — and with
    # them `ξ0` and the quartic roots — come out the same either way. When `P → 0` with
    # `Q > 0`, `sq → |Q|/2` and the `+` branch cancels: `ωp → 0` turns `v` into `0/0`,
    # making every root NaN, and just outside that band `v` is a ratio of rounding errors
    # that leaves the roots with spurious imaginary parts. Take the `-` branch only where
    # the cancellation actually bites, so every other evaluation stays bit-identical.
    base = -Q / T(2) + sq
    base = ifelse(abs(base) < T(1e-6) * (abs(Q) / 2 + abs(sq)), -Q / T(2) - sq, base)
    ωp = ^(base + zero(T)im, T(1 / 3))

    #C = ((-1+0im)^(2/3), (-1+0im)^(4/3), 1) .* ωp
    C = (-T(1 / 2) + T(√3 / 2)im, -T(1 / 2) - T(√3 / 2)im, one(T) + zero(T)im) .* ωp

    v = -P .* inv.(T(3) .* C)

    ξ0 = argmax(real, (C .+ v)) - A / T(3)
    ξ02 = 2ξ0

    predet1 = A2 + ξ02
    predet2 = (√T(2) * B) * inv(sqrt(ξ0))
    det1 = sqrt(-(predet1 - predet2))
    det2 = sqrt(-(predet1 + predet2))

    sqrtξ02 = sqrt(ξ02)

    r1 = (-sqrtξ02 - det1) / 2
    r2 = (-sqrtξ02 + det1) / 2
    r3 = (sqrtξ02 - det2) / 2
    r4 = (sqrtξ02 + det2) / 2

    # All case-3 formulas assume that the two real roots precede the complex
    # conjugate pair.  Near a branch boundary the quartic formula can instead
    # return them in slots 1 and 4, which makes products such as
    # (r3 - r2) * (r4 - r2) spuriously negative after taking `real`.
    roots = (r1, r2, r3, r4)
    if _isreal2(r1) && _isreal2(r4) && !_isreal2(r2) && !_isreal2(r3)
        # Keep the conventional orientation of the conjugate pair as well:
        # positive imaginary part in slot 3 and negative in slot 4.
        roots = (r1, r4, r3, r2)
    end
    return roots
end

@inline function _get_root_diffs(r1, r2, r3, r4)
    return r2 - r1, r3 - r1, r3 - r2, r4 - r1, r4 - r2, r4 - r3
end

"""
Mino time of trajectory between an observer at infinity and point at radius rs

# Arguments 

- `pix` : Pixel information
- `rs` : Emission radius
- `νr` : Sign of the radial velocity at emission.
"""
function mino_time(pix, rs, νr)
    return Ir(pix, νr, rs)[1]
end

"""
Returns the antiderivative \$I_r=\\int\\frac{dr}{\\sqrt{\\mathcal{R(r)}}}\$ evaluated at infinity.
See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `roots`  : Roots of the radial potential
"""
function Ir_inf(metric::Kerr{T}, roots) where {T}
    #root_diffs = _get_root_diffs(roots...)
    numreals = sum(map(_isreal2, roots))

    if numreals == 4 #case 2
        return Ir_inf_case1_and_2(metric, map(real, roots))
    elseif numreals == 2 #case3
        return Ir_inf_case3(metric, roots)
    else #case 4
        return Ir_inf_case4(metric, roots)
    end
end

function Ir_inf_case1_and_2(::Kerr, roots::NTuple{4})
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    k = r32 * r41 / (r31 * r42)
    coef = 2 / √real(r31 * r42)
    Ir_inf = coef * JacobiElliptic.F(asin(√(r31 / r41)), k)

    return Ir_inf
end

function Ir_inf_case3(::Kerr, roots::NTuple{4})
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    r1, r2, r21 = real.((r1, r2, r21))

    A2 = abs2(r32)
    B2 = abs2(r31)
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    coef = 1 * √inv(A * B)
    return coef * JacobiElliptic.F((acos(clamp((A - B) / (A + B), -1, 1))), k3)

end

function Ir_inf_case4(::Kerr{T}, roots::NTuple{4}) where {T}
    _, r2, _, r4 = roots

    a1 = abs(imag(r4))
    a2 = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1-a2)^2 + (b1-b2)^2)
    D = sqrt((a1+a2)^2 + (b1-b2)^2)
    k4 = 4C * D / (C + D)^2


    k4 = T(4) * C * D / (C + D)^2

    go = √max((T(4)a2^2 - (C - D)^2) / ((C + D)^2 - T(4)a2^2), zero(T))

    return 2 / (C + D) * JacobiElliptic.F(T(π / 2) + atan(go), k4)
end

"""
Returns the antiderivative \$I_r=\\int\\frac{dr}{\\sqrt{\\mathcal{R(r)}}}\$.
See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `rs` : Emission radius
- `roots`  : Roots of the radial potential
- `νr` : Radial emission direction (Only necessary for case 1&2 geodesics)
"""
function Ir_s(metric::Kerr{T}, rs, roots, νr) where {T}
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return Ir_s_case1_and_2(metric, rs, real.(roots), νr)
    elseif numreals == 2 #case3
        return Ir_s_case3(metric, rs, roots)
    else #case 4
        return Ir_s_case4(metric, rs, roots)
    end
end

function Ir_s_case1_and_2(::Kerr{T}, rs, roots::NTuple{4}, νr) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    k = r32 * r41 / (r31 * r42)
    x2_s = √abs((rs - r4) / (rs - r3) * r31 / r41)
    coef = 2 / √real(r31 * r42)
    Ir_s = (x2_s > 1) ? T(Inf) : coef * JacobiElliptic.F(asin(x2_s), k)

    return -(-1)^νr * Ir_s
end

function Ir_s_case3(::Kerr, rs, roots::NTuple{4})
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    r1, r2, r21 = real.((r1, r2, r21))

    A2 = abs2(r32)
    B2 = abs2(r31)
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    temprat = B * (rs - r2) / (A * (rs - r1))
    x3_s = clamp(((1 - temprat) / (1 + temprat)), -1, 1)
    coef = 1 * √inv(A * B)
    Ir_s = coef * JacobiElliptic.F((acos(x3_s)), k3)

    return Ir_s
end

function Ir_s_case4(::Kerr{T}, rs, roots::NTuple{4}) where {T}
    _, r2, _, r4 = roots

    a1 = abs(imag(r4))
    a2 = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1-a2)^2 + (b1-b2)^2)
    D = sqrt((a1+a2)^2 + (b1-b2)^2)
    k4 = 4C * D / (C + D)^2


    k4 = T(4) * C * D / (C + D)^2

    go = √max((T(4)a2^2 - (C - D)^2) / ((C + D)^2 - T(4)a2^2), zero(T))
    x4_s = (rs + b1) / a2
    coef = 2 / (C + D)
    Ir_s = coef * JacobiElliptic.F(atan(x4_s) + atan(go), k4)
    return Ir_s
end

"""
Returns the antiderivative \$I_ϕ=\\int\\frac{a(2Mr-a\\lambda)}{\\sqrt{\\Delta\\mathcal{R(r)}}}dr\$ evaluated at infinity
without I0 terms.

See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `roots` : Radial roots
- `λ`  : Reduced azimuthal angular momentum
"""
function Iϕ_inf(metric::Kerr{T}, roots, λ) where {T}
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return Iϕ_inf_case2(metric, real.(roots), λ)
    elseif numreals == 2 #case3
        return Iϕ_inf_case3(metric, roots, λ)
    else #case4
        return Iϕ_inf_case4(metric, roots, λ)
    end
end

function Iϕ_inf_case2(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, r43 = _get_root_diffs(roots...)
    a = metric.spin
    a2 = a * a
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)
    rp3 = rp - r3
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4
    #FIXME: This is a hack to avoid the division by zero
    if (rp3 == zero(T))
        rp3 = eps(T)
    end

    k = r32 * r41 / (r31 * r42)
    x2_o = √(r31 / r41)

    Πp_o = _Pi_horizon_combined(r31, r42, r43, r41, rp3, rp4, asin(x2_o), k)
    Πm_o = _Pi_horizon_combined(r31, r42, r43, r41, rm3, rm4, asin(x2_o), k)

    Ipo_inf_m_I0_terms = -Πp_o
    Imo_inf_m_I0_terms = -Πm_o

    return 2a / (rp - rm) *
           ((rp - a * λ / 2) * Ipo_inf_m_I0_terms - (rm - a * λ / 2) * Imo_inf_m_I0_terms)
end

function Iϕ_inf_case3(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r1, r2, r21 = real.((r1, r2, r21))

    a = metric.spin
    a2 = a * a
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = abs2(r32)
    B2 = abs2(r31)
    if A2 < zero(T) || B2 < zero(T)
        return T(Inf)
    end
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)

    x3_o = clamp((A - B) / (A + B), -1, 1)
    φ_o = acos(x3_o)

    (isnan(φ_o)) && return T(NaN)

    C1p_o = _R1_alpha_combined(A, B, rp1, rp2, r21, φ_o, k3)
    C1m_o = _R1_alpha_combined(A, B, rm1, rm2, r21, φ_o, k3)

    Ip = -inv(B * rp2 + A * rp1) * C1p_o
    Im = -inv(B * rm2 + A * rm1) * C1m_o

    return (2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im))
end

function Iϕ_inf_case4(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    _, r2, _, r4 = roots
    a = metric.spin
    a2 = a * a
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)

    a1 = abs(imag(r4))
    a2 = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1-a2)^2 + (b1-b2)^2)
    D = sqrt((a1+a2)^2 + (b1-b2)^2)
    k4 = 4C * D / (C + D)^2

    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - 1) / (go + x4_p)
    gm = (go * x4_m - 1) / (go + x4_m)

    (isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    # Every elliptic integral in this case is evaluated at this one amplitude, and five of
    # them are `F` with identical arguments, so evaluate it once and thread it through.
    φ4_o = T(π / 2) + atan(go)
    F4_o, E4_o = JacobiElliptic.CarlsonAlg._F_E(φ4_o, k4)
    S1p_o = S1(gp, φ4_o, k4, F4_o)
    S1m_o = S1(gm, φ4_o, k4, F4_o)

    Ip =
        go / (a2 * (1 - go * x4_p)) *
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1p_o))
    Im =
        go / (a2 * (1 - go * x4_m)) *
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_m))) * (S1m_o))

    return 2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
end

"""
Returns the antiderivative \$I_ϕ=\\int\\frac{a(2Mr-a\\lambda)}{\\sqrt{\\Delta\\mathcal{R(r)}}}dr\$ with full I0 terms.

See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `τ`: Minotime 
- `roots` : Radial roots
- `νr` : Radial emission direction (Only necessary for case 1&2 geodesics)
- `λ`  : Reduced azimuthal angular momentum
"""
function Iϕ_w_I0_terms(metric::Kerr{T}, rs, τ, roots, νr, λ) where {T}
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return Iϕ_w_I0_terms_case2(metric, rs, τ, real.(roots), νr, λ)
    elseif numreals == 2 #case3
        return Iϕ_w_I0_terms_case3(metric, rs, τ, roots, λ)
    else #case4
        return Iϕ_w_I0_terms_case4(metric, rs, τ, roots, λ)
    end
end

function Iϕ_w_I0_terms_case2(metric::Kerr{T}, rs, τ, roots::NTuple{4}, νr, λ) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, r43 = _get_root_diffs(roots...)
    a = metric.spin
    a2 = a * a
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)
    rp3 = rp - r3
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4
    #FIXME: This is a hack to avoid the division by zero
    if (rp3 == zero(T))
        rp3 = eps(T)
    end

    k = r32 * r41 / (r31 * r42)
    x2_s = √max((rs - r4) / (rs - r3) * r31 / r41, zero(T))
    #!(-1 < x2_s < 1) 
    !(x2_s < 1) && return T(Inf)

    Πp_s = _Pi_horizon_combined(r31, r42, r43, r41, rp3, rp4, asin(x2_s), k)
    Πm_s = _Pi_horizon_combined(r31, r42, r43, r41, rm3, rm4, asin(x2_s), k)

    Ip = -τ / rp3 - (-1)^νr * Πp_s
    Im = -τ / rm3 - (-1)^νr * Πm_s

    return -2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
end

function Iϕ_w_I0_terms_case3(metric::Kerr{T}, rs, τ, roots::NTuple{4}, λ) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r1, r2, r21 = real.((r1, r2, r21))

    a = metric.spin
    a2 = a * a
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = abs2(r32)
    B2 = abs2(r31)
    if A2 < zero(T) || B2 < zero(T)
        return T(Inf)
    end
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)

    temprat = B * (rs - r2) * real(inv(A * (rs - r1)))
    x3_s = clamp(((1 - temprat) * real(inv(1 + temprat))), -1, 1)
    φ_s = acos(x3_s)

    (isnan(φ_s)) && return T(NaN)

    C1p_s = _R1_alpha_combined(A, B, rp1, rp2, r21, φ_s, k3)
    C1m_s = _R1_alpha_combined(A, B, rm1, rm2, r21, φ_s, k3)

    Ip = -inv(B * rp2 + A * rp1) * ((B + A) * τ - C1p_s)
    Im = -inv(B * rm2 + A * rm1) * ((B + A) * τ - C1m_s)

    return -(2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im))
end

function Iϕ_w_I0_terms_case4(metric::Kerr{T}, rs, τ, roots::NTuple{4}, λ) where {T}
    _, r2, _, r4 = roots
    a = metric.spin
    a2 = a * a
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)

    a1 = abs(imag(r4))
    a2 = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1-a2)^2 + (b1-b2)^2)
    D = sqrt((a1+a2)^2 + (b1-b2)^2)
    k4 = 4C * D / (C + D)^2

    x4_s = (rs + b1) / a2
    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - 1) / (go + x4_p)
    gm = (go * x4_m - 1) / (go + x4_m)

    (isnan(x4_s) || isnan(x4_p) || isnan(x4_m)) && return T(Inf)
    # One amplitude for every elliptic integral here; `F` at it is shared (see
    # [`regularized_Pi`](@ref)).
    φ4_s = atan(x4_s) + atan(go)
    F4_s, E4_s = JacobiElliptic.CarlsonAlg._F_E(φ4_s, k4)
    S1p_s = S1(gp, φ4_s, k4, F4_s)
    S1m_s = S1(gm, φ4_s, k4, F4_s)

    Ip =
        go / (a2 * (1 - go * x4_p)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (-S1p_s))
    Im =
        go / (a2 * (1 - go * x4_m)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_m))) * (-S1m_s))

    return -2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
end

"""
Returns the antiderivative \$I_t=\\int\\frac{r^2\\Delta+2Mr(r^2+a^2-a\\lambda)}{\\sqrt{\\Delta\\mathcal{R(r)}}}dr\$ 
evaluated at infinity without I0 terms.

See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `roots` : Radial roots
- `λ`  : Reduced azimuthal angular momentum
"""
function It_inf(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return It_inf_case2(metric, real.(roots), λ)
    elseif numreals == 2 #case3
        return It_inf_case3(metric, roots, λ)
    else #case4
        return It_inf_case4(metric, roots, λ)
    end
end

function It_inf_case2(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r43 = r4 - r3
    a = metric.spin
    rp = 1 + √(1 - a^2)
    rm = 1 - √(1 - a^2)
    rp3 = rp - r3
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4
    #FIXME: This is a hack to avoid the division by zero
    if (rp3 == zero(T))
        rp3 = eps(T)
    end

    k = r32 * r41 / (r31 * r42)
    x2_o = √abs(r31 / r41)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    # `inv(√n)` and `x2_o` are the same turning point written two ways and differed by a
    # few ulp; use the one amplitude, and take `F` and `E` from the single Carlson `RF`
    # pass they share (`_F_E` falls back to the separate entry points off the principal
    # arch, so it is safe for any argument).
    φ2_o = asin(x2_o)
    F2_o, E2_o = JacobiElliptic.CarlsonAlg._F_E(φ2_o, k)
    E_o = √(r31 * r42) * E2_o

    #equation B37
    I1_total = log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, φ2_o, k, F2_o))# Removed the logarithmic divergence
    #equation B38
    I2_total = r3 - E_o # asymptotic Divergent piece is not included

    Πp_o = _Pi_horizon_combined(r31, r42, r43, r41, rp3, rp4, asin(x2_o), k)
    Πm_o = _Pi_horizon_combined(r31, r42, r43, r41, rm3, rm4, asin(x2_o), k)

    Ip_total = -Πp_o
    Im_total = -Πm_o

    return -(
        4 / (rp - rm) *
        (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
        2 * I1_total +
        I2_total
    )
end

function It_inf_case3(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r21 = real(r21)
    r2 = real(r2)
    r1 = real(r1)
    a = metric.spin
    rp = 1 + √(1 - a^2)
    rm = 1 - √(1 - a^2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = abs2(r32)
    B2 = abs2(r31)
    if A2 < zero(T) || B2 < zero(T)
        return T(Inf)
    end
    A, B = √A2, √B2

    k3 = real(((A + B)^2 - r21^2) / (4 * A * B))

    x3_o = min((A - B) / (A + B), 1)
    φ_o = acos(x3_o)

    (isnan(φ_o)) && return T(NaN)

    αo = (B + A) / (B - A)

    # `regularized_R1` and `regularized_R2` both want `F(φ_o, k3)`, and `R2` also wants
    # `E` at that amplitude; one `_F_E` pass covers all three.
    F3_o, E3_o = JacobiElliptic.CarlsonAlg._F_E(φ_o, k3)
    Π1_o = 2 * r21 * √real(A * B) / (B2 - A2) * regularized_R1(αo, φ_o, k3, F3_o) # Divergence is removed, will be added back in the end
    Π2_o = ((2 * r21 * √(A * B) / (B2 - A2))^2) * regularized_R2(αo, φ_o, k3, F3_o, E3_o)# Divergence is removed, will be added back in the end

    C1p_o = _R1_alpha_combined(A, B, rp1, rp2, r21, φ_o, k3)
    C1m_o = _R1_alpha_combined(A, B, rm1, rm2, r21, φ_o, k3)

    # Removed logarithmic divergence
    I1_total = Π1_o + log(16 * r21^2 / ((A2 - B2)^2 + 4 * A * B * r21^2)) / 2
    # Removed linear divergence
    I2_total = (-√(A * B) * (Π2_o)) + (B * r2 + A * r1) / (A + B)
    Ip_total = -inv(B * rp2 + A * rp1) * C1p_o
    Im_total = -inv(B * rm2 + A * rm1) * C1m_o

    return -(
        4 / (rp - rm) *
        (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
        2 * I1_total +
        I2_total
    )
end

function It_inf_case4(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    oneT = one(T)
    a = metric.spin
    _, r2, _, r4 = roots

    rp = one(T) + √(oneT - a^2)
    rm = one(T) - √(oneT - a^2)

    a1 = abs(imag(r4))
    a2 = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1-a2)^2 + (b1-b2)^2)
    D = sqrt((a1+a2)^2 + (b1-b2)^2)
    k4 = 4C * D / (C + D)^2

    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - 1) / (go + x4_p)
    gm = (go * x4_m - 1) / (go + x4_m)

    (isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    # Every elliptic integral in this case is evaluated at this one amplitude, and five of
    # them are `F` with identical arguments, so evaluate it once and thread it through.
    φ4_o = T(π / 2) + atan(go)
    F4_o, E4_o = JacobiElliptic.CarlsonAlg._F_E(φ4_o, k4)
    S1p_o = S1(gp, φ4_o, k4, F4_o)
    S1m_o = S1(gm, φ4_o, k4, F4_o)

    S1reg_o = regularizedS1(go, φ4_o, k4, F4_o)
    # Divergences are removed here and added back at the end.
    Π1_o = 2 / (C + D) * (a2 / go * (1 + go^2)) * S1reg_o
    Π2_o =
        2 / (C + D) *
        (a2 / go * (1 + go^2))^2 *
        regularizedS2(go, φ4_o, k4, S1reg_o, F4_o, E4_o)

    I1_total =
        -(Π1_o) +
        log(
            (
                T(16) *
                (oneT + go^2 - sqrt((oneT + go^2) * (oneT + go^2 - k4))) *
                (oneT + go^2 - k4)
            ) / (
                (C + D)^2 *
                ((oneT + go^2)^2 - k4) *
                k4 *
                (oneT + sqrt(oneT - k4 / (oneT + go^2)))
            ),
        ) / T(2)

    # Removed linear divergence
    I2_total =
        -2(a2 / go - b1) * (Π1_o) + (Π2_o) - (
            (
                16 * a2^4 + (C^2 - D^2)^2 - 8 * (a2^2) * (C^2 + D^2) +
                8 * (a2^3) * (C + D - 2 * b1 * go) +
                2 * a2 * (C + D) * (-(C - D)^2 + 2 * b1 * (C + D) * go)
            ) / (4 * a2 * (4 * a2^2 - (C + D)^2) * go)
        )
    Ip_total =
        go / (a2 * (1 - go * x4_p)) *
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1p_o))
    Im_total =
        go / (a2 * (1 - go * x4_m)) *
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_m))) * (S1m_o))

    return -(
        T(4) / (rp - rm) *
        (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
        2 * I1_total +
        I2_total
    )
end

"""
Returns the antiderivative \$I_t=\\int\\frac{r^2\\Delta+2Mr(r^2+a^2-a\\lambda)}{\\sqrt{\\Delta\\mathcal{R(r)}}}dr\$ with
 I0 terms.

See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `roots` : Radial roots
- `λ`  : Reduced azimuthal angular momentum
- `νr` : Radial emission direction (Only necessary for case 1&2 geodesics)
"""
@inline function It_w_I0_terms(metric::Kerr{T}, rs, τ, roots::NTuple{4}, λ, νr) where {T}
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return It_w_I0_terms_case2(metric, rs, τ, real.(roots), λ, νr)
    elseif numreals == 2 #case3
        return It_w_I0_terms_case3(metric, rs, τ, roots, λ)
    else #case4
        return It_w_I0_terms_case4(metric, rs, τ, roots, λ)
    end
end

@inline function It_w_I0_terms_case2(
    metric::Kerr{T},
    rs,
    τ,
    roots::NTuple{4},
    λ,
    νr,
) where {T}
    r1, r2, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r43 = r4 - r3
    a = metric.spin
    rp = 1 + √(1 - a^2)
    rm = 1 - √(1 - a^2)
    rp3 = rp - r3
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4
    #FIXME: This is a hack to avoid the division by zero
    if (rp3 == zero(T))
        rp3 = eps(T)
    end

    k = r32 * r41 / (r31 * r42)
    x2_s = √abs((rs - r4) / (rs - r3) * r31 / r41)
    #!(-1 < x2_s < 1)
    !(x2_s < 1) && return T(Inf)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    E_s = √(r31 * r42) * JacobiElliptic.E(asin(x2_s), k)
    Π1_s = coef * JacobiElliptic.Pi(n, asin(x2_s), k)

    I0_total = τ

    poly_coefs = (
        r1 * r2 * r3 * r4,
        -r2 * r3 * r4 + r1 * (r2 * (-r3 - r4) - r3 * r4),
        r3 * r4 + r2 * (r3 + r4) + r1 * (r2 + r3 + r4),
        -r1 - r2 - r3 - r4,
        1,
    )

    #equation B37
    I1_total = r3 * I0_total + r43 * (-1)^νr * Π1_s# Removed the logarithmic divergence
    #equation B38
    I2_s = √(evalpoly(rs, poly_coefs)) / (rs - r3) - E_s
    I2_total = -(r1 * r4 + r2 * r3) / 2 * τ + (-1)^νr * I2_s# asymptotic Divergent piece is not included

    Πp_s = _Pi_horizon_combined(r31, r42, r43, r41, rp3, rp4, asin(x2_s), k)
    Πm_s = _Pi_horizon_combined(r31, r42, r43, r41, rm3, rm4, asin(x2_s), k)

    Ip_total = -τ / rp3 - (-1)^νr * Πp_s
    Im_total = -τ / rm3 - (-1)^νr * Πm_s

    return (
        4 / (rp - rm) *
        (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
        4 * I0_total +
        2 * I1_total +
        I2_total
    )
end

@inline function It_w_I0_terms_case3(metric::Kerr{T}, rs, τ, roots::NTuple{4}, λ) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r21 = real(r21)
    r2 = real(r2)
    r1 = real(r1)
    a = metric.spin
    rp = 1 + √(1 - a^2)
    rm = 1 - √(1 - a^2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = abs2(r32)
    B2 = abs2(r31)
    if A2 < zero(T) || B2 < zero(T)
        return T(Inf)
    end
    A, B = √A2, √B2

    k3 = real(((A + B)^2 - r21^2) / (4 * A * B))

    temprat = B * (rs - r2) * inv(A * (rs - r1))
    x3_s = clamp(real(((1 - temprat) * inv(1 + temprat))), -1, 1)
    φ_s = acos(x3_s)

    (isnan(φ_s)) && return T(NaN)

    αo = (B + A) / (B - A)

    Π1_s = 2 * r21 * √real(A * B) / (B2 - A2) * R1(αo, φ_s, k3)
    Π2_s = ((2 * r21 * √(A * B) / (B2 - A2))^2) * R2(αo, φ_s, k3)

    C1p_s = _R1_alpha_combined(A, B, rp1, rp2, r21, φ_s, k3)
    C1m_s = _R1_alpha_combined(A, B, rm1, rm2, r21, φ_s, k3)

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = (B * r2 + A * r1) / (B + A) * I0_total - Π1_s
    # Removed linear divergence
    I2_total = (
        (((B * r2 + A * r1) / (B + A))^2) * I0_total -
        2 * (B * r2 + A * r1) / (B + A) * (-Π1_s) - √(A * B) * (-Π2_s)
    ) #+ (B * r2 + A * r1) / (A + B)
    Ip_total = -inv(B * rp2 + A * rp1) * ((B + A) * I0_total - C1p_s)
    Im_total = -inv(B * rm2 + A * rm1) * ((B + A) * I0_total - C1m_s)

    return (
        4 / (rp - rm) *
        (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
        4 * I0_total +
        2 * I1_total +
        I2_total
    )
end

@inline function It_w_I0_terms_case4(metric::Kerr{T}, rs, τ, roots::NTuple{4}, λ) where {T}
    a = metric.spin
    _, r2, _, r4 = roots
    rp = 1 + √(1 - a^2)
    rm = 1 - √(1 - a^2)

    a1 = abs(imag(r4))
    a2 = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1-a2)^2 + (b1-b2)^2)
    D = sqrt((a1+a2)^2 + (b1-b2)^2)
    k4 = 4C * D / (C + D)^2

    x4_s = (rs + b1) / a2
    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - 1) / (go + x4_p)
    gm = (go * x4_m - 1) / (go + x4_m)

    (isnan(x4_s) || isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    # One amplitude for every elliptic integral here; `F` at it is shared (see
    # [`regularized_Pi`](@ref)).
    φ4_s = atan(x4_s) + atan(go)
    F4_s, E4_s = JacobiElliptic.CarlsonAlg._F_E(φ4_s, k4)
    S1p_s = S1(gp, φ4_s, k4, F4_s)
    S1m_s = S1(gm, φ4_s, k4, F4_s)

    Π1_s = 2 / (C + D) * (a2 / go * (1 + go^2)) * S1(go, φ4_s, k4, F4_s)
    Π2_s = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * S2(go, φ4_s, k4, F4_s, E4_s)

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = (a2 / go - b1) * τ - (-Π1_s)

    # Removed linear divergence
    I2_total = ((a2 / go - b1)^2) * τ - 2(a2 / go - b1) * (-Π1_s) + (-Π2_s)
    Ip_total =
        go / (a2 * (1 - go * x4_p)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (-S1p_s))
    Im_total =
        go / (a2 * (1 - go * x4_m)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_m))) * (-S1m_s))

    return (
        4 / (rp - rm) *
        (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
        4 * I0_total +
        2 * I1_total +
        I2_total
    )# + (logdiv + lineardiv)
end

@inline function radial_inf_integrals(met::Kerr{T}, roots::NTuple{4}) where {T}
    numreals = sum(_isreal2.(roots))
    if numreals == 4
        I1, I2, Ip, Im = Krang.radial_inf_integrals_case2(met, roots)
    elseif numreals == 2
        I1, I2, Ip, Im = Krang.radial_inf_integrals_case3(met, roots)
    else
        I1, I2, Ip, Im = Krang.radial_inf_integrals_case4(met, roots)
    end
    return I1, I2, Ip, Im
end

"""
Returns the radial integrals for the case where there are four real roots in the radial potential, with roots outside the horizon.
"""
@inline function radial_inf_integrals_case2(metric::Kerr{T}, roots::NTuple{4}) where {T}
    roots = real.(roots)
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r43 = r4 - r3
    a = metric.spin
    rp = 1 + √(1 - a^2)
    rm = 1 - √(1 - a^2)
    rp3 = rp - r3
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4
    #FIXME: This is a hack to avoid the division by zero
    if (rp3 == zero(T))
        rp3 = eps(T)
    end

    k = r32 * r41 / (r31 * r42)
    x2_o = √abs(r31 / r41)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    # `inv(√n)` and `x2_o` are the same turning point written two ways and differed by a
    # few ulp; use the one amplitude, and take `F` and `E` from the single Carlson `RF`
    # pass they share (`_F_E` falls back to the separate entry points off the principal
    # arch, so it is safe for any argument).
    φ2_o = asin(x2_o)
    F2_o, E2_o = JacobiElliptic.CarlsonAlg._F_E(φ2_o, k)
    E_o = √(r31 * r42) * E2_o

    #equation B37
    I1o_m_I0_terms =
        log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, φ2_o, k, F2_o))
    #equation B38
    I2o_m_I0_terms = r3 - E_o# asymptotic Divergent piece is not included

    Ipo_m_I0_terms = -_Pi_horizon_combined(r31, r42, r43, r41, rp3, rp4, asin(x2_o), k)
    Imo_m_I0_terms = -_Pi_horizon_combined(r31, r42, r43, r41, rm3, rm4, asin(x2_o), k)

    return I1o_m_I0_terms, I2o_m_I0_terms, Ipo_m_I0_terms, Imo_m_I0_terms

end

"""
Returns the radial integrals for the case where there are two real roots in the radial potential
"""
@inline function radial_inf_integrals_case3(metric::Kerr{T}, roots::NTuple{4}) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    r1, r2, r21 = real.((r1, r2, r21))
    a = metric.spin
    a2 = a * a
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)
    rp1 = rp - r1
    rp2 = rp - r2
    rm1 = rm - r1
    rm2 = rm - r2

    A2 = abs2(r32)
    B2 = abs2(r31)
    A, B = √A2, √B2
    k3 = ((A + B)^2 - r21^2) / (4 * A * B)

    αo = (B + A) / (B - A)

    φ_o = acos((A - B) / (A + B))
    F3_o, E3_o = JacobiElliptic.CarlsonAlg._F_E(φ_o, k3)
    Π1_o = 2 * r21 * √real(A * B) / (B2 - A2) * regularized_R1(αo, φ_o, k3, F3_o)
    Π2_o = ((2 * r21 * √(A * B) / (B2 - A2))^2) * regularized_R2(αo, φ_o, k3, F3_o, E3_o)

    C1p_o = _R1_alpha_combined(A, B, rp1, rp2, r21, φ_o, k3)
    C1m_o = _R1_alpha_combined(A, B, rm1, rm2, r21, φ_o, k3)

    I1o_m_I0_terms = Π1_o + log(16 * r21^2 / ((A2 - B2)^2 + 4 * A * B * r21^2)) / 2
    # Removed linear divergence
    I2o_m_I0_terms = -√(A * B) * Π2_o + (B * r2 + A * r1) / (A + B)
    Ipo_m_I0_terms = -inv(B * rp2 + A * rp1) * C1p_o
    Imo_m_I0_terms = -inv(B * rm2 + A * rm1) * C1m_o

    return I1o_m_I0_terms, I2o_m_I0_terms, Ipo_m_I0_terms, Imo_m_I0_terms
end

"""
Returns the radial integrals for the case where there are no real roots in the radial potential
"""
@inline function radial_inf_integrals_case4(metric::Kerr{T}, roots::NTuple{4}) where {T}
    _, r2, _, r4 = roots
    a2 = metric.spin^2
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)

    a1 = abs(imag(r4))
    a2 = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1-a2)^2 + (b1-b2)^2)
    D = sqrt((a1+a2)^2 + (b1-b2)^2)
    k4 = 4C * D / (C + D)^2

    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - 1) / (go + x4_p)
    gm = (go * x4_m - 1) / (go + x4_m)

    # Every elliptic integral in this case is evaluated at this one amplitude, and five of
    # them are `F` with identical arguments, so evaluate it once and thread it through.
    φ4_o = T(π / 2) + atan(go)
    F4_o, E4_o = JacobiElliptic.CarlsonAlg._F_E(φ4_o, k4)
    S1p_o = S1(gp, φ4_o, k4, F4_o)
    S1m_o = S1(gm, φ4_o, k4, F4_o)
    S1reg_o = regularizedS1(go, φ4_o, k4, F4_o)
    # Divergences are removed here and added back at the end.
    Π1_o = 2 / (C + D) * (a2 / go * (1 + go^2)) * S1reg_o
    Π2_o =
        2 / (C + D) *
        (a2 / go * (1 + go^2))^2 *
        regularizedS2(go, φ4_o, k4, S1reg_o, F4_o, E4_o)

    # Removed logarithmic divergence
    I1o_num = (16 * (1 + go^2 - sqrt((1 + go^2) * (1 + go^2 - k4))) * (1 + go^2 - k4))
    I1o_den = ((C + D)^2 * ((1 + go^2)^2 - k4) * k4 * (1 + sqrt(1 - k4 / (1 + go^2))))
    I1o_m_I0_terms = -Π1_o + log(I1o_num/I1o_den) / 2

    # Removed linear divergence
    I2o_m_I0_terms =
        -2(a2 / go - b1) * Π1_o + Π2_o - (
            (
                16 * a2^4 + (C^2 - D^2)^2 - 8 * (a2^2) * (C^2 + D^2) +
                8 * (a2^3) * (C + D - 2 * b1 * go) +
                2 * a2 * (C + D) * (-(C - D)^2 + 2 * b1 * (C + D) * go)
            ) / (4 * a2 * (4 * a2^2 - (C + D)^2) * go)
        )
    Ipo_m_I0_terms =
        go / (a2 * (1 - go * x4_p)) *
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * S1p_o)
    Imo_m_I0_terms =
        go / (a2 * (1 - go * x4_m)) *
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_m))) * S1m_o)

    return I1o_m_I0_terms, I2o_m_I0_terms, Ipo_m_I0_terms, Imo_m_I0_terms
end

@inline function radial_w_I0_terms_integrals(
    met::Kerr{T},
    rs,
    roots::NTuple{4},
    τ,
    νr,
) where {T}
    numreals = sum(_isreal2.(roots))
    if numreals == 4
        I1, I2, Ip, Im =
            Krang.radial_w_I0_terms_integrals_case2(met, rs, real.(roots), τ, νr)
    elseif numreals == 2
        I1, I2, Ip, Im = Krang.radial_w_I0_terms_integrals_case3(met, rs, roots, τ)
    else
        I1, I2, Ip, Im = Krang.radial_w_I0_terms_integrals_case4(met, rs, roots, τ)
    end
    return I1, I2, Ip, Im
end

"""
Returns the radial integrals for the case where there are four real roots in the radial potential, with roots outside the horizon.
"""
@inline function radial_w_I0_terms_integrals_case2(
    metric::Kerr{T},
    rs,
    roots::NTuple{4},
    τ,
    νr,
) where {T}
    r1, r2, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r43 = r4 - r3
    a = metric.spin
    rp = 1 + √(1 - a^2)
    rm = 1 - √(1 - a^2)
    rp3 = rp - r3
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4

    #FIXME: This is a hack to avoid the division by zero
    if (rp3 == zero(T))
        rp3 = eps(T)
    end

    k = r32 * r41 / (r31 * r42)
    x2_s = √abs((rs - r4) / (rs - r3) * r31 / r41)
    !(x2_s < 1) && return zero(T), zero(T), zero(T), zero(T), zero(T)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    E_s = √(r31 * r42) * JacobiElliptic.E(asin(x2_s), k)
    Π1_s = coef * JacobiElliptic.Pi(n, asin(x2_s), k)

    I0_total = τ

    poly_coefs = (
        r1 * r2 * r3 * r4,
        -r2 * r3 * r4 + r1 * (r2 * (-r3 - r4) - r3 * r4),
        r3 * r4 + r2 * (r3 + r4) + r1 * (r2 + r3 + r4),
        -r1 - r2 - r3 - r4,
        1,
    )

    #equation B37
    I1_total = -r3 * I0_total - r43 * (-1)^νr * Π1_s# Removed the logarithmic divergence
    #equation B38
    I2_s = √abs(evalpoly(rs, poly_coefs)) / (rs - r3) - E_s
    I2_total = (r1 * r4 + r2 * r3) / 2 * τ - (-1)^νr * I2_s# asymptotic Divergent piece is not included

    Πp_s = _Pi_horizon_combined(r31, r42, r43, r41, rp3, rp4, asin(x2_s), k)
    Πm_s = _Pi_horizon_combined(r31, r42, r43, r41, rm3, rm4, asin(x2_s), k)

    Ip_total = τ / rp3 + (-1)^νr * Πp_s
    Im_total = τ / rm3 + (-1)^νr * Πm_s

    return I1_total, I2_total, Ip_total, Im_total

end

"""
Returns the radial integrals for the case where there are two real roots in the radial potential
"""
@inline function radial_w_I0_terms_integrals_case3(
    metric::Kerr{T},
    rs,
    roots::NTuple{4},
    τ,
) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    r1, r2, r21 = real.((r1, r2, r21))
    a = metric.spin
    a2 = a * a
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)
    rp1 = rp - r1
    rp2 = rp - r2
    rm1 = rm - r1
    rm2 = rm - r2

    A2 = abs2(r32)
    B2 = abs2(r31)
    A, B = √A2, √B2
    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    temprat = B * (rs - r2) * inv(A * (rs - r1))
    x3_s = real((1 - temprat) * inv(1 + temprat))

    abs(x3_s) > 1 && return zero(T), zero(T), zero(T), zero(T)

    φ_s = acos(x3_s)
    αo = (B + A) / (B - A)
    coef = 2 * r21 * √(A * B) / (B2 - A2)

    Π1_s = coef * R1(αo, φ_s, k3)
    Π2_s = (coef^2) * R2(αo, φ_s, k3)

    C1p_s = _R1_alpha_combined(A, B, rp1, rp2, r21, φ_s, k3)
    C1m_s = _R1_alpha_combined(A, B, rm1, rm2, r21, φ_s, k3)

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = -(B * r2 + A * r1) / (B + A) * I0_total + Π1_s
    # Removed linear divergence
    I2_total = -(
        (((B * r2 + A * r1) / (B + A))^2) * I0_total -
        2 * (B * r2 + A * r1) / (B + A) * (-Π1_s) - √(A * B) * (-Π2_s)
    )
    Ip_total = inv(B * rp2 + A * rp1) * ((B + A) * I0_total - C1p_s)
    Im_total = inv(B * rm2 + A * rm1) * ((B + A) * I0_total - C1m_s)

    return I1_total, I2_total, Ip_total, Im_total
end

"""
Returns the radial integrals for the case where there are no real roots in the radial potential
"""
@inline function radial_w_I0_terms_integrals_case4(
    metric::Kerr{T},
    rs,
    roots::NTuple{4},
    τ,
) where {T}
    _, r2, _, r4 = roots
    a = metric.spin
    a2 = a * a
    rp = 1 + √(1 - a2)
    rm = 1 - √(1 - a2)

    a1 = abs(imag(r4))
    a2 = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1-a2)^2 + (b1-b2)^2)
    D = sqrt((a1+a2)^2 + (b1-b2)^2)
    k4 = 4C * D / (C + D)^2

    x4_s = (rs + b1) / a2
    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - 1) / (go + x4_p)
    gm = (go * x4_m - 1) / (go + x4_m)

    (isnan(x4_s) || isnan(x4_p) || isnan(x4_m)) && return zero(T), zero(T), zero(T), zero(T)
    # One amplitude for every elliptic integral here; `F` at it is shared (see
    # [`regularized_Pi`](@ref)).
    φ4_s = atan(x4_s) + atan(go)
    F4_s, E4_s = JacobiElliptic.CarlsonAlg._F_E(φ4_s, k4)
    S1p_s = S1(gp, φ4_s, k4, F4_s)
    S1m_s = S1(gm, φ4_s, k4, F4_s)
    #Fr_s = 2 / (C + D) * JacobiElliptic.F(atan(x4_s) + atan(go), k4)
    Π1_s = 2 / (C + D) * (a2 / go * (1 + go^2)) * S1(go, φ4_s, k4, F4_s)
    Π2_s = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * S2(go, φ4_s, k4, F4_s, E4_s)

    # Removed logarithmic divergence
    I1_total = -(a2 / go - b1) * τ + (-Π1_s)

    # Removed linear divergence
    I2_total = -((a2 / go - b1)^2) * τ + 2(a2 / go - b1) * (-Π1_s) - (-Π2_s)
    Ip_total =
        -go / (a2 * (1 - go * x4_p)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (-S1p_s))
    Im_total =
        -go / (a2 * (1 - go * x4_m)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_m))) * (-S1m_s))

    return I1_total, I2_total, Ip_total, Im_total
end

"""
Returns the maximum mino time that can be accrued along a ray.

# Arguments

- `metric` : Kerr metric
- `roots` : Roots of the radial potential
"""
function total_mino_time(metric::Kerr{T}, roots::NTuple{4}) where {T}
    numreals = unsafe_trunc(Int, sum(Krang._isreal2, roots))
    I0_inf = Ir_inf(metric, roots)

    if numreals == 4
        τf = 2I0_inf
    else
        rh = Krang.horizon(metric)
        τf = I0_inf - Krang.Ir_s(metric, rh, roots, true)
    end
    return τf
end

##----------------------------------------------------------------------------------------------------------------------
# Pixel utility functions
##----------------------------------------------------------------------------------------------------------------------

"""
Returns the antiderivative \$I_r=\\int\\frac{dr}{\\sqrt{\\mathcal{R(r)}}}\$.
See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `pix`  : Pixel information
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
- `rs` : Emission radius
"""
function Ir(pix::AbstractPixel, νr::Bool, rs)
    return I0_inf(pix) - Ir_s(metric(pix), rs, roots(pix), νr)
end

"""
Returns the antiderivative \$I_ϕ=\\int\\frac{a(2Mr-a\\lambda)}{\\Delta\\sqrt{\\mathcal{R}(r)}}dr\$.
See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `pix`: SlowLightIntensityPixel
- `rs` : Emission radius
- `τ` : Mino time
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
"""
function Iϕ(pix::AbstractPixel, rs, τ, νr)
    metric = pix.metric
    λtemp = λ(metric, pix.screen_coordinate[1], pix.θo)
    tempIϕ_inf = Iϕ_inf(pix)
    tempIϕ_s = Iϕ_w_I0_terms(metric, rs, τ, pix.roots, νr, λtemp)

    return tempIϕ_inf - tempIϕ_s
end

"""
Returns the antiderivative \$I_t=\\int\\frac{r^2\\Delta+2Mr(r^2+a^2-a\\lambda)}{\\Delta\\sqrt{\\mathcal{R}(r)}}dr\$.
See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `pix`: SlowLightIntensityPixel
- `rs` : Emission radius
- `τ` : Mino time
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
"""
function It(pix::AbstractPixel, rs, τ, νr)
    metric = pix.metric
    λtemp = λ(metric, pix.screen_coordinate[1], pix.θo)
    tempIt_inf = It_inf(pix)

    return tempIt_inf - It_w_I0_terms(metric, rs, τ, pix.roots, λtemp, νr)
end

"""
Return the radial integrals

- `pix`: SlowLightIntensityPixel
- `rs` : Emission radius
- `τ` : Mino time
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
"""
@inline function radial_integrals(pix::AbstractPixel, rs, τ, νr)
    met = metric(pix)
    I1_o, I2_o, Ip_o, Im_o = radial_inf_integrals_m_I0_terms(pix)
    I1_s, I2_s, Ip_s, Im_s = radial_w_I0_terms_integrals(met, rs, pix.roots, τ, νr)
    return τ, I1_o - I1_s, I2_o - I2_s, Ip_o - Ip_s, Im_o - Im_s
end

@inline function _rs_case1_and_2(pix::AbstractPixel, rh, τ::T)::Tuple{T,Bool,Bool} where {T}
    radial_roots = real.(roots(pix))
    _, _, r3, r4 = radial_roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(radial_roots...)

    k = r32 * r41 / (r31 * r42)
    if (rh == r3)
        rh += eps(T)
    end
    err_return = (zero(T), false, false)
    x2_s = √abs((rh - r4) / (rh - r3) * r31 / r41)
    coef = 2 / √real(r31 * r42)
    Ir_s = !(x2_s < 1) ? zero(T) : coef * JacobiElliptic.F(asin(x2_s), k)

    fo = I0_inf(pix)

    r4 < rh && τ > (fo - Ir_s) && return err_return# invalid case2

    X2 = √(r31 * r42) * (fo - τ) / 2
    if τ > 2fo
        return err_return
    end
    sn = r41 * JacobiElliptic.sn(X2, k)^2
    return (r31 * r4 - r3 * sn) / (r31 - sn), X2 > zero(T), true
end

@inline function _rs_case3(pix::AbstractPixel, rh, τ::T)::Tuple{T,Bool,Bool} where {T}
    radial_roots = roots(pix)
    r1, r2, _, _ = radial_roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(radial_roots...)
    r1, r2, r21 = real.((r1, r2, r21))

    err_return = (zero(T), false, false)
    fo = I0_inf(pix)
    A = √abs2(r32)
    B = √abs2(r31)
    k = (((A + B)^2 - r21^2) / (4 * A * B))
    temprat = B * (rh - r2) / (A * (rh - r1))
    x3_s = clamp(((1 - temprat) / (1 + temprat)), -1, 1)
    coef = 1 * √inv(A * B)
    Ir_s = coef * JacobiElliptic.F((acos(x3_s)), k)
    τ > (fo - Ir_s) && return err_return

    X3 = √(A * B) * real(fo - τ)
    if X3 < zero(T)
        return err_return
    end
    cn = JacobiElliptic.cn(X3, k)
    num = -A * r1 + B * r2 + (A * r1 + B * r2) * cn
    den = -A + B + (A + B) * cn

    return real(num / den), X3 > zero(T), true
end

@inline function _rs_case4(pix::AbstractPixel, rh, τ::T)::Tuple{T,Bool,Bool} where {T}
    radial_roots = roots(pix)
    _, r2, _, r4 = radial_roots

    a1 = abs(imag(r4))
    a2 = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1-a2)^2 + (b1-b2)^2)
    D = sqrt((a1+a2)^2 + (b1-b2)^2)
    k4 = 4C * D / (C + D)^2

    go = √max((T(4)a2^2 - (C - D)^2) / ((C + D)^2 - T(4)a2^2), zero(T))
    x4_s = (rh + b1) / a2
    coef = 2 / (C + D)
    Ir_s = coef * JacobiElliptic.F(atan(x4_s) + atan(go), k4)

    fo = I0_inf(pix)
    τ > (fo - Ir_s) && return (zero(T), false, false)

    X4 = (C + D) / T(2) * (fo - τ)
    num = go - JacobiElliptic.sc(X4, k4)
    den = 1 + go * JacobiElliptic.sc(X4, k4)

    return -(a2 * num / den + b1), X4 > zero(T), true
end

##----------------------------------------------------------------------------------------------------------------------
# Inclination functions
##----------------------------------------------------------------------------------------------------------------------

"""
Mino time of trajectory between two inclinations for a given screen coordinate

# Arguments 

- `pix` : Pixel information
- `θs` : Emission inclination
- `isindir` : Is the path direct or indirect?
- `n` : nth image in orde of amount of minotime traversed
"""
function mino_time(pix::AbstractPixel, θs, isindir, n)
    return Gθ(pix, θs, isindir, n)[1]
end

"""
Returns the antiderivative \$G_\\theta=\\int\\frac{d\\theta}{\\sqrt{\\Theta(\\theta)}}\$.
See [`θ_potential(x)`](@ref) for an implementation of \$\\Theta(\theta)\$.

# Arguments 

- `pix` : Pixel information
- `θs` : Emission inclination
- `isindir` : Is the path direct or indirect?
- `n` : nth image ordered by minotime
"""
@inline function Gθ(
    pix::AbstractPixel,
    θs::T,
    isindir,
    n,
)::Tuple{T,T,T,T,Bool,Bool} where {T}
    _, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)
    signβ = sign(β)
    ηtemp = η(pix)
    λtemp = λ(pix)

    a = met.spin
    a2 = a^2
    Go, Ghat = absGθo_Gθhat(pix)
    Gs, minotime, isvortical = zero(T), zero(T), ηtemp < zero(T)

    cosθs = cos(θs)
    cosθo = cos(θo)
    isincone = abs(cosθs) < abs(cosθo)
    if isincone && (isindir != ((signβ > 0) ⊻ (θo > T(π / 2))))
        return minotime, Gs, Go, Ghat, isvortical, false
    end
    if ((((signβ < 0) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical) ||
       (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return minotime, Gs, Go, Ghat, isvortical, false
    end

    Δθ = (1 - (ηtemp + λtemp^2) / a2) / 2
    Δθ2 = Δθ^2
    desc = √(Δθ2 + ηtemp / a2)
    up = min(Δθ + desc, 1 - eps(T))
    um = Δθ - desc
    m = up / um

    if isvortical
        args = (cosθs^2 - um) / (up - um)
        argo = (cosθo^2 - um) / (up - um)
        k = 1 - m
        if (!(zero(T) < argo < 1) || !(zero(T) < args < 1))
            return minotime, Gs, Go, Ghat, isvortical, false
        end
        tempfac = 1 / √abs(um * a2)
        Go *= (θs > T(π / 2) ? -1 : 1)
        Gs = (θs > T(π / 2) ? -1 : 1) * tempfac * JacobiElliptic.F(asin(√args), k)
    else
        args = cosθs / √(up)
        argo = cosθo / √(up)
        k = m
        if !(-1 < args < 1) || !(-1 < argo < 1)
            return zero(T), zero(T), zero(T), zero(T), isvortical, false
        end
        tempfac = 1 / √abs(um * a^2)
        Gs = tempfac * JacobiElliptic.F(asin(args), k)
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    minotime = (
        isindir ? (n + 1) * Ghat - signβ * Go - (-1)^νθ * Gs :
        n * Ghat - signβ * Go - (-1)^νθ * Gs
    ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return minotime, Gs, Go, Ghat, isvortical, true
end

function Gs(pix::AbstractPixel, τ::T) where {T}
    α, β = screen_coordinate(pix)
    θo = inclination(pix)
    met = metric(pix)
    ηtemp = η(met, α, β, θo)
    λtemp = λ(met, α, θo)
    signβ = sign(β)

    τ == T(NaN) && return T(NaN)
    a = met.spin

    Go, Ghat = absGθo_Gθhat(pix)
    Gs, isvortical = zero(T), ηtemp < zero(T)

    Δθ = T(0.5) * (1 - (ηtemp + λtemp^2) / a^2)
    up = min(Δθ + √(Δθ^2 + ηtemp / a^2), 1 - eps(T))
    um = Δθ - √(Δθ^2 + ηtemp / a^2)
    m = up / um
    k = m

    argo = zero(T)
    k = zero(T)
    tempfac = 1 / √abs(um * a^2)

    if isvortical
        argo = (cos(θo)^2 - um) / (up - um)
        k = 1 - m
        tempfac = 1 / √abs(um * a^2)
        Δτtemp = (τ % Ghat + (θo > T(π / 2) ? -1 : 1) * signβ * Go)
        n = floor(τ / Ghat)
        Δτ =
            (θo > T(π / 2) ? -1 : 1) *
            abs(argmin(abs, [(-1)^n * signβ * (Ghat - Δτtemp), (-1)^n * signβ * Δτtemp]))
    else
        argo = cos(θo) / √(up)
        k = m
        tempfac = inv(√abs(um * a^2))
        Δτtemp = (τ % Ghat + signβ * Go)
        n = floor(τ / Ghat)
        Δτ = argmin(abs, [(-1)^n * signβ * (Ghat - Δτtemp), (-1)^n * signβ * Δτtemp])
    end

    return Δτ
end

@inline function Gϕ(pix::AbstractPixel, θs::T, isindir, n) where {T}
    _, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)

    signβ = sign(β)
    ηtemp = η(pix)
    λtemp = λ(pix)

    a = met.spin
    a2 = a^2
    Go, Ghat = absGϕo_Gϕhat(pix)
    Gs, ans, isvortical = zero(T), zero(T), ηtemp < zero(T)

    cosθs = cos(θs)
    cosθo = cos(θo)
    isincone = abs(cosθs) < abs(cosθo)
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return ans, Gs, Go, Ghat, isvortical, false
    end
    if (
        (((signβ < zero(T)) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical
    ) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return ans, Gs, Go, Ghat, isvortical, false
    end

    Δθ = (1 - (ηtemp + λtemp^2) / a2) / T(2)
    desc = √(Δθ^2 + ηtemp / a2)
    up = min(Δθ + desc, 1 - eps(T))
    um = Δθ - desc
    m = up / um
    k = m

    if isvortical
        args = (cosθs^2 - um) / (up - um)
        argo = (cosθo^2 - um) / (up - um)
        k = 1 - m
        if (!(zero(T) < argo < 1) || !(zero(T) < args < 1))
            return ans, Gs, Go, Ghat, isvortical, false
        end
        tempfac = inv((1 - um) * √abs(um * a2))
        argn = (up - um) / (1 - um)
        Go = ((θs > T(π / 2)) ? -1 : 1) * Go
        Gs = ((θs > T(π / 2)) ? -1 : 1) * tempfac * JacobiElliptic.Pi(argn, asin(√args), k)
    else
        args = cosθs / √(up)
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-1 < args < 1) || !(-1 < argo < 1)
            return ans, Gs, Go, Ghat, isvortical, false
        end
        tempfac = inv(√abs(um * a2))
        Gs = tempfac * JacobiElliptic.Pi(up, asin(args), k)
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    ans = (
        isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs :
        n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs
    ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return ans, Gs, Go, Ghat, isvortical, true
end

@inline function Gt(pix::AbstractPixel, θs::T, isindir, n) where {T}
    _, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)
    signβ = sign(β)

    # The pixel already caches these; recomputing them from (α, β, θo) gives the same bits.
    ηtemp = η(pix)
    λtemp = λ(pix)
    a = met.spin
    a2 = a^2
    Go, Ghat = absGto_Gthat(pix)
    Gs, ans, isvortical = zero(T), zero(T), ηtemp < zero(T)

    cosθs = cos(θs)
    cosθo = cos(θo)
    isincone = abs(cosθs) < abs(cosθo)
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return ans, Gs, Go, Ghat, isvortical, false
    end
    if (
        (((signβ < zero(T)) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical
    ) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return ans, Gs, Go, Ghat, isvortical, false
    end

    Δθ = (1 - (ηtemp + λtemp^2) / a2) / T(2)
    desc = √(Δθ^2 + ηtemp / a2)
    up = min(Δθ + desc, 1 - eps(T))
    um = Δθ - desc
    m = up / um
    k = m

    if isvortical
        args = (cosθs^2 - um) / (up - um)
        argo = (cosθo^2 - um) / (up - um)
        k = 1 - m
        if (!(zero(T) < argo < 1) || !(zero(T) < args < 1))
            return ans, Gs, Go, Ghat, isvortical, false
        end
        tempfac = √abs(um / a2)
        Go *= ((θs > T(π / 2)) ? -1 : 1)
        Gs = ((θs > T(π / 2)) ? -1 : 1) * tempfac * JacobiElliptic.E(asin(√args), k)

    else
        args = cosθs / √(up)
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-1 < args < 1) || !(-1 < argo < 1)
            return ans, Gs, Go, Ghat, isvortical, false
        end
        tempfac = -2 * up * inv(√abs(um * a2))
        # `E` and `F` at one amplitude come out of a single Carlson pass, and the Enzyme
        # rule for `_F_E` likewise builds the whole 2x2 Jacobian from that one pair.
        φs = asin(args)
        Fs, Es = JacobiElliptic.CarlsonAlg._F_E(φs, k)
        Gs = tempfac * (Es - Fs) / (2k)
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    ans = (
        isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs :
        n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs
    ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return ans, Gs, Go, Ghat, isvortical, true
end

"""
Returns `(Gϕ, Gt, issuccess)`: the ϕ and t angular antiderivatives of [`Gϕ`](@ref) and
[`Gt`](@ref) for the same emission inclination.

Both share the whole angular setup — `η`, `λ`, the turning points `up`/`um`, the modulus `k`
and the amplitude `asin(args)` — and differ only in which elliptic integral they evaluate at
that amplitude. The `emission_coordinates` methods always want both, so computing the setup
once avoids doing it twice. Results are bit-identical to calling `Gϕ` and `Gt` separately.
"""
@inline function _Gϕt(pix::AbstractPixel, θs::T, isindir, n)::Tuple{T,T,Bool} where {T}
    _, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)
    signβ = sign(β)
    ηtemp = η(pix)
    λtemp = λ(pix)

    a = met.spin
    a2 = a^2
    isvortical = ηtemp < zero(T)
    z = zero(T)

    cosθs = cos(θs)
    cosθo = cos(θo)
    isincone = abs(cosθs) < abs(cosθo)
    if isincone && (isindir != ((signβ > z) ⊻ (θo > T(π / 2))))
        return z, z, false
    end
    if ((((signβ < z) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical) ||
       (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return z, z, false
    end

    Δθ = (1 - (ηtemp + λtemp^2) / a2) / T(2)
    desc = √(Δθ^2 + ηtemp / a2)
    up = min(Δθ + desc, 1 - eps(T))
    um = Δθ - desc
    m = up / um

    if isvortical
        args = (cosθs^2 - um) / (up - um)
        argo = (cosθo^2 - um) / (up - um)
        k = 1 - m
        if (!(z < argo < 1) || !(z < args < 1))
            return z, z, false
        end
        φs = asin(√args)
        sgn = (θs > T(π / 2)) ? -1 : 1
        Gϕs =
            sgn *
            inv((1 - um) * √abs(um * a2)) *
            JacobiElliptic.Pi((up - um) / (1 - um), φs, k)
        Gts = sgn * √abs(um / a2) * JacobiElliptic.E(φs, k)
        # Fetched only here: on a pixel that does not cache them these are a full
        # recompute, and every early return above would otherwise pay for it.
        Gϕo_t, Gϕhat = absGϕo_Gϕhat(pix)
        Gto_t, Gthat = absGto_Gthat(pix)
        Gϕo = sgn * Gϕo_t
        Gto = sgn * Gto_t
    else
        args = cosθs / √(up)
        argo = cosθo / √(up)
        k = m
        if !(-1 < args < 1) || !(-1 < argo < 1)
            return z, z, false
        end
        φs = asin(args)
        tempfac = inv(√abs(um * a2))
        Gϕs = tempfac * JacobiElliptic.Pi(up, φs, k)
        Fs, Es = JacobiElliptic.CarlsonAlg._F_E(φs, k)
        Gts = (-2 * up * tempfac) * (Es - Fs) / (2k)
        Gϕo, Gϕhat = absGϕo_Gϕhat(pix)
        Gto, Gthat = absGto_Gthat(pix)
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    sgnθ = νθ ? 1 : -1
    norb = isindir ? n + 1 : n
    return norb * Gϕhat - signβ * Gϕo + sgnθ * Gϕs,
    norb * Gthat - signβ * Gto + sgnθ * Gts,
    true
end

##----------------------------------------------------------------------------------------------------------------------
# SlowLightIntensityCachedPixel utility functions
##----------------------------------------------------------------------------------------------------------------------

@inline function _absGθo_Gθhat(metric::Kerr{T}, θo, η, λ)::NTuple{2,T} where {T}
    a = metric.spin
    a2 = a^2
    Go, Ghat, isvortical = zero(T), zero(T), η < zero(T)

    Δθ = (1 - (η + λ^2) / a2) / 2
    Δθ2 = Δθ^2
    desc = √max(Δθ2 + η / a2, zero(T))
    up = min(Δθ + desc, 1 - eps(T))
    um = Δθ - desc
    m = up / um

    cosθo = cos(θo)
    if isvortical
        argo = (cosθo^2 - um) / (up - um)
        k = 1 - m
        if (!(zero(T) < argo < 1))
            return Go, Ghat
        end
        tempfac = 1 / √abs(um * a2)
        Go = tempfac * JacobiElliptic.F(asin(√argo), k)
        Ghat = 2 * tempfac * JacobiElliptic.K(k)
    else
        argo = cosθo / √(up)
        k = m
        if !(-1 < argo < 1)
            return Go, Ghat
        end
        tempfac = 1 / √abs(um * a2)
        Go = tempfac * JacobiElliptic.F(asin(argo), k)
        Ghat = 2 * tempfac * JacobiElliptic.K(k)
    end

    return Go, Ghat
end

@inline function _absGϕo_Gϕhat(metric::Kerr{T}, θo, η, λ)::NTuple{2,T} where {T}

    a = metric.spin
    a2 = a^2
    Go, Ghat, isvortical = zero(T), zero(T), η < zero(T)

    Δθ = (1 - (η + λ^2) / a2) / T(2)
    desc = √(Δθ^2 + η / a2)
    up = min(Δθ + desc, 1 - eps(T))
    um = Δθ - desc
    m = up / um
    k = m

    cosθo = cos(θo)
    if isvortical
        argo = (cosθo^2 - um) / (up - um)
        k = 1 - m
        if (!(zero(T) < argo < 1))
            return Go, Ghat
        end
        tempfac = inv((1 - um) * √abs(um * a2))
        argn = (up - um) / (1 - um)
        Go = tempfac * JacobiElliptic.Pi(argn, asin(√argo), k)
        Ghat = 2tempfac * JacobiElliptic.Pi(argn, k)
    else
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-1 < argo < 1)
            return Go, Ghat
        end
        tempfac = inv(√abs(um * a2))

        Go = tempfac * JacobiElliptic.Pi(up, asin(argo), k)
        Ghat = 2tempfac * JacobiElliptic.Pi(up, k)
    end

    return Go, Ghat
end

@inline function _absGto_Gthat(metric::Kerr{T}, θo, η, λ)::NTuple{2,T} where {T}
    a = metric.spin
    a2 = a^2
    Go, Ghat, isvortical = zero(T), zero(T), η < zero(T)

    Δθ = (1 - (η + λ^2) / a2) / T(2)
    desc = √(Δθ^2 + η / a2)
    up = min(Δθ + desc, 1 - eps(T))
    um = Δθ - desc
    m = up / um
    k = m

    cosθo = cos(θo)
    if isvortical
        argo = (cosθo^2 - um) / (up - um)
        k = 1 - m
        if (!(zero(T) < argo < 1))
            return Go, Ghat
        end
        tempfac = √abs(um / a2)
        Go = tempfac * JacobiElliptic.E(asin(√argo), k)
        Ghat = 2tempfac * JacobiElliptic.E(k)

    else
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-1 < argo < 1)
            return Go, Ghat
        end
        tempfac = -2 * up * inv(√abs(um * a2))
        φo = asin(argo)
        Fo, Eo = JacobiElliptic.CarlsonAlg._F_E(φo, k)
        Go = tempfac * (Eo - Fo) / (2k)
        # `K(k)` and `E(k)` each run the whole AGM pass and discard the other half.
        # `ellipke` returns both from one pass and, since JacobiElliptic gained Enzyme
        # rules for it, costs the same under `Reverse` as the two separate calls.
        Kc, Ec = JacobiElliptic.ArithmeticGeometricMeanAlg.ellipke(k)
        Ghat = tempfac * (Ec - Kc) / k
    end

    return Go, Ghat
end
