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
Checks if a complex number is real to some tolerance
"""
function _isreal2(num::Complex{T}) where {T}
    ren, imn = reim(num)
    ren2 = ren^2
    imn2 = imn^2
    return imn2 / (imn2 + ren2) < eps(T)
end

"""
Regularized elliptic integral of the third kind
    
# Arguments

- `n`: Parameter
- `ϕ`: Arguments
- `k`: Parameter
"""
function regularized_Pi(n, ϕ, k)
    ω = k / n
    return JacobiElliptic.F(ϕ, k) - JacobiElliptic.Pi(ω, ϕ, k)
end

function regularized_R1(α::T, φ::T, j) where {T}
    n = α^2 / (α^2 - one(T))
    return one(T) / (one(T) - α^2) * (regularized_Pi(n, φ, j))
end

function regularized_R2(α::T, φ::T, j) where {T}

    return one(T) / (one(T) - α^2) * (
        JacobiElliptic.F(φ, j) - α^2 / (j + (one(T) - j) * α^2) * (JacobiElliptic.E(φ, j)
        #-α*sin(φ)*sqrt(one(T)-j*sin(φ)^2)/(one(T)+α*cos(φ)) #Linear Divergent term
        )
    ) #+ 
    #inv(j+(one(T)-j)*α^2)*(2*j-α^2/(α^2-one(T)))*regularized_R1(α, φ, j)
end

function regularizedS1(α, φ, j)
    α2 = α * α
    return inv(1 + α2) * (JacobiElliptic.F(φ, j) + α2 * regularized_Pi(1 + α2, φ, j))# - α*f2(α, sin(φ), j)) # logarithmic divergence removed
end

function regularizedS2(α, φ, j)

    return -inv((1 + α^2) * (1 - j + α^2)) *
           ((1 - j) * JacobiElliptic.F(φ, j) + α^2 * JacobiElliptic.E(φ, j)) +# + α^2*√(1-j*sin(φ)^2)-α^3) +
           (inv(1 + α^2) + (1 - j) / (1 - j + α^2)) * regularizedS1(α, φ, j)
end

function p1(α::T, j::T) where {T}
    √abs((α^2 - one(T)) / (j + (one(T) - j) * α^2))
end

function f1(α::T, sinφ::T, j) where {T}
    p1temp = p1(α, j)
    tempsinφ = √(one(T) - j * sinφ^2)
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

function R1(α::T, φ::T, j) where {T}
    #FIXME: This function is undefined when n=1 in Pi(n, ϕ, m) and when α^2 =1
    denom = ((one(T) - α^2) + eps(T))
    return one(T) / denom *
           (JacobiElliptic.Pi(-α^2 / denom + eps(T), φ, j) - α * f1(α, sin(φ), j))
end

function R2(α::T, φ::T, j) where {T}
    #FIXME: This function is undefined when α*cos(φ) = 1
    denom = ((one(T) - α^2) + eps(T))
    return one(T) / denom * (
        JacobiElliptic.F(φ, j) -
        α^2 / (j + (one(T) - j) * α^2) * (
            JacobiElliptic.E(φ, j) -
            α * sin(φ) * √(one(T) - j * sin(φ)^2) / ((one(T) + α * cos(φ)) + eps(T))
        )
    ) +
           inv(j + (one(T) - j) * α^2) *
           (2 * j + α^2 / ((one(T) - α^2) + eps(T))) *
           R1(α, φ, j)
end

function S1(α, φ, j)
    α2 = α * α
    return inv(1 + α2) * (
        JacobiElliptic.F(φ, j) + α2 * JacobiElliptic.Pi(1 + α2, φ, j) -
        α * f2(α, sin(φ), j)
    )
end

function S2(α, φ, j)
    return -inv((1 + α^2) * (1 - j + α^2)) * (
        (1 - j) * JacobiElliptic.F(φ, j) +
        α^2 * JacobiElliptic.E(φ, j) +
        α^2 * √(1 - j * sin(φ)^2) * (α - tan(φ)) / (1 + α * tan(φ)) - α^3
    ) + (inv(1 + α^2) + (1 - j) / (1 - j + α^2)) * S1(α, φ, j)
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
- `λ`: Energy reduced Azimuthal angular momentul
- `η`: Energy reduced Carter integral 
- `θo`: Observer inclination
"""
function β(metric::Kerr, λ, η, θo)
    return sqrt(max(η - (α(metric, λ, θo)^2 - metric.spin^2) * cos(θo)^2, 0))
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
@inline function βboundary(metric::Kerr{T}, α, θo, θs) where {T}
    a = metric.spin
    cosθs2 = cos(θs)^2
    return √max(
        (cos(θo)^2 - cosθs2) * (α^2 - a^2 * (one(T) - cosθs2)) / (cosθs2 - one(T)),
        zero(T),
    ) #eq 15 DOI 10.3847/1538-4357/acafe3 
end

"""
Radial potential of spacetime

# Arguments

- `metric`: Kerr{T} metric
- `η`  : Reduced Carter constant
- `λ`  : Reduced azimuthal angular momentum
- `r`  : Boyer Lindquist radius
"""
function r_potential(metric::Kerr{T}, η, λ, r) where {T}
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
function θ_potential(metric::Kerr{T}, η, λ, θ) where {T}
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
@inline function get_radial_roots(metric::Kerr{T}, η, λ) where {T}
    a = metric.spin

    a2 = a * a
    A = a2 - η - λ * λ
    A2 = A + A
    B = T(2) * (η + (λ - a)^2)
    C = -a2 * η

    P = -A * A / T(12) - C
    Q = -A / T(3) * (A * A / T(36) + zero(T)im - C) - B * B / T(8)

    Δ3 = -T(4) * P * P * P - T(27) * Q * Q
    ωp = ^(-Q / T(2) + sqrt(-Δ3 / T(108)) + zero(T)im, T(1 / 3))

    #C = ((-1+0im)^(2/3), (-1+0im)^(4/3), 1) .* ωp
    C = (-T(1 / 2) + T(√3 / 2)im, -T(1 / 2) - T(√3 / 2)im, one(T) + zero(T)im) .* ωp

    v = -P .* inv.(T(3) .* C)

    ξ0 = argmax(real, (C .+ v)) - A / T(3)
    ξ02 = ξ0 + ξ0

    predet1 = A2 + ξ02
    predet2 = (√T(2) * B) * inv(sqrt(ξ0))
    det1 = sqrt(-(predet1 - predet2))
    det2 = sqrt(-(predet1 + predet2))

    sqrtξ02 = sqrt(ξ02)

    r1 = (-sqrtξ02 - det1) / 2
    r2 = (-sqrtξ02 + det1) / 2
    r3 = (sqrtξ02 - det2) / 2
    r4 = (sqrtξ02 + det2) / 2

    roots = (r1, r2, r3, r4)
    if (sum(_isreal2.(roots)) == 2) && (abs(imag(roots[4])) < sqrt(eps(T)))
        roots = (roots[1], roots[4], roots[2], roots[3])
    end
    return roots
end

@inline function _get_root_diffs(r1::T, r2::T, r3::T, r4::T) where {T}
    return r2 - r1, r3 - r1, r3 - r2, r4 - r1, r4 - r2, r4 - r3
end

"""
Mino time of trajectory between an observer at infinity and point at radius rs

# Arguments 

- `pix` : Pixel information
- `rs` : Emission radius
- `isindir` : Is the path direct or indirect?
"""
function mino_time(pix, rs::T, isindir) where {T}
    return Ir(pix, isindir, rs)[1]
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
    return T(Inf)
end

function Ir_inf_case1_and_2(::Kerr{T}, roots::NTuple{4}) where {T}
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    k = r32 * r41 / (r31 * r42)
    coef = 2 / √real(r31 * r42)
    Ir_inf = coef * JacobiElliptic.F(asin(√(r31 / r41)), k)

    return Ir_inf
end

function Ir_inf_case3(::Kerr{T}, roots::NTuple{4}) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    r1, r2, r21 = real.((r1, r2, r21))

    A2 = abs(r32 * r42)
    B2 = abs(r31 * r41)
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    coef = one(T) * √inv(A * B)
    Ir_inf = coef * JacobiElliptic.F((acos(clamp((A - B) / (A + B), -one(T), one(T)))), k3)

    return Ir_inf
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
    coef = 2 / (C + D)

    Ir_inf = coef * JacobiElliptic.F(T(π / 2) + atan(go), k4)
    #return (C+D), coef
    return Ir_inf
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
    return T(Inf)
end

function Ir_s_case1_and_2(::Kerr{T}, rs, roots::NTuple{4}, νr) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    k = r32 * r41 / (r31 * r42)
    x2_s = √abs((rs - r4) / (rs - r3) * r31 / r41)
    coef = 2 / √real(r31 * r42)
    Ir_s = (x2_s > one(T)) ? T(Inf) : coef * JacobiElliptic.F(asin(x2_s), k)

    return -(-1)^νr * Ir_s
end

function Ir_s_case3(::Kerr{T}, rs, roots::NTuple{4}) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    r1, r2, r21 = real.((r1, r2, r21))

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    temprat = B * (rs - r2) / (A * (rs - r1))
    x3_s = clamp(((one(T) - temprat) / (one(T) + temprat)), -one(T), one(T))
    coef = one(T) * √inv(A * B)
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
    end #case4
    return Iϕ_inf_case4(metric, roots, λ)
end

function Iϕ_inf_case2(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, r43 = _get_root_diffs(roots...)
    a = metric.spin
    a2 = a * a
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)
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

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_o = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)

    arg = rm3 * r41 / (rm4 * r31)
    # TODO: Come up with a better solve for this
    if arg == k # Fix for ∂_m Pi blowing up when n = m
        arg += eps(T)
    end
    Πm_o = coef_m * JacobiElliptic.Pi(arg, asin(x2_o), k)

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
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    if A2 < zero(T) || B2 < zero(T)
        return T(Inf)
    end
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)

    x3_o = clamp((A - B) / (A + B), -one(T), one(T))
    φ_o = acos(x3_o)

    (isnan(φ_o)) && return T(NaN)

    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

    R1p_o = R1(αp, φ_o, k3)
    R1m_o = R1(αm, φ_o, k3)

    Ip = -inv(B * rp2 + A * rp1) * (2 * r21 * √(A * B) / (B * rp2 - A * rp1) * R1p_o)
    Im = -inv(B * rm2 + A * rm1) * (2 * r21 * √(A * B) / (B * rm2 - A * rm1) * R1m_o)

    return (2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im))
end

function Iϕ_inf_case4(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    _, r2, _, r4 = roots
    a = metric.spin
    a2 = a * a
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)

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
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    S1p_o = S1(gp, T(π / 2) + atan(go), k4)
    S1m_o = S1(gm, T(π / 2) + atan(go), k4)

    Ip =
        go / (a2 * (1 - go * x4_p)) *
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1p_o))
    Im =
        go / (a2 * (1 - go * x4_m)) *
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1m_o))

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
    end #case4
    return Iϕ_w_I0_terms_case4(metric, rs, τ, roots, λ)
end

function Iϕ_w_I0_terms_case2(metric::Kerr{T}, rs, τ, roots::NTuple{4}, νr, λ) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, r43 = _get_root_diffs(roots...)
    a = metric.spin
    a2 = a * a
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)
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

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_s), k)
    arg = rm3 * r41 / (rm4 * r31)
    # TODO: Come up with a better solve for this
    if arg == k # Fix for ∂_m Pi blowing up when n = m
        arg += eps(T)
    end
    Πm_s = coef_m * JacobiElliptic.Pi(arg, asin(x2_s), k)

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
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    if A2 < zero(T) || B2 < zero(T)
        return T(Inf)
    end
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)

    temprat = B * (rs - r2) * real(inv(A * (rs - r1)))
    x3_s = clamp(((one(T) - temprat) * real(inv(one(T) + temprat))), -one(T), one(T))
    φ_s = acos(x3_s)

    (isnan(φ_s)) && return T(NaN)

    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

    R1p_s = R1(αp, φ_s, k3)
    R1m_s = R1(αm, φ_s, k3)

    Ip =
        -inv(B * rp2 + A * rp1) *
        ((B + A) * τ + 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (-R1p_s))
    Im =
        -inv(B * rm2 + A * rm1) *
        ((B + A) * τ + 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (-R1m_s))

    return -(2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im))
end

function Iϕ_w_I0_terms_case4(metric::Kerr{T}, rs, τ, roots::NTuple{4}, λ) where {T}
    _, r2, _, r4 = roots
    a = metric.spin
    a2 = a * a
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)

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
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_s) || isnan(x4_p) || isnan(x4_m)) && return T(Inf)
    S1p_s = S1(gp, atan(x4_s) + atan(go), k4)
    S1m_s = S1(gm, atan(x4_s) + atan(go), k4)

    Ip =
        go / (a2 * (1 - go * x4_p)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (-S1p_s))
    Im =
        go / (a2 * (1 - go * x4_m)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (-S1m_s))

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
    end #case4

    return It_inf_case4(metric, roots, λ)
end

function It_inf_case2(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r43 = r4 - r3
    a = metric.spin
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)
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
    E_o = √(r31 * r42) * JacobiElliptic.E(asin(x2_o), k)

    #equation B37
    I1_total =
        log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k))# Removed the logarithmic divergence
    #equation B38
    I2_total = r3 - E_o # asymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_o = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)
    arg = rm3 * r41 / (rm4 * r31)
    # TODO: Come up with a better solve for this
    if arg == k # Fix for ∂_m Pi blowing up when n = m
        arg += eps(T)
    end
    Πm_o = coef_m * JacobiElliptic.Pi(arg, asin(x2_o), k) # Seems to cause NaNs

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
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    if A2 < zero(T) || B2 < zero(T)
        return T(Inf)
    end
    A, B = √A2, √B2

    k3 = real(((A + B)^2 - r21^2) / (4 * A * B))

    x3_o = min((A - B) / (A + B), one(T))
    φ_o = acos(x3_o)

    (isnan(φ_o)) && return T(NaN)

    αo = (B + A) / (B - A)
    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

    Π1_o = 2 * r21 * √real(A * B) / (B2 - A2) * regularized_R1(αo, φ_o, k3) # Divergence is removed, will be added back in the end
    Π2_o = ((2 * r21 * √(A * B) / (B2 - A2))^2) * regularized_R2(αo, φ_o, k3)# Divergence is removed, will be added back in the end

    # Removed logarithmic divergence
    I1_total = Π1_o + log(16 * r21^2 / ((A2 - B2)^2 + 4 * A * B * r21^2)) / 2
    # Removed linear divergence
    I2_total = (-√(A * B) * (Π2_o)) + (B * r2 + A * r1) / (A + B)
    Ip_total =
        -inv(B * rp2 + A * rp1) *
        (2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (R1(αp, φ_o, k3)))
    Im_total =
        -inv(B * rm2 + A * rm1) *
        (2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (R1(αm, φ_o, k3)))

    return -(
        4 / (rp - rm) *
        (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
        2 * I1_total +
        I2_total
    )
end

function It_inf_case4(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    a = metric.spin
    _, r2, _, r4 = roots

    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

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
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    S1p_o = S1(gp, T(π / 2) + atan(go), k4)
    S1m_o = S1(gm, T(π / 2) + atan(go), k4)

    Π1_o = 2 / (C + D) * (a2 / go * (1 + go^2)) * regularizedS1(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end
    Π2_o =
        2 / (C + D) * (a2 / go * (1 + go^2))^2 * regularizedS2(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end

    I1_total =
        -(Π1_o) +
        1 / 2 * log(
            (16 * (1 + go^2 - sqrt((1 + go^2) * (1 + go^2 - k4))) * (1 + go^2 - k4)) /
            ((C + D)^2 * ((1 + go^2)^2 - k4) * k4 * (1 + sqrt(1 - k4 / (1 + go^2)))),
        )

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
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1m_o))

    return -(
        4 / (rp - rm) *
        (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
        2 * I1_total +
        I2_total
    )# + (logdiv + lineardiv)
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
    end #case4

    return It_w_I0_terms_case4(metric, rs, τ, roots, λ)
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
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)
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
        one(T),
    )

    #equation B37
    I1_total = r3 * I0_total + r43 * (-1)^νr * Π1_s# Removed the logarithmic divergence
    #equation B38
    I2_s = √(evalpoly(rs, poly_coefs)) / (rs - r3) - E_s
    I2_total = -(r1 * r4 + r2 * r3) / 2 * τ + (-1)^νr * I2_s# asymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_s), k)
    arg = rm3 * r41 / (rm4 * r31)
    # TODO: Come up with a better solve for this
    if arg == k # Fix for ∂_m Pi blowing up when n = m
        arg += eps(T)
    end
    Πm_s = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_s), k)

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
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    if A2 < zero(T) || B2 < zero(T)
        return T(Inf)
    end
    A, B = √A2, √B2

    k3 = real(((A + B)^2 - r21^2) / (4 * A * B))

    temprat = B * (rs - r2) * inv(A * (rs - r1))
    x3_s = clamp(real(((one(T) - temprat) * inv(one(T) + temprat))), -one(T), one(T))
    φ_s = acos(x3_s)

    (isnan(φ_s)) && return T(NaN)

    αo = (B + A) / (B - A)
    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

    Π1_s = 2 * r21 * √real(A * B) / (B2 - A2) * R1(αo, φ_s, k3)
    Π2_s = ((2 * r21 * √(A * B) / (B2 - A2))^2) * R2(αo, φ_s, k3)

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = (B * r2 + A * r1) / (B + A) * I0_total - Π1_s
    # Removed linear divergence
    I2_total = (
        (((B * r2 + A * r1) / (B + A))^2) * I0_total -
        2 * (B * r2 + A * r1) / (B + A) * (-Π1_s) - √(A * B) * (-Π2_s)
    ) #+ (B * r2 + A * r1) / (A + B)
    Ip_total =
        -inv(B * rp2 + A * rp1) *
        ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (-R1(αp, φ_s, k3)))
    Im_total =
        -inv(B * rm2 + A * rm1) *
        ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (-R1(αm, φ_s, k3)))

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
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

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
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_s) || isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    S1p_s = S1(gp, atan(x4_s) + atan(go), k4)
    S1m_s = S1(gm, atan(x4_s) + atan(go), k4)

    Π1_s = 2 / (C + D) * (a2 / go * (1 + go^2)) * S1(go, atan(x4_s) + atan(go), k4)
    Π2_s = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * S2(go, atan(x4_s) + atan(go), k4)

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
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (-S1m_s))

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
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)
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
    E_o = √(r31 * r42) * JacobiElliptic.E(asin(x2_o), k)

    #equation B37
    I1o_m_I0_terms =
        log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k))
    #equation B38
    I2o_m_I0_terms = r3 - E_o# asymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)

    Ipo_m_I0_terms =
        -coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31) + eps(T), asin(x2_o), k)
    arg = rm3 * r41 / (rm4 * r31)
    # TODO: Come up with a better solve for this
    if arg == k # Fix for ∂_m Pi blowing up when n = m
        arg += eps(T)
    end
    Imo_m_I0_terms = -coef_m * JacobiElliptic.Pi(arg, asin(x2_o), k)

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
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)
    rp1 = rp - r1
    rp2 = rp - r2
    rm1 = rm - r1
    rm2 = rm - r2

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    A, B = √A2, √B2
    k3 = ((A + B)^2 - r21^2) / (4 * A * B)

    αo = (B + A) / (B - A)
    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

    φ_o = acos((A - B) / (A + B))
    Π1_o = 2 * r21 * √real(A * B) / (B2 - A2) * regularized_R1(αo, φ_o, k3)
    Π2_o = ((2 * r21 * √(A * B) / (B2 - A2))^2) * regularized_R2(αo, φ_o, k3)

    I1o_m_I0_terms = Π1_o + log(16 * r21^2 / ((A2 - B2)^2 + 4 * A * B * r21^2)) / 2
    # Removed linear divergence
    I2o_m_I0_terms = -√(A * B) * Π2_o + (B * r2 + A * r1) / (A + B)
    Ipo_m_I0_terms =
        -inv(B * rp2 + A * rp1) *
        (2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (R1(αp, φ_o, k3)))
    Imo_m_I0_terms =
        -inv(B * rm2 + A * rm1) *
        (2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (R1(αm, φ_o, k3)))

    return I1o_m_I0_terms, I2o_m_I0_terms, Ipo_m_I0_terms, Imo_m_I0_terms
end

"""
Returns the radial integrals for the case where there are no real roots in the radial potential
"""
@inline function radial_inf_integrals_case4(metric::Kerr{T}, roots::NTuple{4}) where {T}
    _, r2, _, r4 = roots
    a2 = metric.spin^2
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)

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
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    S1p_o = S1(gp, T(π / 2) + atan(go), k4)
    S1m_o = S1(gm, T(π / 2) + atan(go), k4)
    Π1_o = 2 / (C + D) * (a2 / go * (1 + go^2)) * regularizedS1(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end
    Π2_o =
        2 / (C + D) * (a2 / go * (1 + go^2))^2 * regularizedS2(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end

    # Removed logarithmic divergence
    I1o_m_I0_terms =
        -Π1_o +
        1 / 2 * log(
            (16 * (1 + go^2 - sqrt((1 + go^2) * (1 + go^2 - k4))) * (1 + go^2 - k4)) /
            ((C + D)^2 * ((1 + go^2)^2 - k4) * k4 * (1 + sqrt(1 - k4 / (1 + go^2)))),
        )

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
        (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * S1m_o)

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
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)
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
        one(T),
    )

    #equation B37
    I1_total = -r3 * I0_total - r43 * (-1)^νr * Π1_s# Removed the logarithmic divergence
    #equation B38
    I2_s = √abs(evalpoly(rs, poly_coefs)) / (rs - r3) - E_s
    I2_total = (r1 * r4 + r2 * r3) / 2 * τ - (-1)^νr * I2_s# asymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_s), k)
    arg = rm3 * r41 / (rm4 * r31)
    # TODO: Come up with a better solve for this
    if arg == k # Fix for ∂_m Pi blowing up when n = m
        arg += eps(T)
    end
    Πm_s = coef_m * JacobiElliptic.Pi(arg, asin(x2_s), k)

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
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)
    rp1 = rp - r1
    rp2 = rp - r2
    rm1 = rm - r1
    rm2 = rm - r2

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    A, B = √A2, √B2
    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    temprat = B * (rs - r2) * inv(A * (rs - r1))
    x3_s = real((one(T) - temprat) * inv(one(T) + temprat))

    abs(x3_s) > one(T) && return zero(T), zero(T), zero(T), zero(T)

    φ_s = acos(x3_s)
    αo = (B + A) / (B - A)
    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)
    coef = 2 * r21 * √(A * B) / (B2 - A2)

    Π1_s = coef * R1(αo, φ_s, k3)
    Π2_s = (coef^2) * R2(αo, φ_s, k3)

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = -(B * r2 + A * r1) / (B + A) * I0_total + Π1_s
    # Removed linear divergence
    I2_total = -(
        (((B * r2 + A * r1) / (B + A))^2) * I0_total -
        2 * (B * r2 + A * r1) / (B + A) * (-Π1_s) - √(A * B) * (-Π2_s)
    )
    Ip_total =
        inv(B * rp2 + A * rp1) *
        ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (-R1(αp, φ_s, k3)))
    Im_total =
        inv(B * rm2 + A * rm1) *
        ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (-R1(αm, φ_s, k3)))

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
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)

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
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_s) || isnan(x4_p) || isnan(x4_m)) && return zero(T), zero(T), zero(T), zero(T)
    S1p_s = S1(gp, atan(x4_s) + atan(go), k4)
    S1m_s = S1(gm, atan(x4_s) + atan(go), k4)
    #Fr_s = 2 / (C + D) * JacobiElliptic.F(atan(x4_s) + atan(go), k4)
    Π1_s = 2 / (C + D) * (a2 / go * (1 + go^2)) * S1(go, atan(x4_s) + atan(go), k4)
    Π2_s = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * S2(go, atan(x4_s) + atan(go), k4)

    # Removed logarithmic divergence
    I1_total = -(a2 / go - b1) * τ + (-Π1_s)

    # Removed linear divergence
    I2_total = -((a2 / go - b1)^2) * τ + 2(a2 / go - b1) * (-Π1_s) - (-Π2_s)
    Ip_total =
        -go / (a2 * (1 - go * x4_p)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (-S1p_s))
    Im_total =
        -go / (a2 * (1 - go * x4_m)) *
        (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (-S1m_s))

    return I1_total, I2_total, Ip_total, Im_total
end

"""
Returns the maximum mino time that can be accrued along a ray.

# Arguments

- `metric` : Kerr metric
- `roots` : Roots of the radial potential
- `I0_inf` : Mino time at infinity
"""
function total_mino_time(metric::Kerr{T}, roots::NTuple{4}) where {T}
    numreals = unsafe_trunc(Int, sum((Krang._isreal2.(roots))))
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
Returns the antiderivative \$I_ϕ=\\int\\frac{a(2Mr-a\\lambda)}{\\sqrt{\\Delta\\mathcal{R(r)}}}dr\$.
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
Returns the antiderivative \$I_t=\\int\\frac{a(2Mr-a\\lambda)}{\\sqrt{\\Delta\\mathcal{R(r)}}}dr\$.
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
    x2_s = √abs((rh - r4) / (rh - r3) * r31 / r41)
    coef = 2 / √real(r31 * r42)
    Ir_s = !(x2_s < one(T)) ? zero(T) : coef * JacobiElliptic.F(asin(x2_s), k)

    fo = I0_inf(pix)

    r4 < rh && τ > (fo - Ir_s) && return zero(T), false, false# invalid case2

    X2 = √(r31 * r42) * (fo - τ) / 2
    if τ > 2fo
        return zero(T), false, false
    end
    sn = r41 * JacobiElliptic.sn(X2, k)^2
    return (r31 * r4 - r3 * sn) / (r31 - sn), X2 > zero(T), true
end

@inline function _rs_case3(pix::AbstractPixel, rh, τ::T)::Tuple{T,Bool,Bool} where {T}
    radial_roots = roots(pix)
    r1, r2, _, _ = radial_roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(radial_roots...)
    r1, r2, r21 = real.((r1, r2, r21))


    fo = I0_inf(pix)
    A = √abs(r32 * r42)
    B = √abs(r31 * r41)
    k = (((A + B)^2 - r21^2) / (4 * A * B))
    temprat = B * (rh - r2) / (A * (rh - r1))
    x3_s = clamp(((one(T) - temprat) / (one(T) + temprat)), -one(T), one(T))
    coef = one(T) * √inv(A * B)
    Ir_s = coef * JacobiElliptic.F((acos(x3_s)), k)
    τ > (fo - Ir_s) && return zero(T), false, false

    X3 = √(A * B) * real(fo - τ)
    if X3 < zero(T)
        return zero(T), false, false
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


    k4 = T(4) * C * D / (C + D)^2

    go = √max((T(4)a2^2 - (C - D)^2) / ((C + D)^2 - T(4)a2^2), zero(T))
    x4_s = (rh + b1) / a2
    coef = 2 / (C + D)
    Ir_s = coef * JacobiElliptic.F(atan(x4_s) + atan(go), k4)

    fo = I0_inf(pix)
    τ > (fo - Ir_s) && return zero(T), false, false

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

    Δθ = (one(T) - (ηtemp + λtemp^2) / a2) / 2
    Δθ2 = Δθ^2
    desc = √(Δθ2 + ηtemp / a2)
    up = min(Δθ + desc, one(T) - eps(T))
    um = Δθ - desc
    m = up / um
    k = m

    #isvortical = η < 0.
    args = zero(T)
    argo = zero(T)
    k = zero(T)
    if isvortical
        args = (cosθs^2 - um) / (up - um)
        argo = (cosθo^2 - um) / (up - um)
        k = one(T) - m
        if (!(zero(T) < argo < one(T)) || !(zero(T) < args < one(T)))
            return minotime, Gs, Go, Ghat, isvortical, false
        end
        tempfac = one(T) / √abs(um * a2)
        Go *= (θs > T(π / 2) ? -one(T) : one(T))
        Gs = (θs > T(π / 2) ? -one(T) : one(T)) * tempfac * JacobiElliptic.F(asin(√args), k)
    else
        args = cosθs / √(up)
        argo = cosθo / √(up)
        k = m
        if !(-one(T) < args < one(T)) || !(-one(T) < argo < one(T))
            return zero(T), zero(T), zero(T), zero(T), isvortical, false
        end
        tempfac = one(T) / √abs(um * a^2)
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

    Δθ = T(0.5) * (one(T) - (ηtemp + λtemp^2) / a^2)
    up = min(Δθ + √(Δθ^2 + ηtemp / a^2), one(T) - eps(T))
    um = Δθ - √(Δθ^2 + ηtemp / a^2)
    m = up / um
    k = m

    argo = zero(T)
    k = zero(T)
    tempfac = one(T) / √abs(um * a^2)

    if isvortical
        argo = (cos(θo)^2 - um) / (up - um)
        k = one(T) - m
        tempfac = one(T) / √abs(um * a^2)
        Δτtemp = (τ % Ghat + (θo > T(π / 2) ? -one(T) : one(T)) * signβ * Go)
        n = floor(τ / Ghat)
        Δτ =
            (θo > T(π / 2) ? -one(T) : one(T)) * abs(
                argmin(
                    abs,
                    [(-one(T))^n * signβ * (Ghat - Δτtemp), (-one(T))^n * signβ * Δτtemp],
                ),
            )
    else
        argo = cos(θo) / √(up)
        k = m
        tempfac = inv(√abs(um * a^2))
        Δτtemp = (τ % Ghat + signβ * Go)
        n = floor(τ / Ghat)
        Δτ = argmin(
            abs,
            [(-one(T))^n * signβ * (Ghat - Δτtemp), (-one(T))^n * signβ * Δτtemp],
        )
    end

    return Δτ
end

@inline function Gϕ(pix::AbstractPixel, θs::T, isindir, n) where {T}
    α, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)

    signβ = sign(β)
    ηtemp = η(met, α, β, θo)
    λtemp = λ(met, α, θo)

    a = met.spin
    Go, Ghat = absGϕo_Gϕhat(pix)
    Gs, ans, isvortical = zero(T), zero(T), ηtemp < zero(T)

    isincone = abs(cos(θs)) < abs(cos(θo))
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return ans, Gs, Go, Ghat, isvortical, false
    end
    if (
        (((signβ < zero(T)) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical
    ) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return ans, Gs, Go, Ghat, isvortical, false
    end

    Δθ = (1 - (ηtemp + λtemp^2) / a^2) / T(2)
    up = min(Δθ + √(Δθ^2 + ηtemp / a^2), one(T) - eps(T))
    um = Δθ - √(Δθ^2 + ηtemp / a^2)
    m = up / um
    k = m

    #isvortical = η < 0.
    args = zero(T)
    argo = zero(T)
    #k = 0
    if isvortical
        args = (cos(θs)^2 - um) / (up - um)
        argo = (cos(θo)^2 - um) / (up - um)
        k = one(T) - m
        if (!(zero(T) < argo < one(T)) || !(zero(T) < args < one(T)))
            return ans, Gs, Go, Ghat, isvortical, false
        end
        tempfac = inv((1 - um) * √abs(um * a^2))
        argn = (up - um) / (1 - um)
        Go = ((θs > T(π / 2)) ? -1 : 1) * Go
        Gs = ((θs > T(π / 2)) ? -1 : 1) * tempfac * JacobiElliptic.Pi(argn, asin(√args), k)
    else
        args = cos(θs) / √(up)
        argo = cos(θo) / √(up)
        #k = abs(m)
        if !(-one(T) < args < one(T)) || !(-one(T) < argo < one(T))
            return ans, Gs, Go, Ghat, isvortical, false
        end
        tempfac = inv(√abs(um * a^2))
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
    α, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)
    signβ = sign(β)

    ηtemp = η(met, α, β, θo)
    λtemp = λ(met, α, θo)
    a = met.spin
    Go, Ghat = absGto_Gthat(pix)
    Gs, ans, isvortical = zero(T), zero(T), ηtemp < zero(T)

    isincone = abs(cos(θs)) < abs(cos(θo))
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return ans, Gs, Go, Ghat, isvortical, false
    end
    if (
        (((signβ < zero(T)) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical
    ) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return ans, Gs, Go, Ghat, isvortical, false
    end

    Δθ = (1 - (ηtemp + λtemp^2) / a^2) / T(2)
    up = min(Δθ + √(Δθ^2 + ηtemp / a^2), one(T) - eps(T))
    um = Δθ - √(Δθ^2 + ηtemp / a^2)
    m = up / um
    k = m

    #isvortical = η < 0.
    args = zero(T)
    argo = zero(T)
    #k = 0
    cosθo = cos(θo)
    if isvortical
        args = (cos(θs)^2 - um) / (up - um)
        argo = (cosθo^2 - um) / (up - um)
        k = one(T) - m
        if (!(zero(T) < argo < one(T)) || !(zero(T) < args < one(T)))
            return ans, Gs, Go, Ghat, isvortical, false
        end
        tempfac = √abs(um / a^2)
        Go *= ((θs > T(π / 2)) ? -1 : 1)
        Gs = ((θs > T(π / 2)) ? -1 : 1) * tempfac * JacobiElliptic.E(asin(√args), k)

    else
        args = cos(θs) / √(up)
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-one(T) < args < one(T)) || !(-one(T) < argo < one(T))
            return ans, Gs, Go, Ghat, isvortical, false
        end
        tempfac = -2 * up * inv(√abs(um * a^2))
        Gs =
            tempfac * (JacobiElliptic.E(asin(args), k) - JacobiElliptic.F(asin(args), k)) /
            (2k)
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    ans = (
        isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs :
        n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs
    ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return ans, Gs, Go, Ghat, isvortical, true
end

##----------------------------------------------------------------------------------------------------------------------
# SlowLightIntensityCachedPixel utility functions
##----------------------------------------------------------------------------------------------------------------------

@inline function _absGθo_Gθhat(metric::Kerr{T}, θo, η, λ)::NTuple{2,T} where {T}
    a = metric.spin
    a2 = a^2
    Go, Ghat, isvortical = zero(T), zero(T), η < zero(T)

    Δθ = (one(T) - (η + λ^2) / a2) / 2
    Δθ2 = Δθ^2
    desc = √max(Δθ2 + η / a2, 0)
    up = min(Δθ + desc, one(T) - eps(T))
    um = Δθ - desc
    m = up / um
    k = m

    argo = zero(T)
    k = zero(T)

    cosθo = cos(θo)
    if isvortical
        argo = (cosθo^2 - um) / (up - um)
        k = one(T) - m
        if (!(zero(T) < argo < one(T)))
            return Go, Ghat
        end
        tempfac = one(T) / √abs(um * a2)
        Go = tempfac * JacobiElliptic.F(asin(√argo), k)
        Ghat = 2 * tempfac * JacobiElliptic.K(k)
    else
        argo = cosθo / √(up)
        k = m
        if !(-one(T) < argo < one(T))
            return Go, Ghat
        end
        tempfac = one(T) / √abs(um * a^2)
        Go = tempfac * JacobiElliptic.F(asin(argo), k)
        Ghat = 2 * tempfac * JacobiElliptic.K(k)
    end

    return Go, Ghat
end

@inline function _absGϕo_Gϕhat(metric::Kerr{T}, θo, η, λ)::NTuple{2,T} where {T}

    a = metric.spin
    Go, Ghat, isvortical = zero(T), zero(T), η < zero(T)

    Δθ = (1 - (η + λ^2) / a^2) / T(2)
    up = min(Δθ + √(Δθ^2 + η / a^2), one(T) - eps(T))
    um = Δθ - √(Δθ^2 + η / a^2)
    m = up / um
    k = m

    #isvortical = η < 0.
    argo = zero(T)
    #k = 0
    cosθo = cos(θo)
    if isvortical
        argo = (cosθo^2 - um) / (up - um)
        k = one(T) - m
        if (!(zero(T) < argo < one(T)))
            return Go, Ghat
        end
        tempfac = inv((1 - um) * √abs(um * a^2))
        argn = (up - um) / (1 - um)
        Go = tempfac * JacobiElliptic.Pi(argn, asin(√argo), k)
        Ghat = 2tempfac * JacobiElliptic.Pi(argn, k)
    else
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-one(T) < argo < one(T))
            return Go, Ghat
        end
        tempfac = inv(√abs(um * a^2))

        Go = tempfac * JacobiElliptic.Pi(up, asin(argo), k)
        Ghat = 2tempfac * JacobiElliptic.Pi(up, k)
    end

    return Go, Ghat
end

@inline function _absGto_Gthat(metric::Kerr{T}, θo, η, λ)::NTuple{2,T} where {T}
    a = metric.spin
    Go, Ghat, isvortical = zero(T), zero(T), η < zero(T)

    Δθ = (1 - (η + λ^2) / a^2) / T(2)
    up = min(Δθ + √(Δθ^2 + η / a^2), one(T) - eps(T))
    um = Δθ - √(Δθ^2 + η / a^2)
    m = up / um
    k = m

    #isvortical = η < 0.
    argo = zero(T)
    #k = 0
    cosθo = cos(θo)
    if isvortical
        argo = (cosθo^2 - um) / (up - um)
        k = one(T) - m
        if (!(zero(T) < argo < one(T)))
            return Go, Ghat
        end
        tempfac = √abs(um / a^2)
        Go = tempfac * JacobiElliptic.E(asin(√argo), k)
        Ghat = 2tempfac * JacobiElliptic.E(k)

    else
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-one(T) < argo < one(T))
            return Go, Ghat
        end
        tempfac = -2 * up * inv(√abs(um * a^2))
        Go =
            tempfac * (JacobiElliptic.E(asin(argo), k) - JacobiElliptic.F(asin(argo), k)) /
            (2k)
        Ghat = tempfac * (JacobiElliptic.E(k) - JacobiElliptic.K(k)) / k
    end

    return Go, Ghat
end
