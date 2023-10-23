# Miscellaneous functions
##----------------------------------------------------------------------------------------------------------------------
export λ, η, α, β, αboundary, βboundary, r_potential, θ_potential, get_radial_roots, 
Ir,
Gθ

function _pow(z::Complex{T}, i) where {T}
    zabs = abs(z)
    zangle = angle(z)
    return (zabs^i) * (cos(zangle * i) + sin(zangle * i) * one(T)im)
end

function _pow(z::T, i) where {T<:Real}
    zabs = abs(z)
    if sign(z) < zero(T)
        return (zabs^i) * (cos(T(π) * i) + sin(T(π) * i)im)
    end
    return zabs^i + zero(T)im
end

"""
Checks if a complex number is real to some tolerance
"""
function _isreal2(num::Complex{T}) where T
    ren, imn = reim(num)
    return abs(imn / ren) < eps(T)^(T(1 / 5))
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
        JacobiElliptic.F(φ, j)
        -
        α^2 / (j + (one(T) - j) * α^2) * (
            JacobiElliptic.E(φ, j)
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

    return -inv((1 + α^2) * (1 - j + α^2)) * ((1 - j) * JacobiElliptic.F(φ, j) + α^2 * JacobiElliptic.E(φ, j))+# + α^2*√(1-j*sin(φ)^2)-α^3) +
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
    return p2 / 2 * log(abs((1 - p2) / (1 + p2) * (1 + p2 * √(1 - j * sinφ^2)) / (1 - p2 * √(1 - j * sinφ^2))))
end

function R1(α::T, φ::T, j) where {T}
    return one(T) / (one(T) - α^2) * (JacobiElliptic.Pi(α^2 / (α^2 - one(T)), φ, j) - α * f1(α, sin(φ), j))
end

function R2(α::T, φ::T, j) where {T}
    return one(T) / (one(T) - α^2) * (
        JacobiElliptic.F(φ, j)
        -
        α^2 / (j + (one(T) - j) * α^2) * (
            JacobiElliptic.E(φ, j) - α * sin(φ) * √(one(T) - j * sin(φ)^2) / (one(T) + α * cos(φ))
        )
    ) + inv(j + (one(T) - j) * α^2) * (2 * j - α^2 / (α^2 - one(T))) * R1(α, φ, j)
end

function S1(α, φ, j)
    α2 = α * α
    return inv(1 + α2) * (JacobiElliptic.F(φ, j) + α2 * JacobiElliptic.Pi(1 + α2, φ, j) - α * f2(α, sin(φ), j))
end

function S2(α, φ, j)
    return -inv((1 + α^2) * (1 - j + α^2)) * ((1 - j) * JacobiElliptic.F(φ, j) + α^2 * JacobiElliptic.E(φ, j) + α^2 * √(1 - j * sin(φ)^2) * (α - tan(φ)) / (1 + α * tan(φ)) - α^3) +
    (inv(1 + α^2) + (1 - j) / (1 - j + α^2)) * S1(α, φ, j)
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
Horizontal Bardeen Screen Coordinate

# Arguments

- `metric`: Kerr
- `λ`: Energy reduced Azimuthal angular momentul
- `η`: Energy reduced Carter integral 
- `θo`: Observer inclination
"""
function β(metric::Kerr, λ, η, θo)
    return sqrt(η - (α(metric, λ, θo)^2 - metric.spin^2) * cos(θo)^2)
end

"""
Defines a horizontal boundary on the assmyptotic observers screen where emission that originates from θs must fall within.

# Arguments

- `metric`: Kerr metric
- `θs`  : Emission Inclination
"""
function αboundary(metric::Kerr, θs)
    return metric.spin * sin(θs)
end

"""
Defines a vertical boundary on the Assyptotic observers screen where emission that originates from θs must fall within.

# Arguments

- `metric`: Kerr{T} metric
- `α`   : Horizontal Bardeen screen coordinate
- `θo`  : Observer inclination
- `θs`  : Emission Inclination
"""
function βboundary(metric::Kerr{T}, α, θo, θs) where {T}
    a = metric.spin
    cosθs2 = cos(θs)^2
    √max((cos(θo)^2 - cosθs2) * (α^2 - a^2 * (one(T) - cosθs2)) / (cosθs2 - one(T)), zero(T)) #eq 15 DOI 10.3847/1538-4357/acafe3 
end

"""
Radial potential of spacetime

# Arguments

- `metric`: Kerr{T} metric
- `η`  : Reduced Carter constant
- `λ`  : Reduced azimuthal agular momentum
- `r`  : Boyer Lindquist radius
"""
function r_potential(metric::Kerr{T}, η, λ, r) where {T}
    a = metric.spin
    λ2 = λ^2
    a * (a * (r * (r + 2) - η) - 4 * λ * r) + r * ((η + η) + (λ2 + λ2) + r * (-η - λ2 + r^2)) # Eq 7 PhysRevD.101.044032
end

"""
Theta potential of a Kerr blackhole

# Arguments

- `metric`: Kerr{T} metric
- `η`  : Reduced Carter constant
- `λ`  : Reduced azimuthal agular momentum
- `θ`  : Boyer Lindquist inclination
"""
function θ_potential(metric::Kerr{T}, η, λ, θ) where {T}
    a = metric.spin
    η + a^2 * cos(θ)^2 - λ^2 * cot(θ)^2
end

##----------------------------------------------------------------------------------------------------------------------
# Radial functions
##----------------------------------------------------------------------------------------------------------------------
"""
Returns roots of \$r^4 + (a^2-η-λ^2)r^2 + 2(η+(a-λ)^2)r - a^2η\$

# Arguments

- `metric`: Kerr{T} metric
- `η`  : Reduced Carter constant
- `λ`  : Reduced azimuthal agular momentum
"""
function get_radial_roots(metric::Kerr{T}, η, λ) where {T}
    a = metric.spin

    a2 = a * a
    A = a2 - η - λ * λ
    A2 = A + A
    B = T(2) * (η + (λ - a)^2)
    C = -a2 * η

    P = -A * A / T(12) - C
    Q = -A / T(3) * (A * A / T(36) + zero(T)im - C) - B * B / T(8)

    Δ3 = -T(4) * P * P * P - T(27) * Q * Q
    ωp = _pow(-Q / T(2) + _pow(-Δ3 / T(108), T(0.5)) + zero(T)im, T(1 / 3))

    #C = ((-1+0im)^(2/3), (-1+0im)^(4/3), 1) .* ωp
    C = (T(-0.4999999999999999) + T(0.8660254037844387)im, T(-0.5000000000000002) - T(0.8660254037844385)im, one(T) + zero(T)im) .* ωp

    v = -P .* _pow.(T(3) .* C, -one(T))

    ξ0 = argmax(real, (C .+ v)) - A / T(3)
    ξ02 = ξ0 + ξ0

    predet1 = A2 + ξ02
    predet2 = (√T(2) * B) * _pow(ξ0, T(-0.5))
    det1 = _pow(-(predet1 - predet2), T(0.5))
    det2 = _pow(-(predet1 + predet2), T(0.5))

    sqrtξ02 = _pow(ξ02, T(0.5))

    r1 = (-sqrtξ02 - det1) / T(2)
    r2 = (-sqrtξ02 + det1) / T(2)
    r3 = (sqrtξ02 - det2) / T(2)
    r4 = (sqrtξ02 + det2) / T(2)

    return r1, r2, r3, r4
end

_get_root_diffs(r1, r2, r3, r4) = r2 - r1, r3 - r1, r3 - r2, r4 - r1, r4 - r2, r4 - r3

"""
Returns the antiderivative \$I_r=\\int\\frac{dr}{\\sqrt{\\mathcal{R(r)}}}\$

# Arguments

- `metric`: Kerr{T} metric
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
- `θo` : Observer inclination
- `α`  : Horizontal Bardeen screen coordinate
- `β`  : Vertical Bardeen screen coordinate
"""
function Ir(metric::Kerr, νr::Bool, θo, rs, α, β)
    return Ir(metric, νr, rs, η(metric, α, β, θo), λ(metric, α, θo))
end
function Ir(metric::Kerr{T}, νr::Bool, rs, η, λ) where {T}
    roots = get_radial_roots(metric, η, λ)
    root_diffs = _get_root_diffs(roots...)
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return Ir_case1_and_2(metric, real.(roots), rs, νr)[1]
    elseif numreals == 2 #case3
        if abs(imag(roots[4])) < T(1e-10)
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        return Ir_case3(metric, roots, root_diffs, rs)
    else #case 4
        return Ir_case4(metric, roots, root_diffs, rs)
    end
    return T(NaN), T(NaN)
end

function Ir_case1_and_2(::Kerr{T}, roots::NTuple{4}, rs, νr) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    k = r32 * r41 / (r31 * r42)
    x2_s2 = ((rs - r4) / (rs - r3) * r31 / r41)
    x2_s = √abs(x2_s2)
    coef = 2 / √real(r31 * r42)
    Ir_s = !(-one(T) < x2_s2 < one(T)) ? T(NaN) : coef * JacobiElliptic.F(asin(x2_s), k)
    Ir_inf = coef * JacobiElliptic.F(asin(√(r31 / r41)), k)

    return Ir_inf - (νr ? Ir_s : -Ir_s), Ir_inf
end

function Ir_case3(::Kerr{T}, roots::NTuple{4}, root_diffs::NTuple{6}, rs) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = root_diffs

    r1, r2, r21 = real.((r1, r2, r21))

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    #if A2 < zero(0.0) || B2 < zero(T)
    #  return zero(T)
    #end
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    temprat = B * (rs - r2) / (A * (rs - r1))
    x3_s = ((one(T) - temprat) / (one(T) + temprat))
    coef = one(T) * √inv(A * B)
    Ir_s = coef * JacobiElliptic.F((acos(x3_s)), k3)
    Ir_inf = coef * JacobiElliptic.F((acos(((A - B) / (A + B)))), k3)

    return Ir_inf - Ir_s, Ir_inf
end

function Ir_case4(::Kerr{T}, roots::NTuple{4}, root_diffs::NTuple{6}, rs) where {T}
    _, r1, _, r4 = roots
    _, r31, r32, r41, r42, _ = root_diffs

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return zero(T)
    end
    C = √real(r31 * r42)
    D = √real(r32 * r41)
    k4 = 4C * D / (C + D)^2
    a2 = abs(imag(r1))
    b1 = real(r4)

    k4 = T(4) * C * D / (C + D)^2

    go = √max((T(4)a2^2 - (C - D)^2) / ((C + D)^2 - T(4)a2^2), zero(T))
    x4_s = (rs + b1) / a2
    coef = 2 / (C + D)
    Ir_s_coef = JacobiElliptic.F(atan(x4_s) + atan(go), k4)
    Ir_inf_coef = JacobiElliptic.F(T(π / 2) + atan(go), k4)
    #return (C+D), coef
    return coef.*(Ir_inf_coef - Ir_s_coef, Ir_inf_coef)
end

function Iϕ_case2(metric::Kerr{T}, roots::NTuple{4}, root_diffs::NTuple{6}, rs, τ, νr, λ) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, r43 = root_diffs
    a = metric.spin
    a2 = a * a
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)
    rp3 = rp - r3
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4

    k = r32 * r41 / (r31 * r42)
    x2_s = √((rs - r4) / (rs - r3) * r31 / r41)
    x2_o = √(r31 / r41)
    !(-1 < x2_s < 1) && return zero(T)

    #Fo = isnan(τo) ? 2 / √real(r31 * r42) * JacobiElliptic.F(asin(x2_o), k) : τo
    #Fs = 2 / √real(r31 * r42) * JacobiElliptic.F(asin(x2_s), k)

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_s), k)
    Πp_o = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)
    Πm_s = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_s), k)
    Πm_o = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_o), k)

    Ip = -Πp_o - τ / rp3
    Im = -Πm_o - τ / rm3

    if νr
        Ip += Πp_s
        Im += Πm_s
    else
        Ip -= Πp_s
        Im -= Πm_s
    end

    return 2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
end

function Iϕ_case3(metric::Kerr{T}, roots::NTuple{4}, root_diffs::NTuple{6}, rs, τ, λ) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = root_diffs
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
        return zero(T)
    end
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)

    temprat = B * (rs - r2) * real(_pow(A * (rs - r1), -one(T)))
    x3_s = clamp(((one(T) - temprat) * real(_pow(one(T) + temprat, -one(T)))), -one(T), one(T))
    x3_o = clamp((A - B) / (A + B), -one(T), one(T))
    φ_s = acos(x3_s)
    φ_o = acos(x3_o)

    (isnan(φ_s) || isnan(φ_o)) && return T(NaN)

    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

    R1p_s = R1(αp, φ_s, k3)
    R1p_o = R1(αp, φ_o, k3)
    R1m_s = R1(αm, φ_s, k3)
    R1m_o = R1(αm, φ_o, k3)

    Ip = -inv(B * rp2 + A * rp1) * ((B + A) * τ + 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (R1p_o - R1p_s))
    Im = -inv(B * rm2 + A * rm1) * ((B + A) * τ + 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (R1m_o - R1m_s))

    return (2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im))
end

function Iϕ_case4(metric::Kerr{T}, roots::NTuple{4}, root_diffs::NTuple{6}, rs, τ, λ) where {T}
    r1, _, _, r4 = roots
    _, r31, r32, r41, r42, _ = root_diffs
    a = metric.spin
    a2 = a * a
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return zero(T)
    end

    C = √real(r31 * r42)
    D = √real(r32 * r41)
    k4 = T(4) * C * D / (C + D)^2
    a2 = abs(imag(r1))
    b1 = real(r4)

    k4 = 4 * C * D / (C + D)^2

    x4_s = (rs + b1) / a2
    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_s) || isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    S1p_s = S1(gp, atan(x4_s) + atan(go), k4)
    S1p_o = S1(gp, T(π / 2) + atan(go), k4)
    S1m_s = S1(gm, atan(x4_s) + atan(go), k4)
    S1m_o = S1(gm, T(π / 2) + atan(go), k4)

    Ip = go / (a2 * (1 - go * x4_p)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1p_o - S1p_s))
    Im = go / (a2 * (1 - go * x4_m)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1m_o - S1m_s))

    return 2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
end

function It_case2(metric::Kerr{T}, roots::NTuple{4}, root_diffs::NTuple{6}, rs, τ, νr, λ) where {T}
    r1, r2, r3, r4 = roots
    _, r31, r32, r41, r42, _ = root_diffs
    r43 = r4 - r3
    a = metric.spin
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)
    rp3 = rp - r3
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4

    k = r32 * r41 / (r31 * r42)
    x2_s = √abs((rs - r4) / (rs - r3) * r31 / r41)
    x2_o = √abs(r31 / r41)
    !(-1 < x2_s < 1) && return zero(T)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    E_s = √(r31 * r42) * JacobiElliptic.E(asin(x2_s), k)
    E_o = √(r31 * r42) * JacobiElliptic.E(asin(x2_o), k)
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
    I1_total = r3 * I0_total + log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k) + (νr ? -1 : 1) * Π1_s)# Removed the logarithmic divergence
    #equation B38
    I2_s = √(evalpoly(rs, poly_coefs)) / (rs - r3) - E_s
    I2_total = r3 - (r1 * r4 + r2 * r3) / 2 * τ - E_o + (νr ? -1 : 1) * I2_s# Assymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_s), k)
    Πp_o = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)
    Πm_s = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_s), k)
    Πm_o = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_o), k)

    Ip_total = -Πp_o - τ / rp3
    Im_total = -Πm_o - τ / rm3

    if νr
        Ip_total += Πp_s
        Im_total += Πm_s
    else
        Ip_total -= Πp_s
        Im_total -= Πm_s
    end

    return -(4 / (rp - rm) * (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) + 4 * I0_total + 2 * I1_total + I2_total)
end

function It_case3(metric::Kerr{T}, roots::NTuple{4}, root_diffs::NTuple{6}, rs, τ, λ) where {T}
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = root_diffs
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
        return zero(T)
    end
    A, B = √A2, √B2

    k3 = real(((A + B)^2 - r21^2) / (4 * A * B))

    temprat = B * (rs - r2) * _pow(A * (rs - r1), -one(T))
    x3_s = clamp(real(((one(T) - temprat) * _pow(one(T) + temprat, -one(T)))), -one(T), one(T))
    x3_o = min((A - B) / (A + B), one(T))
    φ_s = acos(x3_s)
    φ_o = acos(x3_o)

    (isnan(φ_s) || isnan(φ_o)) && return T(NaN)

    αo = (B + A) / (B - A)
    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

    Π1_s = 2 * r21 * √real(A * B) / (B2 - A2) * R1(αo, φ_s, k3)
    Π1_o = 2 * r21 * √real(A * B) / (B2 - A2) * regularized_R1(αo, φ_o, k3) # Divergence is removed, will be added back in the end
    Π2_s = ((2 * r21 * √(A * B) / (B2 - A2))^2) * R2(αo, φ_s, k3)
    Π2_o = ((2 * r21 * √(A * B) / (B2 - A2))^2) * regularized_R2(αo, φ_o, k3)# Divergence is removed, will be added back in the end

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = (B * r2 + A * r1) / (B + A) * I0_total + Π1_o - Π1_s + log(16 * r21^2 / ((A2 - B2)^2 + 4 * A * B * r21^2)) / 2 
    # Removed linear divergence
    I2_total = ((((B * r2 + A * r1) / (B + A))^2) * I0_total - 2 * (B * r2 + A * r1) / (B + A) * (-Π1_s) - √(A * B) * (Π2_o - Π2_s)) + (B * r2 + A * r1) / (A + B)
    Ip_total = -inv(B * rp2 + A * rp1) * ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (R1(αp, φ_o, k3) - R1(αp, φ_s, k3)))
    Im_total = -inv(B * rm2 + A * rm1) * ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (R1(αm, φ_o, k3) - R1(αm, φ_s, k3)))

    return -(4 / (rp - rm) * (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) + 4 * I0_total + 2 * I1_total + I2_total)
end

function It_case4(metric::Kerr{T}, roots::NTuple{4}, root_diffs::NTuple{6}, rs, τ, λ) where {T}
    a = metric.spin
    r1, _, _, r4 = roots
    _, r31, r32, r41, r42, _ = root_diffs
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return zero(T)
    end

    C = √real(r31 * r42)
    D = √real(r32 * r41)
    k4 = T(4) * C * D / (C + D)^2
    a2 = abs(imag(r1))
    b1 = real(r4)

    k4 = 4 * C * D / (C + D)^2

    x4_s = (rs + b1) / a2
    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_s) || isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    S1p_s = S1(gp, atan(x4_s) + atan(go), k4)
    S1p_o = S1(gp, T(π / 2) + atan(go), k4)
    S1m_s = S1(gm, atan(x4_s) + atan(go), k4)
    S1m_o = S1(gm, T(π / 2) + atan(go), k4)


    Π1_s = 2 / (C + D) * (a2 / go * (1 + go^2)) * S1(go, atan(x4_s) + atan(go), k4)
    Π1_o = 2 / (C + D) * (a2 / go * (1 + go^2)) * regularizedS1(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end
    Π2_s = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * S2(go, atan(x4_s) + atan(go), k4)
    Π2_o = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * regularizedS2(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = (a2 / go - b1) * τ - (Π1_o - Π1_s) + 1/2*log(
        (16*(1 + go^2 - sqrt((1 + go^2)*(1 + go^2 - k4)))*(1 + go^2 - 
   k4))/((C + D)^2*((1 + go^2)^2 - k4)*k4*(1 + sqrt(
   1 - k4/(1 + go^2)))))

    # Removed linear divergence
    I2_total =
        ((a2 / go - b1)^2) * τ - 2(a2 / go - b1) * (Π1_o - Π1_s) + (Π2_o - Π2_s) -
        ((16 * a2^4 + (C^2 - D^2)^2 - 8 * (a2^2) * (C^2 + D^2) + 8 * (a2^3) * (C + D - 2 * b1 * go) + 2 * a2 * (C + D) * (-(C - D)^2 + 2 * b1 * (C + D) * go)) / (4 * a2 * (4 * a2^2 - (C + D)^2) * go))
    Ip_total = go / (a2 * (1 - go * x4_p)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1p_o - S1p_s))
    Im_total = go / (a2 * (1 - go * x4_m)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1m_o - S1m_s))


    #logdiv = zero(T)
    #lineardiv = -((16 * a2^4 + (C^2 - D^2)^2 - 8 * (a2^2) * (C^2 + D^2) + 8 * (a2^3) * (C + D - 2 * b1 * go) + 2 * a2 * (C + D) * (-(C - D)^2 + 2 * b1 * (C + D) * go)) / (4 * a2 * (4 * a2^2 - (C + D)^2) * go))

    return -(4 / (rp - rm) * (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) + 4 * I0_total + 2 * I1_total + I2_total)# + (logdiv + lineardiv)
end

"""
Returns the radial integrals for the case where there are four real roots in the radial potential, with roots outside the horizon.
"""
function radial_integrals_case2(metric::Kerr{T}, rs, roots::NTuple{4}, τ, νr) where {T}
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

    k = r32 * r41 / (r31 * r42)
    x2_s = √abs((rs - r4) / (rs - r3) * r31 / r41)
    x2_o = √abs(r31 / r41)
    !(-1 < x2_s < 1) && return T(NaN), T(NaN), T(NaN), T(NaN), T(NaN)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    E_s = √(r31 * r42) * JacobiElliptic.E(asin(x2_s), k)
    E_o = √(r31 * r42) * JacobiElliptic.E(asin(x2_o), k)
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
    I1_total = r3 * I0_total + log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k) + (νr ? -1 : 1) * Π1_s)# Removed the logarithmic divergence
    #equation B38
    I2_s = √abs(evalpoly(rs, poly_coefs)) / (rs - r3) - E_s
    I2_total = r3 - (r1 * r4 + r2 * r3) / 2 * τ - E_o + (νr ? -1 : 1) * I2_s# Assymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_s), k)
    Πp_o = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)
    Πm_s = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_s), k)
    Πm_o = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_o), k)

    Ip_total = -Πp_o - τ / rp3
    Im_total = -Πm_o - τ / rm3

    if νr
        Ip_total += Πp_s
        Im_total += Πm_s
    else
        Ip_total -= Πp_s
        Im_total -= Πm_s
    end

    return I0_total, I1_total, I2_total, Ip_total, Im_total

end

"""
Returns the radial integrals for the case where there are two real roots in the radial potential
"""
function radial_integrals_case3(metric::Kerr{T}, rs, roots::NTuple{4}, τ) where {T}
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
    temprat = B * (rs - r2) * _pow(A * (rs - r1), -one(T))
    x3_s = real((one(T) - temprat) * _pow(one(T) + temprat, -one(T)))

    abs(x3_s) > one(T) && return T(NaN), T(NaN), T(NaN), T(NaN), T(NaN)

    φ_s = acos(x3_s)
    αo = (B + A) / (B - A)
    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)
    coef = 2 * r21 * √(A * B) / (B2 - A2)

    Π1_s = coef * R1(αo, φ_s, k3)
    Π2_s = (coef^2) * R2(αo, φ_s, k3)

    φ_o = acos((A - B) / (A + B))
    Π1_o = 2 * r21 * √real(A * B) / (B2 - A2) * regularized_R1(αo, φ_o, k3)
    Π2_o = ((2 * r21 * √(A * B) / (B2 - A2))^2) * regularized_R2(αo, φ_o, k3)

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = (B * r2 + A * r1) / (B + A) * I0_total + Π1_o - Π1_s + log(16 * r21^2 / ((A2 - B2)^2 + 4 * A * B * r21^2)) / 2
    # Removed linear divergence
    I2_total = ((((B * r2 + A * r1) / (B + A))^2) * I0_total - 2 * (B * r2 + A * r1) / (B + A) * (-Π1_s) - √(A * B) * (Π2_o - Π2_s)) + (B * r2 + A * r1) / (A + B)
    Ip_total = -inv(B * rp2 + A * rp1) * ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (R1(αp, φ_o, k3) - R1(αp, φ_s, k3)))
    Im_total = -inv(B * rm2 + A * rm1) * ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (R1(αm, φ_o, k3) - R1(αm, φ_s, k3)))

    return I0_total, I1_total, I2_total, Ip_total, Im_total
end
"""
Returns the radial integrals for the case where there are no real roots in the radial potential
"""
function radial_integrals_case4(metric::Kerr{T}, rs, roots::NTuple{4}, τ) where {T}
    r1, _, _, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    a = metric.spin
    a2 = a * a
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)

    C = √real(r31 * r42)
    D = √real(r32 * r41)
    a2 = abs(imag(r1))
    b1 = real(r4)
    k4 = 4 * C * D / (C + D)^2

    x4_s = (rs + b1) / a2
    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_s) || isnan(x4_p) || isnan(x4_m)) && return T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN)
    S1p_s = S1(gp, atan(x4_s) + atan(go), k4)
    S1m_s = S1(gm, atan(x4_s) + atan(go), k4)
    #Fr_s = 2 / (C + D) * JacobiElliptic.F(atan(x4_s) + atan(go), k4)
    Π1_s = 2 / (C + D) * (a2 / go * (1 + go^2)) * S1(go, atan(x4_s) + atan(go), k4)
    Π2_s = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * S2(go, atan(x4_s) + atan(go), k4)

    S1p_o = S1(gp, T(π / 2) + atan(go), k4)
    S1m_o = S1(gm, T(π / 2) + atan(go), k4)
    #Fr_o = 2 / (C + D) * JacobiElliptic.F(T(π/2) + atan(go), k4)
    Π1_o = 2 / (C + D) * (a2 / go * (1 + go^2)) * regularizedS1(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end
    Π2_o = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * regularizedS2(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = (a2 / go - b1) * τ - (Π1_o - Π1_s) + 1/2*log(
        (16*(1 + go^2 - sqrt((1 + go^2)*(1 + go^2 - k4)))*(1 + go^2 - 
   k4))/((C + D)^2*((1 + go^2)^2 - k4)*k4*(1 + sqrt(
   1 - k4/(1 + go^2)))))

    # Removed linear divergence
    I2_total =
        ((a2 / go - b1)^2) * τ - 2(a2 / go - b1) * (Π1_o - Π1_s) + (Π2_o - Π2_s) -
        ((16 * a2^4 + (C^2 - D^2)^2 - 8 * (a2^2) * (C^2 + D^2) + 8 * (a2^3) * (C + D - 2 * b1 * go) + 2 * a2 * (C + D) * (-(C - D)^2 + 2 * b1 * (C + D) * go)) / (4 * a2 * (4 * a2^2 - (C + D)^2) * go))
    Ip_total = go / (a2 * (1 - go * x4_p)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1p_o - S1p_s))
    Im_total = go / (a2 * (1 - go * x4_m)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1m_o - S1m_s))

    return I0_total, I1_total, I2_total, Ip_total, Im_total
end

##----------------------------------------------------------------------------------------------------------------------
# Inclination functions
##----------------------------------------------------------------------------------------------------------------------
"""
Mino time of trajectory between two inclinations for a given screen coordinate

# Arguments 

- `α` : Horizontal Bardeen screen coordinate
- `β` : Vertical Bardeen screen coordinate 
- `a` : Blackhole angular Momentum
- `θs` : Emission inclination
- `θo` : Observer inclination
- `isindir` : Is the path direct or indirect?
- `n` : nth image in orde of amount of minotime traversed
"""
function mino_time(metric::Kerr{T}, α, β, θs, θo, isindir, n) where {T}
    return Gθ(metric, α, β, θs, θo, isindir, n)[1]
end

"""
Returns the antiderivative \$G_\\theta=\\int\\frac{d\\theta}{\\sqrt{\\Theta(\\theta)}}\$

# Arguments 

- `metric`: Kerr{T} metric
- `α` : Horizontal Bardeen screen coordinate
- `β` : Vertical Bardeen screen coordinate
- `θs` : Emission inclination
- `θo` : Observer inclination
- `isindir` : Is the path direct or indirect?
- `n` : nth image ordered by minotime
"""
function Gθ(metric::Kerr{T}, α, β, θs, θo, isindir, n) where {T}
    return _Gθ(metric::Kerr{T}, sign(β), θs, θo, isindir, n, η(metric, α, β, θo), λ(metric, α, θo))
end

function _Gθ(metric::Kerr{T}, signβ, θs, θo, isindir, n, η, λ) where {T}
    a = metric.spin
    a2 = a^2
    Go, Gs, Ghat, minotime, isvortical = zero(T), zero(T), zero(T), zero(T), η < zero(T)

    isincone = abs(cos(θs)) < abs(cos(θo))
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
    end
    if ((((signβ < 0) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
    end

    Δθ = (one(T) - (η + λ^2) / a2) / 2
    Δθ2 = Δθ^2
    desc = √(Δθ2 + η / a2)
    up = Δθ + desc
    um = Δθ - desc
    m = up / um
    k = m

    #isvortical = η < 0.
    args = zero(T)
    argo = zero(T)
    k = zero(T)
    if isvortical
        args = (cos(θs)^2 - um) / (up - um)
        argo = (cos(θo)^2 - um) / (up - um)
        k = one(T) - m
        if (!(zero(T) < argo < one(T)) || !(zero(T) < args < one(T)))
            return T(NaN), Gs, Go, Ghat, isvortical
        end
        tempfac = one(T) / √abs(um * a2)
        Go = ((θs > T(π / 2)) ? -one(T) : one(T)) * tempfac * JacobiElliptic.F(asin(√argo), k)
        Gs = ((θs > T(π / 2)) ? -one(T) : one(T)) * tempfac * JacobiElliptic.F(asin(√args), k)
        Ghat = 2 * tempfac * JacobiElliptic.K(k)
    else
        args = cos(θs) / √(up)
        argo = cos(θo) / √(up)
        k = m
        if !(-one(T) < args < one(T)) || !(-one(T) < argo < one(T))
            return T(NaN), Gs, Go, Ghat, isvortical
        end
        tempfac = one(T) / √abs(um * a^2)
        Go = tempfac * JacobiElliptic.F(asin(argo), k)
        Gs = tempfac * JacobiElliptic.F(asin(args), k)
        Ghat = 2 * tempfac * JacobiElliptic.K(k)
    end


    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    minotime = (isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs : n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return minotime, Gs, Go, Ghat, isvortical
end

function Gs(metric::Kerr{T}, α, β, θo, τ) where {T}
    return _Gs(metric, sign(β), θo, η(metric, α, β, θo), λ(metric, α, θo), τ)
end

function _Gs(metric::Kerr{T}, signβ, θo, η, λ, τ) where {T}
    τ == T(Inf) && return T(Inf)
    a = metric.spin
    Go, Ghat, Ghat_2, isvortical = zero(T), zero(T), zero(T), η < zero(T)

    Δθ = T(0.5) * (one(T) - (η + λ^2) / a^2)
    up = Δθ + √(Δθ^2 + η / a^2)
    um = Δθ - √(Δθ^2 + η / a^2)
    m = up / um
    k = m

    argo = zero(T)
    k = zero(T)
    tempfac = one(T) / √abs(um * a^2)

    if isvortical
        argo = (cos(θo)^2 - um) / (up - um)
        k = one(T) - m
        tempfac = one(T) / √abs(um * a^2)
        Go = tempfac * JacobiElliptic.F(asin(√argo), k)
        Ghat_2 = tempfac * JacobiElliptic.K(k)
        Ghat = 2Ghat_2
        Δτtemp = (τ % Ghat + (θo > T(π / 2) ? -one(T) : one(T)) * signβ * Go)
        n = floor(τ / Ghat)
        Δτ = (θo > T(π / 2) ? -one(T) : one(T)) * abs(argmin(abs, [(-one(T))^n * signβ * (Ghat - Δτtemp), (-one(T))^n * signβ * Δτtemp]))
    else
        argo = cos(θo) / √(up)
        k = m
        tempfac = inv(√abs(um * a^2))
        Go = tempfac * JacobiElliptic.F(asin(argo), k)
        Ghat_2 = tempfac * JacobiElliptic.K(k)
        Ghat = Ghat_2 + Ghat_2
        Δτtemp = (τ % Ghat + signβ * Go)
        n = floor(τ / Ghat)
        Δτ = argmin(abs, [(-one(T))^n * signβ * (Ghat - Δτtemp), (-one(T))^n * signβ * Δτtemp])
    end

    return Δτ
end

function Gϕ(metric::Kerr{T}, α, β, θs, θo, isindir, n) where {T}
    return _Gϕ(metric::Kerr{T}, sign(β), θs, θo, isindir, n, η(metric, α, β, θo), λ(metric, α, θo))
end

function _Gϕ(metric::Kerr{T}, signβ, θs, θo, isindir, n, η, λ) where {T}
    a = metric.spin
    Go, Gs, Ghat, ans, isvortical = zero(T), zero(T), zero(T), zero(T), η < zero(T)

    isincone = abs(cos(θs)) < abs(cos(θo))
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), isvortical
    end
    if ((((signβ < zero(T)) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), isvortical
    end

    Δθ = (1 - (η + λ^2) / a^2) / T(2)
    up = Δθ + √(Δθ^2 + η / a^2)
    um = Δθ - √(Δθ^2 + η / a^2)
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
            return T(NaN), Gs, Ghat, isvortical
        end
        tempfac = inv((1 - um) * √abs(um * a^2))
        argn = (up - um) / (1 - um)
        Go = ((θs > T(π / 2)) ? -1 : 1) * tempfac * JacobiElliptic.Pi(argn, asin(√argo), k)
        Gs = ((θs > T(π / 2)) ? -1 : 1) * tempfac * JacobiElliptic.Pi(argn, asin(√args), k)
        Ghat = 2tempfac * JacobiElliptic.Pi(argn, k)
    else
        args = cos(θs) / √(up)
        argo = cos(θo) / √(up)
        #k = abs(m)
        if !(-one(T) < args < one(T)) || !(-one(T) < argo < one(T))
            return T(Inf), Gs, Ghat, isvortical
        end
        tempfac = inv(√abs(um * a^2))
        Go = tempfac * JacobiElliptic.Pi(up, asin(argo), k)
        Gs = tempfac * JacobiElliptic.Pi(up, asin(args), k)
        Ghat = 2tempfac * JacobiElliptic.Pi(up, k)
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    ans = (isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs : n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return ans, Gs, Ghat, isvortical

end

function Gt(metric::Kerr{T}, α, β, θs, θo, isindir, n) where {T}
    return _Gt(metric::Kerr{T}, sign(β), θs, θo, isindir, n, η(metric, α, β, θo), λ(metric, α, θo))
end

function _Gt(metric::Kerr{T}, signβ, θs, θo, isindir, n, η, λ) where {T}
    a = metric.spin
    Go, Gs, Ghat, minotime, isvortical = zero(T), zero(T), zero(T), zero(T), η < zero(T)

    isincone = abs(cos(θs)) < abs(cos(θo))
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return T(Inf), T(Inf), T(Inf), isvortical
    end
    if ((((signβ < zero(T)) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return T(Inf), T(Inf), T(Inf), isvortical
    end

    Δθ = (1 - (η + λ^2) / a^2) / T(2)
    up = Δθ + √(Δθ^2 + η / a^2)
    um = Δθ - √(Δθ^2 + η / a^2)
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
            return T(Inf), Gs, Ghat, isvortical
        end
        tempfac = √abs(um / a^2)
        argn = (up - um) / (1 - um)
        Go = ((θs > T(π / 2)) ? -1 : 1) * tempfac * JacobiElliptic.E(asin(√argo), k)
        Gs = ((θs > T(π / 2)) ? -1 : 1) * tempfac * JacobiElliptic.E(asin(√args), k)
        Ghat = 2tempfac * JacobiElliptic.E(k)

    else
        args = cos(θs) / √(up)
        argo = cos(θo) / √(up)
        #k = abs(m)
        if !(-one(T) < args < one(T)) || !(-one(T) < argo < one(T))
            return T(Inf), Gs, Ghat, isvortical
        end
        tempfac = -2 * up * inv(√abs(um * a^2))
        Go = tempfac * (JacobiElliptic.E(asin(argo), k) - JacobiElliptic.F(asin(argo), k)) / (2k)
        Gs = tempfac * (JacobiElliptic.E(asin(args), k) - JacobiElliptic.F(asin(argo), k)) / (2k)
        Ghat = 2tempfac * (JacobiElliptic.E(k) - JacobiElliptic.K(k)) / k
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    minotime = real(isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs : n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return minotime, Gs, Ghat, isvortical
end
