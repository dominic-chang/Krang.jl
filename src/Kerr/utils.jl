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
    return abs(imn / ren) < sqrt(eps(T))
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
- `λ`  : Reduced azimuthal angular momentum
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
- `λ`  : Reduced azimuthal angular momentum
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
- `λ`  : Reduced azimuthal angular momentum
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

function _get_root_diffs(r1::T, r2::T, r3::T, r4::T) where T 
    return r2 - r1, r3 - r1, r3 - r2, r4 - r1, r4 - r2, r4 - r3
end

"""
Mino time of trajectory between an observer at infinity and point at radius rs

# Arguments 

- `metric`: Kerr{T} metric
- `pix` : Pixel information
- `rs` : Emission radius
- `isindir` : Is the path direct or indirect?
"""
function mino_time(metric::Kerr{T}, pix, rs, isindir) where {T}
    return Ir(metric, isindir, rs, pix)[1]
end

"""
Returns the antiderivative \$I_r=\\int\\frac{dr}{\\sqrt{\\mathcal{R(r)}}}\$.
See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
- `rs` : Emission radius
- `pix`  : Pixel information
"""
function Ir(metric::Kerr{T}, νr::Bool, rs, pix::BasicPixel) where T
    roots = pix.roots
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return Ir_case1_and_2(metric, pix, rs, νr)[1]
    elseif numreals == 2 #case3
        if abs(imag(roots[4])) < T(1e-10)
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        return Ir_case3(metric, pix, rs)
    else #case 4
        return Ir_case4(metric, pix, rs)
    end
    return T(NaN), T(NaN)
end

function Ir_case1_and_2(::Kerr{T}, pix::BasicPixel, rs, νr) where {T} 
    roots = real.(pix.roots)
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    k = r32 * r41 / (r31 * r42)
    x2_s2 = ((rs - r4) / (rs - r3) * r31 / r41)
    x2_s = √abs(x2_s2)
    coef = 2 / √real(r31 * r42)
    Ir_s = !(-one(T) < x2_s2 < one(T)) ? T(NaN) : coef * JacobiElliptic.F(asin(x2_s), k)
    #Ir_inf = coef * JacobiElliptic.F(asin(√(r31 / r41)), k)
    Ir_inf = pix.Ir

    return Ir_inf - (νr ? Ir_s : -Ir_s), Ir_inf
end

function Ir_case3(::Kerr{T}, pix::BasicPixel, rs) where {T}
    roots = pix.roots
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    r1, r2, r21 = real.((r1, r2, r21))

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    temprat = B * (rs - r2) / (A * (rs - r1))
    x3_s = ((one(T) - temprat) / (one(T) + temprat))
    coef = one(T) * √inv(A * B)
    Ir_s = coef * JacobiElliptic.F((acos(x3_s)), k3)
    Ir_inf = pix.Ir

    return Ir_inf - Ir_s, Ir_inf
end

function Ir_case4(::Kerr{T}, pix::BasicPixel, rs) where {T}
    roots = pix.roots
    _, r1, _, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return T(NaN)
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
    #Ir_inf_coef = JacobiElliptic.F(T(π / 2) + atan(go), k4)
    Ir_inf = pix.Ir
    #return (C+D), coef
    #return coef.*(Ir_inf_coef - Ir_s_coef, Ir_inf_coef)
    return Ir_inf - coef*Ir_s_coef, Ir_inf
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
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return Ir_inf_case1_and_2(metric, real.(roots))[1]
    elseif numreals == 2 #case3
        if abs(imag(roots[4])) < T(1e-10)
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        return Ir_inf_case3(metric, roots)
    else #case 4
        return Ir_inf_case4(metric, roots)
    end
    return T(NaN)
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

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    coef = one(T) * √inv(A * B)
    Ir_inf = coef * JacobiElliptic.F((acos(clamp((A - B) / (A + B), -one(T), one(T)))), k3)

    return Ir_inf
end

function Ir_inf_case4(::Kerr{T}, roots::NTuple{4}) where {T}
    _, r1, _, _ = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return T(NaN)
    end
    C = √real(r31 * r42)
    D = √real(r32 * r41)
    k4 = 4C * D / (C + D)^2
    a2 = abs(imag(r1))

    k4 = T(4) * C * D / (C + D)^2

    go = √max((T(4)a2^2 - (C - D)^2) / ((C + D)^2 - T(4)a2^2), zero(T))
    coef = 2 / (C + D)

    Ir_inf = coef*JacobiElliptic.F(T(π / 2) + atan(go), k4)
    #return (C+D), coef
    return Ir_inf
end

"""
Returns the antiderivative \$I_r=\\int\\frac{dr}{\\sqrt{\\mathcal{R(r)}}}\$ evaluated at infinity.
See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `roots`  : Roots of the radial potential
"""
function Ir_o(metric::Kerr{T}, roots) where {T}
    #root_diffs = _get_root_diffs(roots...)
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return Ir_o_case1_and_2(metric, real.(roots))[1]
    elseif numreals == 2 #case3
        if abs(imag(roots[4])) < T(1e-10)
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        return Ir_o_case3(metric, roots)
    else #case 4
        return Ir_o_case4(metric, roots)
    end
    return T(NaN)
end

function Ir_o_case1_and_2(::Kerr{T}, roots::NTuple{4}) where {T}
    roots = real.(pix.roots)
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    k = r32 * r41 / (r31 * r42)
    x2_o = √abs((rs - r4) / (rs - r3) * r31 / r41)
    coef = 2 / √real(r31 * r42)
    Ir_o = (x2_o > one(T)) ? T(NaN) : coef * JacobiElliptic.F(asin(x2_o), k)

    return Ir_o
end

function Ir_o_case3(::Kerr{T}, roots::NTuple{4}) where {T}
    roots = pix.roots
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    r1, r2, r21 = real.((r1, r2, r21))

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    temprat = B * (rs - r2) / (A * (rs - r1))
    x3_o = ((one(T) - temprat) / (one(T) + temprat))
    coef = one(T) * √inv(A * B)
    Ir_o = coef * JacobiElliptic.F((acos(x3_o)), k3)

    return Ir_o
end

function Ir_o_case4(::Kerr{T}, roots::NTuple{4}) where {T}
    roots = pix.roots
    _, r1, _, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return T(NaN)
    end
    C = √real(r31 * r42)
    D = √real(r32 * r41)
    k4 = 4C * D / (C + D)^2
    a2 = abs(imag(r1))
    b1 = real(r4)

    k4 = T(4) * C * D / (C + D)^2

    go = √max((T(4)a2^2 - (C - D)^2) / ((C + D)^2 - T(4)a2^2), zero(T))
    x4_o = (rs + b1) / a2
    coef = 2 / (C + D)
    Ir_o_coef = JacobiElliptic.F(atan(x4_o) + atan(go), k4)
    return coef*Ir_o_coef
end

"""
Returns the antiderivative \$I_ϕ=\\int\\frac{a(2Mr-a\\lambda)}{\\sqrt{\\Delta\\mathcal{R(r)}}}dr\$.
See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `pix`: BasicPixel
- `rs` : Emission radius
- `θs` : Emission inclination
- `τ` : Mino time
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
"""
function Iϕ(metric::Kerr{T}, pix::BasicPixel, rs, θo, τ, νr) where T
    λtemp = λ(metric, pix.location[1], θo)
    roots = pix.roots
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return Iϕ_case2(metric, pix, rs, τ, νr, λtemp)
    elseif numreals == 2 #case3
        return Iϕ_case3(metric, pix, rs, τ, λtemp)
    end #case4

    return Iϕ_case4(metric, pix, rs, τ, λtemp)
end

function Iϕ_case2(metric::Kerr{T}, pix::BasicPixel, rs, τ, νr, λ) where {T}
    roots = real.(pix.roots)
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

    k = r32 * r41 / (r31 * r42)
    x2_s = √((rs - r4) / (rs - r3) * r31 / r41)
    !(-1 < x2_s < 1) && return T(NaN)

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_s), k)
    Πm_s = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_s), k)

    Ip = - τ / rp3 - (-1)^νr*Πp_s
    Im = - τ / rm3 - (-1)^νr*Πm_s

    return pix.Iϕ + 2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
end

function Iϕ_case3(metric::Kerr{T}, pix::BasicPixel, rs, τ, λ) where {T}
    roots = pix.roots
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
        return T(NaN)
    end
    A, B = √A2, √B2

    k3 = ((A + B)^2 - r21^2) / (4 * A * B)

    temprat = B * (rs - r2) * real(_pow(A * (rs - r1), -one(T)))
    x3_s = clamp(((one(T) - temprat) * real(_pow(one(T) + temprat, -one(T)))), -one(T), one(T))
    φ_s = acos(x3_s)

    (isnan(φ_s)) && return T(NaN)

    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

    R1p_s = R1(αp, φ_s, k3)
    R1m_s = R1(αm, φ_s, k3)

    Ip = -inv(B * rp2 + A * rp1) * ((B + A) * τ + 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (- R1p_s))
    Im = -inv(B * rm2 + A * rm1) * ((B + A) * τ + 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (- R1m_s))

    return pix.Iϕ + (2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im))
end

function Iϕ_case4(metric::Kerr{T}, pix::BasicPixel, rs, τ, λ) where {T}
    roots = pix.roots
    r1, _, _, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    a = metric.spin
    a2 = a * a
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return T(NaN)
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
    S1m_s = S1(gm, atan(x4_s) + atan(go), k4)

    Ip = go / (a2 * (1 - go * x4_p)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * ( - S1p_s))
    Im = go / (a2 * (1 - go * x4_m)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * ( - S1m_s))

    return pix.Iϕ + 2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
end

"""
Returns the antiderivative \$I_ϕ=\\int\\frac{a(2Mr-a\\lambda)}{\\sqrt{\\Delta\\mathcal{R(r)}}}dr\$.
See [`r_potential(x)`](@ref) for an implementation of \$\\mathcal{R}(r)\$.

# Arguments

- `metric`: Kerr{T} metric
- `roots` : Radial roots
- `λ`  : Reduced azimuthal angular momentum
"""
function Iϕ_inf(metric::Kerr{T}, roots, λ) where T
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

    k = r32 * r41 / (r31 * r42)
    x2_o = √(r31 / r41)

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_o = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)
    Πm_o = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_o), k)

    Ipo_inf_m_I0_terms = -Πp_o
    Imo_inf_m_I0_terms = -Πm_o 

    return 2a / (rp - rm) * ((rp - a * λ / 2) * Ipo_inf_m_I0_terms - (rm - a * λ / 2) * Imo_inf_m_I0_terms)
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
        return T(NaN)
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

    Ip = -inv(B * rp2 + A * rp1) * ( 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * R1p_o)
    Im = -inv(B * rm2 + A * rm1) * ( 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * R1m_o)

    return (2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im))
end

function Iϕ_inf_case4(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    r1, _, _, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    a = metric.spin
    a2 = a * a
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return T(NaN)
    end

    C = √real(r31 * r42)
    D = √real(r32 * r41)
    k4 = T(4) * C * D / (C + D)^2
    a2 = abs(imag(r1))
    b1 = real(r4)

    k4 = 4 * C * D / (C + D)^2

    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    S1p_o = S1(gp, T(π / 2) + atan(go), k4)
    S1m_o = S1(gm, T(π / 2) + atan(go), k4)

    Ip = go / (a2 * (1 - go * x4_p)) * ( - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1p_o))
    Im = go / (a2 * (1 - go * x4_m)) * ( - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1m_o))

    return 2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
end

function It(metric::Kerr{T}, pix::BasicPixel, rs, θo, τ, νr) where {T}
    λtemp = λ(metric, pix.location[1], θo)
    roots = pix.roots
    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 2
        return It_case2(metric, pix, rs, τ, νr, λtemp)
    elseif numreals == 2 #case3
        return It_case3(metric, pix, rs, τ, λtemp)
    end #case4

    return It_case4(metric, pix, rs, τ, λtemp)
end

function It_case2(metric::Kerr{T}, pix::BasicPixel, rs, τ, νr, λ) where {T}
    roots = real.(pix.roots)
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
    !(-1 < x2_s < 1) && return T(NaN)

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
    I1_total = r3 * I0_total + r43 * ((νr ? -1 : 1) * Π1_s)# Removed the logarithmic divergence
    #equation B38
    I2_s = √(evalpoly(rs, poly_coefs)) / (rs - r3) - E_s
    I2_total = - (r1 * r4 + r2 * r3) / 2 * τ + (νr ? -1 : 1) * I2_s# Assymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_s), k)
    Πm_s = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_s), k)

    Ip_total = - τ / rp3 - (-1)^νr*Πp_s
    Im_total = - τ / rm3 - (-1)^νr*Πm_s 

    return pix.It -(4 / (rp - rm) * (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) + 4 * I0_total + 2 * I1_total + I2_total)
end

function It_case3(metric::Kerr{T}, pix::BasicPixel, rs, τ, λ) where {T}
    roots = pix.roots
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
        return T(NaN)
    end
    A, B = √A2, √B2

    k3 = real(((A + B)^2 - r21^2) / (4 * A * B))

    temprat = B * (rs - r2) * _pow(A * (rs - r1), -one(T))
    x3_s = clamp(real(((one(T) - temprat) * _pow(one(T) + temprat, -one(T)))), -one(T), one(T))
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
    I2_total = ((((B * r2 + A * r1) / (B + A))^2) * I0_total - 2 * (B * r2 + A * r1) / (B + A) * (-Π1_s) - √(A * B) * ( - Π2_s)) #+ (B * r2 + A * r1) / (A + B)
    Ip_total = -inv(B * rp2 + A * rp1) * ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (- R1(αp, φ_s, k3)))
    Im_total = -inv(B * rm2 + A * rm1) * ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (- R1(αm, φ_s, k3)))

    return pix.It -(4 / (rp - rm) * (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) + 4 * I0_total + 2 * I1_total + I2_total)
end

function It_case4(metric::Kerr{T}, pix::BasicPixel, rs, τ, λ) where {T}
    roots = pix.roots
    a = metric.spin
    r1, _, _, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return T(NaN)
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
    S1m_s = S1(gm, atan(x4_s) + atan(go), k4)

    Π1_s = 2 / (C + D) * (a2 / go * (1 + go^2)) * S1(go, atan(x4_s) + atan(go), k4)
    Π2_s = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * S2(go, atan(x4_s) + atan(go), k4)

    I0_total = τ
    # Removed logarithmic divergence
    I1_total = (a2 / go - b1) * τ - ( - Π1_s) 

    # Removed linear divergence
    I2_total =
        ((a2 / go - b1)^2) * τ - 2(a2 / go - b1) * (- Π1_s) + (- Π2_s) 
    Ip_total = go / (a2 * (1 - go * x4_p)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * ( - S1p_s))
    Im_total = go / (a2 * (1 - go * x4_m)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * ( - S1m_s))

    return pix.It -(4 / (rp - rm) * (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) + 4 * I0_total + 2 * I1_total + I2_total)# + (logdiv + lineardiv)
end

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

    k = r32 * r41 / (r31 * r42)
    x2_o = √abs(r31 / r41)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    E_o = √(r31 * r42) * JacobiElliptic.E(asin(x2_o), k)

    #equation B37
    I1_total = log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k))# Removed the logarithmic divergence
    #equation B38
    I2_total = r3 - E_o # Assymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_o = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)
    Πm_o = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_o), k)

    Ip_total = -Πp_o
    Im_total = -Πm_o

    return -(4 / (rp - rm) * (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) + 2 * I1_total + I2_total)
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
        return T(NaN)
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
    I1_total =  Π1_o + log(16 * r21^2 / ((A2 - B2)^2 + 4 * A * B * r21^2)) / 2 
    # Removed linear divergence
    I2_total = ( - √(A * B) * (Π2_o)) + (B * r2 + A * r1) / (A + B)
    Ip_total = -inv(B * rp2 + A * rp1) * ( 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (R1(αp, φ_o, k3) ))
    Im_total = -inv(B * rm2 + A * rm1) * ( 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (R1(αm, φ_o, k3) ))

    return -(4 / (rp - rm) * (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) + 2 * I1_total + I2_total)
end

function It_inf_case4(metric::Kerr{T}, roots::NTuple{4}, λ) where {T}
    a = metric.spin
    r1, _, _, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    if real(r32 * r41) < zero(T) || real(r31 * r42) < zero(T)
        return T(NaN)
    end

    C = √real(r31 * r42)
    D = √real(r32 * r41)
    k4 = T(4) * C * D / (C + D)^2
    a2 = abs(imag(r1))
    b1 = real(r4)

    k4 = 4 * C * D / (C + D)^2

    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    (isnan(x4_p) || isnan(x4_m)) && return T(NaN)
    S1p_o = S1(gp, T(π / 2) + atan(go), k4)
    S1m_o = S1(gm, T(π / 2) + atan(go), k4)

    Π1_o = 2 / (C + D) * (a2 / go * (1 + go^2)) * regularizedS1(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end
    Π2_o = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * regularizedS2(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end

    I1_total = - (Π1_o) + 1/2*log(
        (16*(1 + go^2 - sqrt((1 + go^2)*(1 + go^2 - k4)))*(1 + go^2 - 
   k4))/((C + D)^2*((1 + go^2)^2 - k4)*k4*(1 + sqrt(
   1 - k4/(1 + go^2)))))

    # Removed linear divergence
    I2_total = - 2(a2 / go - b1) * (Π1_o) + (Π2_o) -
        (
            (16 * a2^4 + (C^2 - D^2)^2 - 8 * (a2^2) * (C^2 + D^2) + 8 * (a2^3) * (C + D - 2 * b1 * go) + 
            2 * a2 * (C + D) * (-(C - D)^2 + 2 * b1 * (C + D) * go)) / (4 * a2 * (4 * a2^2 - (C + D)^2) * go)
        )
    Ip_total = go / (a2 * (1 - go * x4_p)) * ( - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1p_o))
    Im_total = go / (a2 * (1 - go * x4_m)) * ( - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (S1m_o))

    return -(4 / (rp - rm) * (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total)  + 2 * I1_total + I2_total)# + (logdiv + lineardiv)
end

"""
Returns the radial integrals for the case where there are four real roots in the radial potential, with roots outside the horizon.
"""
function radial_integrals_case2(metric::Kerr{T}, rs, pix::BasicPixel, τ, νr) where {T}
    roots = real.(pix.roots)
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
    !(-1 < x2_s < 1) && return T(NaN), T(NaN), T(NaN), T(NaN), T(NaN)

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
    I1_total = pix.I1o_m_I0_terms + r3 * I0_total + r43 * (-1)^νr * Π1_s# Removed the logarithmic divergence
    #equation B38
    I2_s = √abs(evalpoly(rs, poly_coefs)) / (rs - r3) - E_s
    I2_total = pix.I2o_m_I0_terms - (r1 * r4 + r2 * r3) / 2 * τ + (-1)^νr * I2_s# Assymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_s), k)
    Πm_s = coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_s), k)

    Ip_total = pix.Ipo_m_I0_terms - τ / rp3 - (-1)^νr*Πp_s
    Im_total = pix.Imo_m_I0_terms - τ / rm3 - (-1)^νr*Πm_s 

    return I0_total, I1_total, I2_total, Ip_total, Im_total

end

"""
Returns the radial integrals for the case where there are two real roots in the radial potential
"""
function radial_integrals_case3(metric::Kerr{T}, rs, pix::BasicPixel, τ) where {T}
    roots = pix.roots
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


    I0_total = τ
    # Removed logarithmic divergence
    I1_total = pix.I1o_m_I0_terms + (B * r2 + A * r1) / (B + A) * I0_total - Π1_s 
    # Removed linear divergence
    I2_total = pix.I2o_m_I0_terms + ((((B * r2 + A * r1) / (B + A))^2) * I0_total - 2 * (B * r2 + A * r1) / (B + A) * (-Π1_s) - √(A * B) * (- Π2_s)) 
    Ip_total = pix.Ipo_m_I0_terms + -inv(B * rp2 + A * rp1) * ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (- R1(αp, φ_s, k3)))
    Im_total = pix.Imo_m_I0_terms + -inv(B * rm2 + A * rm1) * ((B + A) * I0_total + 2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (- R1(αm, φ_s, k3)))

    return I0_total, I1_total, I2_total, Ip_total, Im_total
end
"""
Returns the radial integrals for the case where there are no real roots in the radial potential
"""
function radial_integrals_case4(metric::Kerr{T}, rs, pix::BasicPixel, τ) where {T}
    roots = pix.roots
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
    I1_total = pix.I1o_m_I0_terms + (a2 / go - b1) * τ - (- Π1_s) 

    # Removed linear divergence
    I2_total = pix.I2o_m_I0_terms +
        ((a2 / go - b1)^2) * τ - 2(a2 / go - b1) * ( - Π1_s) + ( - Π2_s)
    Ip_total = pix.Ipo_m_I0_terms + go / (a2 * (1 - go * x4_p)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (- S1p_s))
    Im_total = pix.Imo_m_I0_terms + go / (a2 * (1 - go * x4_m)) * (τ - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * (- S1m_s))

    return I0_total, I1_total, I2_total, Ip_total, Im_total
end

"""
Returns the radial integrals for the case where there are four real roots in the radial potential, with roots outside the horizon.
"""
function radial_inf_integrals_case2(metric::Kerr{T}, roots::NTuple{4}) where {T}
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

    k = r32 * r41 / (r31 * r42)
    x2_o = √abs(r31 / r41)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    E_o = √(r31 * r42) * JacobiElliptic.E(asin(x2_o), k)

    #equation B37
    I1o_m_I0_terms =  log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k) )
    #equation B38
    I2o_m_I0_terms = r3 - E_o# Assymptotic Divergent piece is not included

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)

    Ipo_m_I0_terms = -coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)
    Imo_m_I0_terms = -coef_m * JacobiElliptic.Pi(rm3 * r41 / (rm4 * r31), asin(x2_o), k)

    return I1o_m_I0_terms, I2o_m_I0_terms, Ipo_m_I0_terms, Imo_m_I0_terms

end

"""
Returns the radial integrals for the case where there are two real roots in the radial potential
"""
function radial_inf_integrals_case3(metric::Kerr{T}, roots::NTuple{4}) where {T}
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

    I1o_m_I0_terms =  Π1_o + log(16 * r21^2 / ((A2 - B2)^2 + 4 * A * B * r21^2)) / 2
    # Removed linear divergence
    I2o_m_I0_terms = - √(A * B) * Π2_o + (B * r2 + A * r1) / (A + B)
    Ipo_m_I0_terms = -inv(B * rp2 + A * rp1) * (2 * r21 * √(A * B) / (B * rp2 - A * rp1) * (R1(αp, φ_o, k3)))
    Imo_m_I0_terms = -inv(B * rm2 + A * rm1) * (2 * r21 * √(A * B) / (B * rm2 - A * rm1) * (R1(αm, φ_o, k3)))

    return I1o_m_I0_terms, I2o_m_I0_terms, Ipo_m_I0_terms, Imo_m_I0_terms
end
"""
Returns the radial integrals for the case where there are no real roots in the radial potential
"""
function radial_inf_integrals_case4(metric::Kerr{T}, roots::NTuple{4}) where {T}
    r1, _, _, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    a2 = metric.spin^2
    rp = one(T) + √(one(T) - a2)
    rm = one(T) - √(one(T) - a2)

    C = √real(r31 * r42)
    D = √real(r32 * r41)
    a2 = abs(imag(r1))
    b1 = real(r4)
    k4 = 4 * C * D / (C + D)^2

    x4_p = (rp + b1) / a2
    x4_m = (rm + b1) / a2

    go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    gp = (go * x4_p - one(T)) / (go + x4_p)
    gm = (go * x4_m - one(T)) / (go + x4_m)

    S1p_o = S1(gp, T(π / 2) + atan(go), k4)
    S1m_o = S1(gm, T(π / 2) + atan(go), k4)
    Π1_o = 2 / (C + D) * (a2 / go * (1 + go^2)) * regularizedS1(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end
    Π2_o = 2 / (C + D) * (a2 / go * (1 + go^2))^2 * regularizedS2(go, T(π / 2) + atan(go), k4) # Divergence is removed, will be added back in the end

    # Removed logarithmic divergence
    I1o_m_I0_terms =  -Π1_o + 1/2*log(
        (16*(1 + go^2 - sqrt((1 + go^2)*(1 + go^2 - k4)))*(1 + go^2 - 
        k4))/((C + D)^2*((1 + go^2)^2 - k4)*k4*(1 + sqrt(
        1 - k4/(1 + go^2))))
    )

    # Removed linear divergence
    I2o_m_I0_terms = - 2(a2 / go - b1) * Π1_o + Π2_o -
        ((16 * a2^4 + 
            (C^2 - D^2)^2 - 8 * (a2^2) * (C^2 + D^2) + 
                8 * (a2^3) * (C + D - 2 * b1 * go) + 
                    2 * a2 * (C + D) * (-(C - D)^2 + 2 * b1 * (C + D) * go)
        ) / (4 * a2 * (4 * a2^2 - (C + D)^2) * go))
    Ipo_m_I0_terms = go / (a2 * (1 - go * x4_p)) * ( - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * S1p_o)
    Imo_m_I0_terms = go / (a2 * (1 - go * x4_m)) * ( - 2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * S1m_o)

    return I1o_m_I0_terms, I2o_m_I0_terms, Ipo_m_I0_terms, Imo_m_I0_terms
end

function radial_inf_integrals(met::Kerr{T}, roots::NTuple{4}) where T
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

##----------------------------------------------------------------------------------------------------------------------
# Inclination functions
##----------------------------------------------------------------------------------------------------------------------
"""
Mino time of trajectory between two inclinations for a given screen coordinate

# Arguments 

- `pix` : Pixel information
- `θs` : Emission inclination
- `θo` : Observer inclination
- `isindir` : Is the path direct or indirect?
- `n` : nth image in orde of amount of minotime traversed
"""
function mino_time(metric::Kerr{T}, pix, θs, θo, isindir, n) where {T}
    return Gθ(metric, pix, θs, θo, isindir, n)[1]
end

"""
Returns the antiderivative \$G_\\theta=\\int\\frac{d\\theta}{\\sqrt{\\Theta(\\theta)}}\$.
See [`θ_potential(x)`](@ref) for an implementation of \$\\Theta(\theta)\$.

# Arguments 

- `metric`: Kerr{T} metric
- `pix` : Pixel information
- `θs` : Emission inclination
- `θo` : Observer inclination
- `isindir` : Is the path direct or indirect?
- `n` : nth image ordered by minotime
"""
function Gθ(metric::Kerr{T}, pix::BasicPixel, θs, θo, isindir, n) where {T}
    α, β = pix.location
    signβ = sign(β)
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)

    a = metric.spin
    a2 = a^2
    Go, Gs, Ghat, minotime, isvortical = pix.absGθ, T(NaN), pix.Ghat, T(NaN), ηtemp < zero(T)

    cosθs = cos(θs)
    cosθo = cos(θo)
    isincone = abs(cosθs) < abs(cosθo)
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
    end
    if ((((signβ < 0) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
    end

    Δθ = (one(T) - (ηtemp + λtemp^2) / a2) / 2
    Δθ2 = Δθ^2
    desc = √(Δθ2 + ηtemp / a2)
    up = Δθ + desc
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
            return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
        end
        tempfac = one(T) / √abs(um * a2)
        Go *= (-one(T))^((θs > T(π / 2)))
        Gs = (-one(T))^((θs > T(π / 2))) * tempfac * JacobiElliptic.F(asin(√args), k)
    else
        args = cosθs / √(up)
        argo = cosθo / √(up)
        k = m
        if !(-one(T) < args < one(T)) || !(-one(T) < argo < one(T))
            return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
        end
        tempfac = one(T) / √abs(um * a^2)
        Gs = tempfac * JacobiElliptic.F(asin(args), k)
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    minotime = (isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs : n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return minotime, Gs, Go, Ghat, isvortical
end

function _absGθo_Gθhat(metric::Kerr{T}, θo, η, λ) where {T}
    a = metric.spin
    a2 = a^2
    Go, Ghat, isvortical = T(NaN), T(NaN), η < zero(T)

    Δθ = (one(T) - (η + λ^2) / a2) / 2
    Δθ2 = Δθ^2
    desc = √(Δθ2 + η / a2)
    up = Δθ + desc
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
            return Go
        end
        tempfac = one(T) / √abs(um * a2)
        Go = tempfac * JacobiElliptic.F(asin(√argo), k)
        Ghat = 2 * tempfac * JacobiElliptic.K(k)
    else
        argo = cosθo / √(up)
        k = m
        if !(-one(T) < argo < one(T))
            return Go
        end
        tempfac = one(T) / √abs(um * a^2)
        Go = tempfac * JacobiElliptic.F(asin(argo), k)
        Ghat = 2 * tempfac * JacobiElliptic.K(k)
    end

    return Go, Ghat
end

function Gs(metric::Kerr{T}, pix, θo, τ) where {T}
    α, β = pix.location
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    signβ = sign(β)

    τ == T(NaN) && return T(NaN)
    a = metric.spin
    Gs, Go, Ghat, isvortical = T(NaN), T(NaN), T(NaN), ηtemp < zero(T)

    Δθ = T(0.5) * (one(T) - (ηtemp + λtemp^2) / a^2)
    up = Δθ + √(Δθ^2 + ηtemp / a^2)
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
        Go = pix.absGθ
        Ghat = pix.Ghat
        Δτtemp = (τ % Ghat + (θo > T(π / 2) ? -one(T) : one(T)) * signβ * Go)
        n = floor(τ / Ghat)
        Δτ = (θo > T(π / 2) ? -one(T) : one(T)) * abs(argmin(abs, [(-one(T))^n * signβ * (Ghat - Δτtemp), (-one(T))^n * signβ * Δτtemp]))
    else
        argo = cos(θo) / √(up)
        k = m
        tempfac = inv(√abs(um * a^2))
        Go = pix.absGθ
        Ghat = pix.Ghat
        Δτtemp = (τ % Ghat + signβ * Go)
        n = floor(τ / Ghat)
        Δτ = argmin(abs, [(-one(T))^n * signβ * (Ghat - Δτtemp), (-one(T))^n * signβ * Δτtemp])
    end

    return Δτ
end

function Gϕ(metric::Kerr{T}, pix::BasicPixel, θs, θo, isindir, n) where {T}
    α, β = pix.location
    signβ = sign(β)
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)

    a = metric.spin
    Go, Gs, Ghat, ans, isvortical = pix.absGϕ, zero(T), pix.Gϕhat, zero(T), ηtemp < zero(T)

    isincone = abs(cos(θs)) < abs(cos(θo))
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
    end
    if ((((signβ < zero(T)) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
    end

    Δθ = (1 - (ηtemp + λtemp^2) / a^2) / T(2)
    up = Δθ + √(Δθ^2 + ηtemp / a^2)
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
            return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
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
            return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
        end
        tempfac = inv(√abs(um * a^2))
        Gs = tempfac * JacobiElliptic.Pi(up, asin(args), k)
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    ans = (isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs : n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return ans, Gs, Go, Ghat, isvortical
end

function _absGϕo_Gϕhat(metric::Kerr{T}, θo, η, λ) where {T}

    a = metric.spin
    Go, Ghat, isvortical = T(NaN), T(NaN), η < zero(T)

    Δθ = (1 - (η + λ^2) / a^2) / T(2)
    up = Δθ + √(Δθ^2 + η / a^2)
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
            return T(NaN), T(NaN)
        end
        tempfac = inv((1 - um) * √abs(um * a^2))
        argn = (up - um) / (1 - um)
        Go = tempfac * JacobiElliptic.Pi(argn, asin(√argo), k)
        Ghat = 2tempfac * JacobiElliptic.Pi(argn, k)
    else
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-one(T) < argo < one(T))
            return T(NaN), T(NaN)
        end
        tempfac = inv(√abs(um * a^2))
        Go = tempfac * JacobiElliptic.Pi(up, asin(argo), k)
        Ghat = 2tempfac * JacobiElliptic.Pi(up, k)
    end

    return Go, Ghat
end

function Gt(metric::Kerr{T}, pix::BasicPixel, θs, θo, isindir, n) where {T}
    α, β = pix.location
    signβ = sign(β)

    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    a = metric.spin
    Go, Gs, Ghat, ans, isvortical = pix.absGt, T(NaN), pix.Gthat, T(NaN), ηtemp < zero(T)

    isincone = abs(cos(θs)) < abs(cos(θo))
    if isincone && (isindir != ((signβ > zero(T)) ⊻ (θo > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
    end
    if ((((signβ < zero(T)) ⊻ (θs > T(π / 2))) ⊻ (n % 2 == 1)) && !isincone && !isvortical) || (isvortical && ((θo >= T(π / 2)) ⊻ (θs > T(π / 2))))
        return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
    end

    Δθ = (1 - (ηtemp + λtemp^2) / a^2) / T(2)
    up = Δθ + √(Δθ^2 + ηtemp / a^2)
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
            return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
        end
        tempfac = √abs(um / a^2)
        Go *= ((θs > T(π / 2)) ? -1 : 1) 
        Gs = ((θs > T(π / 2)) ? -1 : 1) * tempfac * JacobiElliptic.E(asin(√args), k)

    else
        args = cos(θs) / √(up)
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-one(T) < args < one(T)) || !(-one(T) < argo < one(T))
            return T(NaN), T(NaN), T(NaN), T(NaN), isvortical
        end
        tempfac = -2 * up * inv(√abs(um * a^2))
        Gs = tempfac * (JacobiElliptic.E(asin(args), k) - JacobiElliptic.F(asin(args), k)) / (2k)
    end

    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    ans = (isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs : n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    return ans, Gs, Go, Ghat, isvortical
end

function _absGto_Gthat(metric::Kerr{T}, θo, ηtemp, λtemp) where {T}
    a = metric.spin
    Go, Ghat, isvortical = T(NaN), T(NaN), ηtemp < zero(T)

    Δθ = (1 - (ηtemp + λtemp^2) / a^2) / T(2)
    up = Δθ + √(Δθ^2 + ηtemp / a^2)
    um = Δθ - √(Δθ^2 + ηtemp / a^2)
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
            return T(NaN), T(NaN)
        end
        tempfac = √abs(um / a^2)
        Go =  tempfac * JacobiElliptic.E(asin(√argo), k)
        Ghat = 2tempfac * JacobiElliptic.E(k)

    else
        argo = cosθo / √(up)
        #k = abs(m)
        if !(-one(T) < argo < one(T))
            return T(NaN), T(NaN)
        end
        tempfac = -2 * up * inv(√abs(um * a^2))
        Go = tempfac * (JacobiElliptic.E(asin(argo), k) - JacobiElliptic.F(asin(argo), k)) / (2k)
        Ghat = tempfac * (JacobiElliptic.E(k) - JacobiElliptic.K(k)) / k
    end

    return Go, Ghat
end