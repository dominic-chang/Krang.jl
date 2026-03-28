module KrangReactantExt

using Krang
import Krang.JacobiElliptic as JacobiElliptic
using Reactant
import Reactant: TracedRNumber

_isreal2 = Krang._isreal2
Ir_inf = Krang.Ir_inf
horizon = Krang.horizon
Ir_s = Krang.Ir_s
Iϕ_inf_case2 = Krang.Iϕ_inf_case2
Iϕ_inf_case3 = Krang.Iϕ_inf_case3
Iϕ_inf_case4 = Krang.Iϕ_inf_case4
It_inf_case2 = Krang.It_inf_case2
It_inf_case3 = Krang.It_inf_case3
It_inf_case4 = Krang.It_inf_case4
radial_inf_integrals_case2 = Krang.radial_inf_integrals_case2
radial_inf_integrals_case3 = Krang.radial_inf_integrals_case3
radial_inf_integrals_case4 = Krang.radial_inf_integrals_case4
_get_root_diffs = Krang._get_root_diffs
R1 = Krang.R1
S1 = Krang.S1
regularized_Pi = Krang.regularized_Pi
regularized_R1 = Krang.regularized_R1
regularized_R2 = Krang.regularized_R2
regularizedS1 = Krang.regularizedS1
regularizedS2 = Krang.regularizedS2
ellF = JacobiElliptic.F
ellK = JacobiElliptic.K
ellE = JacobiElliptic.E
ellPi = JacobiElliptic.Pi

function _regularize_zero(x, ::Type{T}) where {T}
    return ifelse(x == zero(T), eps(T), x)
end

function _regularize_equal(x, y, ::Type{T}) where {T}
    return ifelse(x == y, x + eps(T), x)
end

function _argmax_real3(x1, x2, x3)
    ans = x1
    @trace if real(x1) >= real(x2)
        ans = Base.ifelse(real(x1) >= real(x3), x1, x3)
    else
        ans = Base.ifelse(real(x2) >= real(x3), x2, x3)
    end
    return ans 
end

function _reactant_get_radial_roots(metric::Krang.Kerr, η, λ)
    TT = typeof(metric.spin)
    a = metric.spin

    a2 = a * a
    A = a2 - η - λ * λ
    A2 = A + A
    B = TT(2) * (η + (λ - a)^2)
    C0 = -a2 * η

    P = -A * A / TT(12) - C0
    Q = -A / TT(3) * (A * A / TT(36) + zero(TT)im - C0) - B * B / TT(8)

    Δ3 = -TT(4) * P * P * P - TT(27) * Q * Q
    ωp = (-Q / TT(2) + sqrt(-Δ3 / TT(108)) + zero(TT)im)^(TT(1 / 3))

    C = (
        (-TT(1 / 2) + TT(√3 / 2)im) * ωp,
        (-TT(1 / 2) - TT(√3 / 2)im) * ωp,
        (one(TT) + zero(TT)im) * ωp,
    )
    v = (
        -P / (TT(3) * C[1]),
        -P / (TT(3) * C[2]),
        -P / (TT(3) * C[3]),
    )
    ξ0 = _argmax_real3(C[1] + v[1], C[2] + v[2], C[3] + v[3]) - A / TT(3)
    ξ02 = ξ0 + ξ0

    predet1 = A2 + ξ02
    predet2 = (√TT(2) * B) * inv(sqrt(ξ0))
    det1 = sqrt(-(predet1 - predet2))
    det2 = sqrt(-(predet1 + predet2))

    sqrtξ02 = sqrt(ξ02)

    r1 = (-sqrtξ02 - det1) / 2
    r2 = (-sqrtξ02 + det1) / 2
    r3 = (sqrtξ02 - det2) / 2
    r4 = (sqrtξ02 + det2) / 2

    roots = (r1, r2, r3, r4)
    numreals = sum(_isreal2, roots) 
    return Base.ifelse((numreals == 2) & (abs(imag(r4)) < sqrt(eps(TT))), (r1, r4, r2, r3), roots)
end

Krang.get_radial_roots(metric::Krang.Kerr, η::TracedRNumber, λ) =
    _reactant_get_radial_roots(metric, η, λ)

Krang.get_radial_roots(metric::Krang.Kerr, η, λ::TracedRNumber) =
    _reactant_get_radial_roots(metric, Reactant.promote_to(TracedRNumber{typeof(metric.spin)}, η), λ)

Krang.get_radial_roots(metric::Krang.Kerr, η::TracedRNumber, λ::TracedRNumber) =
    _reactant_get_radial_roots(metric, η, λ)

function Krang.Ir_inf(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber})
    numreals = sum(_isreal2, roots)
    func1 = Krang.Ir_inf_case1_and_2
    func2 = Krang.Ir_inf_case3
    func3 = Krang.Ir_inf_case4
    result = roots[1]
    Reactant.@trace if numreals == 4
        result = func1(metric, Base.real.(roots))
    elseif numreals == 2
        result = func2(metric, roots)
    else
        result = func3(metric, roots)
    end
    return result
end

function Krang.Ir_s(metric::Krang.Kerr, rs, roots::NTuple{4,<:TracedRNumber}, νr)
    func1 = Krang.Ir_s_case1_and_2
    func2 = Krang.Ir_s_case3
    func3 = Krang.Ir_s_case4
    numreals = sum(_isreal2, roots) 
    result = roots[1]
    Reactant.@trace if numreals == 4
        result = func1(metric, rs, Base.real.(roots), νr)
    elseif numreals == 2
        result = func2(metric, rs, roots)
    else
        result = func3(metric, rs, roots)
    end
    return result
end

function Krang.Ir_s_case1_and_2(metric::Krang.Kerr, rs, roots::NTuple{4,<:TracedRNumber}, νr)
    TT = typeof(metric.spin)
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = Krang._get_root_diffs(roots...)

    k = r32 * r41 / (r31 * r42)
    x2_s = √abs((rs - r4) / (rs - r3) * r31 / r41)
    coef = 2 / √real(r31 * r42)
    Ir_s = Base.ifelse(x2_s > one(TT), TT(Base.Inf), coef * JacobiElliptic.F(asin(x2_s), k))

    return -(-1)^νr * Ir_s
end

function Krang.Iϕ_inf(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
    numreals = sum(_isreal2, roots)
    result = λ
    Reactant.@trace if numreals == 4
        result = Iϕ_inf_case2(metric, real.(roots), λ)
    elseif numreals == 2
        result = Iϕ_inf_case3(metric, roots, λ)
    else
        result = Iϕ_inf_case4(metric, roots, λ)
    end
    return result
end

function Krang.Iϕ_inf_case2(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
    TT = typeof(metric.spin)
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, r43 = Krang._get_root_diffs(roots...)
    a = metric.spin
    a2 = a * a
    rp = one(TT) + √(one(TT) - a2)
    rm = one(TT) - √(one(TT) - a2)
    rp3 = _regularize_zero(rp - r3, TT)
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4

    k = r32 * r41 / (r31 * r42)
    x2_o = √(r31 / r41)

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_o = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)

    arg = _regularize_equal(rm3 * r41 / (rm4 * r31), k, TT)
    Πm_o = coef_m * JacobiElliptic.Pi(arg, asin(x2_o), k)

    Ipo_inf_m_I0_terms = -Πp_o
    Imo_inf_m_I0_terms = -Πm_o

    return 2a / (rp - rm) *
           ((rp - a * λ / 2) * Ipo_inf_m_I0_terms - (rm - a * λ / 2) * Imo_inf_m_I0_terms)
end

function Krang.Iϕ_inf_case3(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
    TT = typeof(metric.spin)
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r1, r2, r21 = real.((r1, r2, r21))

    a = metric.spin
    a2 = a * a
    rp = one(TT) + √(one(TT) - a2)
    rm = one(TT) - √(one(TT) - a2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    ans = A2
    @trace if (A2 < zero(TT)) | (B2 < zero(TT))
        ans = TT(Base.Inf)
    else
        A, B = √A2, √B2
        k3 = ((A + B)^2 - r21^2) / (4 * A * B)
        x3_o = clamp((A - B) / (A + B), -one(TT), one(TT))
        φ_o = acos(x3_o)

        if isnan(φ_o)
            ans = TT(Base.NaN)
        else
            αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
            αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

            R1p_o = R1(αp, φ_o, k3)
            R1m_o = R1(αm, φ_o, k3)

            Ip = -inv(B * rp2 + A * rp1) * (2 * r21 * √(A * B) / (B * rp2 - A * rp1) * R1p_o)
            Im = -inv(B * rm2 + A * rm1) * (2 * r21 * √(A * B) / (B * rm2 - A * rm1) * R1m_o)

            ans = 2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
        end
    end
    return ans 
end

function Krang.Iϕ_inf_case4(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
    TT = typeof(metric.spin)
    _, r2, _, r4 = roots
    a = metric.spin
    a2 = a * a
    rp = one(TT) + √(one(TT) - a2)
    rm = one(TT) - √(one(TT) - a2)

    a1 = abs(imag(r4))
    a2r = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1 - a2r)^2 + (b1 - b2)^2)
    D = sqrt((a1 + a2r)^2 + (b1 - b2)^2)
    k4 = 4C * D / (C + D)^2

    x4_p = (rp + b1) / a2r
    x4_m = (rm + b1) / a2r

    go = √max((4a2r^2 - (C - D)^2) / ((C + D)^2 - TT(4) * a2r^2), zero(TT))
    gp = (go * x4_p - one(TT)) / (go + x4_p)
    gm = (go * x4_m - one(TT)) / (go + x4_m)

    result = x4_p
    Reactant.@trace if isnan(x4_p) | isnan(x4_m)
        result = TT(Base.NaN)
    else
        S1p_o = S1(gp, TT(Base.pi / 2) + atan(go), k4)
        S1m_o = S1(gm, TT(Base.pi / 2) + atan(go), k4)

        Ip = go / (a2r * (1 - go * x4_p)) *
             (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * S1p_o)
        Im = go / (a2r * (1 - go * x4_m)) *
             (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * S1m_o)

        result = 2a / (rp - rm) * ((rp - a * λ / 2) * Ip - (rm - a * λ / 2) * Im)
    end
    return result
end

function Krang.It_inf(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
    numreals = sum(_isreal2, roots) 
    result = λ
    Reactant.@trace if numreals == 4
        result = It_inf_case2(metric, real.(roots), λ)
    elseif numreals == 2
        result = It_inf_case3(metric, roots, λ)
    else
        result = It_inf_case4(metric, roots, λ)
    end
    return result
end

function Krang.It_inf_case2(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
    TT = typeof(metric.spin)
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = Krang._get_root_diffs(roots...)
    r43 = r4 - r3
    a = metric.spin
    rp = one(TT) + √(one(TT) - a^2)
    rm = one(TT) - √(one(TT) - a^2)
    rp3 = _regularize_zero(rp - r3, TT)
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4

    k = r32 * r41 / (r31 * r42)
    x2_o = √abs(r31 / r41)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    E_o = √(r31 * r42) * JacobiElliptic.E(asin(x2_o), k)
    I1_total = log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k))
    I2_total = r3 - E_o

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_o = coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)
    arg = _regularize_equal(rm3 * r41 / (rm4 * r31), k, TT)
    Πm_o = coef_m * JacobiElliptic.Pi(arg, asin(x2_o), k)

    Ip_total = -Πp_o
    Im_total = -Πm_o

    return -(4 / (rp - rm) *
             (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
             2 * I1_total +
             I2_total)
end

function Krang.It_inf_case3(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
    TT = typeof(metric.spin)
    r1, r2, _, _ = roots
    r21, r31, r32, r41, r42, _ = Krang._get_root_diffs(roots...)
    r21 = real(r21)
    r2 = real(r2)
    r1 = real(r1)
    a = metric.spin
    rp = one(TT) + √(one(TT) - a^2)
    rm = one(TT) - √(one(TT) - a^2)
    rp1 = real(rp - r1)
    rp2 = real(rp - r2)
    rm1 = real(rm - r1)
    rm2 = real(rm - r2)

    A2 = real(r32 * r42)
    B2 = real(r31 * r41)
    result = A2
    Reactant.@trace if (A2 < zero(TT)) | (B2 < zero(TT))
        result = TT(Base.Inf)
    else
        A, B = √A2, √B2
        k3 = real(((A + B)^2 - r21^2) / (4 * A * B))
        x3_o = min((A - B) / (A + B), one(TT))
        φ_o = acos(x3_o)

        if isnan(φ_o)
            result = TT(Base.NaN)
        else
            αo = (B + A) / (B - A)
            αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
            αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)

            Π1_o = 2 * r21 * √real(A * B) / (B2 - A2) * regularized_R1(αo, φ_o, k3)
            Π2_o = ((2 * r21 * √(A * B) / (B2 - A2))^2) * regularized_R2(αo, φ_o, k3)

            I1_total = Π1_o + log(16 * r21^2 / ((A2 - B2)^2 + 4 * A * B * r21^2)) / 2
            I2_total = (-√(A * B) * Π2_o) + (B * r2 + A * r1) / (A + B)
            Ip_total = -inv(B * rp2 + A * rp1) *
                       (2 * r21 * √(A * B) / (B * rp2 - A * rp1) * R1(αp, φ_o, k3))
            Im_total = -inv(B * rm2 + A * rm1) *
                       (2 * r21 * √(A * B) / (B * rm2 - A * rm1) * R1(αm, φ_o, k3))

            result = -(4 / (rp - rm) *
                       (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
                       2 * I1_total +
                       I2_total)
        end
    end
    return result
end

function Krang.It_inf_case4(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
    TT = typeof(metric.spin)
    a = metric.spin
    _, r2, _, r4 = roots

    rp = one(TT) + √(one(TT) - a^2)
    rm = one(TT) - √(one(TT) - a^2)

    a1 = abs(imag(r4))
    a2r = abs(imag(r2))
    b1 = real(r4)
    b2 = real(r2)
    C = sqrt((a1 - a2r)^2 + (b1 - b2)^2)
    D = sqrt((a1 + a2r)^2 + (b1 - b2)^2)
    k4 = 4C * D / (C + D)^2

    x4_p = (rp + b1) / a2r
    x4_m = (rm + b1) / a2r

    go = √max((4a2r^2 - (C - D)^2) / ((C + D)^2 - TT(4) * a2r^2), zero(TT))
    gp = (go * x4_p - one(TT)) / (go + x4_p)
    gm = (go * x4_m - one(TT)) / (go + x4_m)

    result = x4_p
    Reactant.@trace if isnan(x4_p) | isnan(x4_m)
        result = TT(Base.NaN)
    else
        S1p_o = S1(gp, TT(Base.pi / 2) + atan(go), k4)
        S1m_o = S1(gm, TT(Base.pi / 2) + atan(go), k4)

        Π1_o = 2 / (C + D) * (a2r / go * (1 + go^2)) * regularizedS1(go, TT(Base.pi / 2) + atan(go), k4)
        Π2_o = 2 / (C + D) * (a2r / go * (1 + go^2))^2 * regularizedS2(go, TT(Base.pi / 2) + atan(go), k4)

        I1_total = -Π1_o + 1 / 2 * log(
            (16 * (1 + go^2 - sqrt((1 + go^2) * (1 + go^2 - k4))) * (1 + go^2 - k4)) /
            ((C + D)^2 * ((1 + go^2)^2 - k4) * k4 * (1 + sqrt(1 - k4 / (1 + go^2))))
        )

        I2_total = -2(a2r / go - b1) * Π1_o + Π2_o - (
            (16 * a2r^4 + (C^2 - D^2)^2 - 8 * (a2r^2) * (C^2 + D^2) +
             8 * (a2r^3) * (C + D - 2 * b1 * go) +
             2 * a2r * (C + D) * (-(C - D)^2 + 2 * b1 * (C + D) * go)) /
            (4 * a2r * (4 * a2r^2 - (C + D)^2) * go)
        )
        Ip_total = go / (a2r * (1 - go * x4_p)) *
                   (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * S1p_o)
        Im_total = go / (a2r * (1 - go * x4_m)) *
                   (-2 / (C + D) * ((1 + go^2) / (go * (go + x4_p))) * S1m_o)

        result = -(4 / (rp - rm) *
                   (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
                   2 * I1_total +
                   I2_total)
    end
    return result
end

function Krang.radial_inf_integrals(met::Krang.Kerr, roots::NTuple{4,<:TracedRNumber})
    result = roots
    func1 = radial_inf_integrals_case2
    func2 = radial_inf_integrals_case3
    func3 = radial_inf_integrals_case4
    numreals = sum(_isreal2, roots) 
    Reactant.@trace if numreals == 4
        result = func1(met, roots)
    elseif numreals == 2
        result = func2(met, roots)
    else
        result = func3(met, roots)
    end
    return result
end

function Krang.radial_inf_integrals_case2(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber})
    TT = typeof(metric.spin)
    roots = real.(roots)
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = Krang._get_root_diffs(roots...)
    r43 = r4 - r3
    a = metric.spin
    rp = one(TT) + √(one(TT) - a^2)
    rm = one(TT) - √(one(TT) - a^2)
    rp3 = _regularize_zero(rp - r3, TT)
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4

    k = r32 * r41 / (r31 * r42)
    x2_o = √abs(r31 / r41)

    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    E_o = √(r31 * r42) * JacobiElliptic.E(asin(x2_o), k)

    I1o_m_I0_terms = log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k))
    I2o_m_I0_terms = r3 - E_o

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)

    Ipo_m_I0_terms = -coef_p * JacobiElliptic.Pi(rp3 * r41 / (rp4 * r31) + eps(TT), asin(x2_o), k)
    arg = _regularize_equal(rm3 * r41 / (rm4 * r31), k, TT)
    Imo_m_I0_terms = -coef_m * JacobiElliptic.Pi(arg, asin(x2_o), k)

    return I1o_m_I0_terms, I2o_m_I0_terms, Ipo_m_I0_terms, Imo_m_I0_terms
end

function Krang.total_mino_time(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber})
    numreals = sum(_isreal2, roots) 
    I0_inf = Ir_inf(metric, roots)
    rh = horizon(metric)
    τf = roots[1]
    Reactant.@trace if numreals == 4
        τf =  2 * I0_inf
    else
        τf = I0_inf - Ir_s(metric, rh, roots, true)
    end
    return τf
end

@inline function Krang._absGθo_Gθhat(metric::Krang.Kerr, θo, η::TracedRNumber, λ)
    T = typeof(metric.spin)
    a = metric.spin
    a2 = a^2
    Go, Ghat = zero(T), zero(T)
    isvortical = η < zero(T)

    Δθ = (one(T) - (η + λ^2) / a2) / T(2)
    Δθ2 = Δθ^2
    desc = √max(Δθ2 + η / a2, zero(T))
    up = min(Δθ + desc, one(T) - eps(T))
    um = Δθ - desc
    m = up / um
    cosθo = cos(θo)

    Reactant.@trace if isvortical
        argo = (cosθo^2 - um) / (up - um)
        if (zero(T) < argo) & (argo < one(T))
            tempfac = inv(√abs(um * a2))
            Go = tempfac * ellF(asin(√argo), one(T) - m)
            Ghat = T(2) * tempfac * ellK(one(T) - m)
        end
    else
        argo = cosθo / √up
        if (-one(T) < argo) & (argo < one(T))
            tempfac = inv(√abs(um * a2))
            Go = tempfac * ellF(asin(argo), m)
            Ghat = T(2) * tempfac * ellK(m)
        end
    end

    return Go, Ghat
end

@inline function Krang._absGϕo_Gϕhat(metric::Krang.Kerr, θo, η::TracedRNumber, λ)
    T = typeof(metric.spin)
    a = metric.spin
    a2 = a^2
    Go, Ghat = zero(T), zero(T)
    isvortical = η < zero(T)

    Δθ = (one(T) - (η + λ^2) / a2) / T(2)
    up = min(Δθ + √(Δθ^2 + η / a2), one(T) - eps(T))
    um = Δθ - √(Δθ^2 + η / a2)
    m = up / um
    cosθo = cos(θo)

    Reactant.@trace if isvortical
        argo = (cosθo^2 - um) / (up - um)
        if (zero(T) < argo) & (argo < one(T))
            tempfac = inv((one(T) - um) * √abs(um * a2))
            argn = (up - um) / (one(T) - um)
            Go = tempfac * ellPi(argn, asin(√argo), one(T) - m)
            Ghat = T(2) * tempfac * ellPi(argn, one(T) - m)
        end
    else
        argo = cosθo / √up
        if (-one(T) < argo) & (argo < one(T))
            tempfac = inv(√abs(um * a2))
            Go = tempfac * ellPi(up, asin(argo), m)
            Ghat = T(2) * tempfac * ellPi(up, m)
        end
    end

    return Go, Ghat
end

@inline function Krang._absGto_Gthat(metric::Krang.Kerr, θo, η::TracedRNumber, λ)
    T = typeof(metric.spin)
    a = metric.spin
    a2 = a^2
    Go, Ghat = zero(T), zero(T)
    isvortical = η < zero(T)

    Δθ = (one(T) - (η + λ^2) / a2) / T(2)
    up = min(Δθ + √(Δθ^2 + η / a2), one(T) - eps(T))
    um = Δθ - √(Δθ^2 + η / a2)
    m = up / um
    cosθo = cos(θo)

    Reactant.@trace if isvortical
        argo = (cosθo^2 - um) / (up - um)
        if (zero(T) < argo) & (argo < one(T))
            tempfac = √abs(um / a2)
            Go = tempfac * ellE(asin(√argo), one(T) - m)
            Ghat = T(2) * tempfac * ellE(one(T) - m)
        end
    else
        argo = cosθo / √up
        if (-one(T) < argo) & (argo < one(T))
            tempfac = -T(2) * up * inv(√abs(um * a2))
            Go = tempfac * (ellE(asin(argo), m) - ellF(asin(argo), m)) / (T(2) * m)
            Ghat = tempfac * (ellE(m) - ellK(m)) / m
        end
    end

    return Go, Ghat
end

function _reactant_slow_light_intensity_pixel(met::Krang.Kerr, α, β, θo)
    TT = typeof(met.spin)
    tempη = Krang.η(met, α, β, θo)
    tempλ = Krang.λ(met, α, θo)
    roots = Krang.get_radial_roots(met, tempη, tempλ)
    r1, r2, r3, r4 = roots
    numreals = sum(_isreal2, roots) 

    roots = Base.ifelse((numreals == 2) & (abs(imag(r4)) < sqrt(eps(TT))), (r1, r4, r2, r3), roots)

    I1, I2, Ip, Im = Krang.radial_inf_integrals(met, roots)
    I0_inf = Krang.Ir_inf(met, roots)
    τ_total = Krang.total_mino_time(met, roots)
    Iϕ_inf_temp = Krang.Iϕ_inf(met, roots, tempλ)
    It_inf_temp = Krang.It_inf(met, roots, tempλ)
    absGθo_Gθhat = Krang._absGθo_Gθhat(met, θo, tempη, tempλ)
    absGϕo_Gϕhat = Krang._absGϕo_Gϕhat(met, θo, tempη, tempλ)
    absGto_Gthat = Krang._absGto_Gthat(met, θo, tempη, tempλ)

    return Krang.SlowLightIntensityPixel(
        met,
        (α, β),
        roots,
        I0_inf,
        τ_total,
        Iϕ_inf_temp,
        It_inf_temp,
        I1,
        I2,
        Ip,
        Im,
        absGθo_Gθhat,       
        absGϕo_Gϕhat,       
        absGto_Gthat,       
        θo,
        tempη,
        tempλ,
    )
end

Ts = (TracedRNumber, Any)
for (αT, βT, θT) in Iterators.product(Ts, Ts, Ts)
    if αT == βT == θT == Any
        continue
    end

    @eval function Krang.SlowLightIntensityPixel(
        met::Krang.Kerr,
        α::$αT,
        β::$βT,
        θo::$θT,
    ) 
        return _reactant_slow_light_intensity_pixel(met, α, β, θo)
    end
end

end
