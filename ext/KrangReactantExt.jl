module KrangReactantExt
#include("SlowLightIntensityPixel.jl")

using Krang
import Krang.JacobiElliptic as JacobiElliptic
using Reactant
import Reactant: TracedRNumber


types = (TracedRNumber, Any)
Ts = Union{TracedRNumber, Any}
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
radial_w_I0_terms_integrals_case2 = Krang.radial_w_I0_terms_integrals_case2 
radial_w_I0_terms_integrals_case3 = Krang.radial_w_I0_terms_integrals_case3
radial_w_I0_terms_integrals_case4 = Krang.radial_w_I0_terms_integrals_case4
_get_root_diffs = Krang._get_root_diffs
screen_coordinate = Krang.screen_coordinate
metric = Krang.metric
inclination = Krang.inclination
roots = Krang.roots
I0_inf = Krang.I0_inf
absGθo_Gθhat = Krang.absGθo_Gθhat
αboundary = Krang.αboundary
βboundary = Krang.βboundary
Intersection = Krang.Intersection
R1 = Krang.R1
R2 = Krang.R2
S1 = Krang.S1
S2 = Krang.S2
regularized_Pi = Krang.regularized_Pi
regularized_R1 = Krang.regularized_R1
regularized_R2 = Krang.regularized_R2
regularizedS1 = Krang.regularizedS1
regularizedS2 = Krang.regularizedS2
ellF = JacobiElliptic.F
ellK = JacobiElliptic.K
ellPi = JacobiElliptic.Pi
ellsn = JacobiElliptic.sn
ellcn = JacobiElliptic.cn
ellsc = JacobiElliptic.sc
_rs_case1_and_2 = Krang._rs_case1_and_2
_rs_case3 = Krang._rs_case3
_rs_case4 = Krang._rs_case4

@inline _ellE(m) = JacobiElliptic.ArithmeticGeometricMeanAlg.E(m)
@inline _ellE(φ, m) = JacobiElliptic.CarlsonAlg.E(φ, m)

function _regularize_zero(x, ::Type{T}) where {T}
    return Base.ifelse(x == zero(T), eps(T), x)
end

function _regularize_equal(x, y, ::Type{T}) where {T}
    return Base.ifelse(x == y, x + eps(T), x)
end

@inline function _argmax_real3(x1, x2, x3)
    return Base.ifelse(real(x1) >= real(x2),
        Base.ifelse(real(x1) >= real(x3), x1, x3),
        Base.ifelse(real(x2) >= real(x3), x2, x3)
    )
end

Reactant.@reactant_overlay function Krang.get_radial_roots(metric::Krang.Kerr, η, λ) 
    a = metric.spin
    TT = typeof(a)
    a2 = a * a
    A = a2 - η - λ * λ
    A2 = A + A
    B = TT(2) * (η + (λ - a)^2)
    C0 = -a2 * η

    P = -A * A / TT(12) - C0
    Q = -A / TT(3) * (A * A / TT(36) + zero(TT)im - C0) - B * B / TT(8)

    Δ3 = -TT(4) * P * P * P - TT(27) * Q * Q
    ωp = (-Q / TT(2) + sqrt(-Δ3 / TT(108)) + zero(TT)im)^(TT(1 / 3))

    
    C1 =(complex(-TT(1 / 2), TT(√3 / 2))) * ωp
    C2 = (complex(-TT(1 / 2), - TT(√3 / 2))) * ωp
    C3 = (complex(one(TT), zero(TT))) * ωp
    V1 = -P / (TT(3) * C1)
    V2 = -P / (TT(3) * C2)
    V3 = -P / (TT(3) * C3)
    ξ0 = _argmax_real3(C1 + V1, C2 + V2, C3 + V3) - A / TT(3)
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

    numreals = sum(_isreal2, (r1,r2,r3,r4)) 
    check = (numreals == 2) & (abs(imag(r4)) < sqrt(eps(TT)))
    return NTuple{4, typeof(r1)}((r1, Base.ifelse(check, r4, r2), Base.ifelse(check, r2, r3), Base.ifelse(check, r3, r4)))
end

#for (ηT, λT) in Iterators.product(types, types)
#    if ηT == λT == Any
#        continue
#    end
#
#    @eval function Krang.get_radial_roots(
#        met::Krang.Kerr,
#        η::$ηT,
#        λ::$λT,
#    ) 
#        return _reactant_get_radial_roots(met, η, λ)
#    end
#end

Reactant.@reactant_overlay function Krang.Ir_inf(metric::Krang.Kerr, roots)
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

Reactant.@reactant_overlay function Krang.Ir_s(metric::Krang.Kerr, rs, roots, νr)
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

Reactant.@reactant_overlay function Krang.Ir_s_case1_and_2(metric::Krang.Kerr, rs, roots, νr)
    TT = typeof(metric.spin)
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = Krang._get_root_diffs(roots...)

    k = r32 * r41 / (r31 * r42)
    x2_s = √abs((rs - r4) / (rs - r3) * r31 / r41)
    coef = 2 / √real(r31 * r42)
    Ir_s = Base.ifelse(x2_s > one(TT), TT(Base.Inf), coef * ellF(asin(x2_s), k))

    return -(-1)^νr * Ir_s
end

Reactant.@reactant_overlay function Krang.Iϕ_inf(metric::Krang.Kerr, roots, λ)
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

Reactant.@reactant_overlay function Krang.Iϕ_inf_case2(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
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
    Πp_o = coef_p * ellPi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)

    arg = _regularize_equal(rm3 * r41 / (rm4 * r31), k, TT)
    Πm_o = coef_m * ellPi(arg, asin(x2_o), k)

    Ipo_inf_m_I0_terms = -Πp_o
    Imo_inf_m_I0_terms = -Πm_o

    return 2a / (rp - rm) *
           ((rp - a * λ / 2) * Ipo_inf_m_I0_terms - (rm - a * λ / 2) * Imo_inf_m_I0_terms)
end

Reactant.@reactant_overlay function Krang.Iϕ_inf_case3(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
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

Reactant.@reactant_overlay function Krang.Iϕ_inf_case4(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
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

Reactant.@reactant_overlay function Krang.It_inf(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
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

Reactant.@reactant_overlay function Krang.It_inf_case2(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
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
    E_o = √(r31 * r42) * _ellE(asin(x2_o), k)
    I1_total = log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k))
    I2_total = r3 - E_o

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_o = coef_p * ellPi(rp3 * r41 / (rp4 * r31), asin(x2_o), k)
    arg = _regularize_equal(rm3 * r41 / (rm4 * r31), k, TT)
    Πm_o = coef_m * ellPi(arg, asin(x2_o), k)

    Ip_total = -Πp_o
    Im_total = -Πm_o

    return -(4 / (rp - rm) *
             (rp * (rp - a * λ / 2) * Ip_total - rm * (rm - a * λ / 2) * Im_total) +
             2 * I1_total +
             I2_total)
end

Reactant.@reactant_overlay function Krang.It_inf_case3(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
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

Reactant.@reactant_overlay function Krang.It_inf_case4(metric::Krang.Kerr, roots::NTuple{4,<:TracedRNumber}, λ)
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

Reactant.@reactant_overlay function Krang.radial_inf_integrals(met::Krang.Kerr, roots)
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

Reactant.@reactant_overlay function Krang.radial_inf_integrals_case2(metric::Krang.Kerr, roots)
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
    E_o = √(r31 * r42) * _ellE(asin(x2_o), k)

    I1o_m_I0_terms = log(16 / (r31 + r42)^2) / 2 + r43 * (coef * regularized_Pi(n, asin(inv(√n)), k))
    I2o_m_I0_terms = r3 - E_o

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)

    Ipo_m_I0_terms = -coef_p * ellPi(rp3 * r41 / (rp4 * r31) + eps(TT), asin(x2_o), k)
    arg = _regularize_equal(rm3 * r41 / (rm4 * r31), k, TT)
    Imo_m_I0_terms = -coef_m * ellPi(arg, asin(x2_o), k)

    return I1o_m_I0_terms, I2o_m_I0_terms, Ipo_m_I0_terms, Imo_m_I0_terms
end

Reactant.@reactant_overlay function Krang.total_mino_time(metric::Krang.Kerr, roots)
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

@inline Reactant.@reactant_overlay function Krang._absGθo_Gθhat(metric::Krang.Kerr, θo::A, η::B, λ::C) where {A,B,C}
    T = typeof(metric.spin)
    a = metric.spin
    a2 = a^2
    Go, Ghat = zero(T), zero(T)
    isvortical = η < zero(η)

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

#for (θT, ηT, λT) in Iterators.product(types, types, types)
#    if θT == ηT == λT == Any
#        continue
#    end
#
#    @eval function Krang._absGθo_Gθhat(
#        met::Krang.Kerr,
#        θ::$θT,
#        η::$ηT,
#        λ::$λT,
#    ) 
#        return _reactant_absGθo_Gθhat(met, θ, η, λ)
#    end
#end

Reactant.@reactant_overlay function Krang._absGϕo_Gϕhat(metric::Krang.Kerr, θo::A, η::B, λ::C) where {A,B,C}
    T = typeof(metric.spin)
    a = metric.spin
    a2 = a^2
    Go, Ghat = zero(T), zero(T)
    isvortical = η < zero(η)

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

#for (θT, ηT, λT) in Iterators.product(types, types, types)
#    if θT == ηT == λT == Any
#        continue
#    end
#
#    @eval function Krang._absGϕo_Gϕhat(
#        met::Krang.Kerr,
#        θ::$θT,
#        η::$ηT,
#        λ::$λT,
#    ) 
#        return _reactant_absGϕo_Gϕhat(met, θ, η, λ)
#    end
#end

Reactant.@reactant_overlay function Krang._absGto_Gthat(metric::Krang.Kerr, θo, η::TracedRNumber, λ)
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
            Go = tempfac * _ellE(asin(√argo), one(T) - m)
            Ghat = T(2) * tempfac * _ellE(one(T) - m)
        end
    else
        argo = cosθo / √up
        if (-one(T) < argo) & (argo < one(T))
            tempfac = -T(2) * up * inv(√abs(um * a2))
            Go = tempfac * (_ellE(asin(argo), m) - ellF(asin(argo), m)) / (T(2) * m)
            Ghat = tempfac * (_ellE(m) - ellK(m)) / m
        end
    end

    return Go, Ghat
end

@inline function _reactant_bad_angular_branch(signβ, θs, θo, isindir, n, isvortical, isincone, TT)
    halfpi = TT(Base.pi / 2)
    return (isincone & (isindir != ((signβ > zero(TT)) ⊻ (θo > halfpi)))) |
           (((((signβ < zero(TT)) ⊻ (θs > halfpi)) ⊻ (n % 2 == 1)) & !isincone & !isvortical) |
            (isvortical & ((θo >= halfpi) ⊻ (θs > halfpi))))
end

function _reactant_Gθ_vortical_minotime(Go, Ghat, θs, θo, signβ, n, isindir, um, up, m, a2, isincone, TT)
    args = (cos(θs)^2 - um) / (up - um)
    argo = (cos(θo)^2 - um) / (up - um)
    k = one(TT) - m
    invalid = !(((zero(TT) < argo) & (argo < one(TT))) & ((zero(TT) < args) & (args < one(TT))))
    safe_args = clamp(args, eps(TT), one(TT) - eps(TT))
    tempfac = inv(√abs(um * a2))
    signθs = Base.ifelse(θs > TT(Base.pi / 2), -one(TT), one(TT))
    Go2 = signθs * Go
    Gs = signθs * tempfac * ellF(asin(√safe_args), k)
    νθ = Base.ifelse(isincone, (n % 2 == 1) ⊻ (θo > θs), !isindir ⊻ (θs > TT(Base.pi / 2)))
    signs = Base.ifelse(νθ, -one(TT), one(TT))
    τcand = Base.ifelse(
        isindir,
        (n + 1) * Ghat - signβ * Go2 - signs * Gs,
        n * Ghat - signβ * Go2 - signs * Gs,
    )
    τ = Base.ifelse(invalid, zero(Go), τcand)
    validi = Base.ifelse(invalid, zero(Go), one(Go))

    return τ, validi
end

function _reactant_Gθ_nonvortical_minotime(Go, Ghat, θs, θo, signβ, n, isindir, um, up, m, a2, isincone, TT)
    args = cos(θs) / √up
    argo = cos(θo) / √up
    invalid = !(((-one(TT) < args) & (args < one(TT))) & ((-one(TT) < argo) & (argo < one(TT))))
    safe_args = clamp(args, -one(TT) + eps(TT), one(TT) - eps(TT))
    tempfac = inv(√abs(um * a2))
    Gs = tempfac * ellF(asin(safe_args), m)
    νθ = Base.ifelse(isincone, (n % 2 == 1) ⊻ (θo > θs), !isindir ⊻ (θs > TT(Base.pi / 2)))
    signs = Base.ifelse(νθ, -one(TT), one(TT))
    τcand = Base.ifelse(
        isindir,
        (n + 1) * Ghat - signβ * Go - signs * Gs,
        n * Ghat - signβ * Go - signs * Gs,
    )
    τ = Base.ifelse(invalid, zero(Go), τcand)
    validi = Base.ifelse(invalid, zero(Go), one(Go))

    return τ, validi
end

function _reactant_Gϕ_minotime(pix::Krang.SlowLightIntensityPixel, θs, isindir, n)
    α, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)
    signβ = sign(β)
    ηtemp = Krang.η(pix)
    λtemp = Krang.λ(pix)
    TT = typeof(met.spin)
    a = met.spin
    Go, Ghat = Krang.absGϕo_Gϕhat(pix)
    ans = zero(TT)
    isvortical = ηtemp < zero(TT)
    isincone = abs(cos(θs)) < abs(cos(θo))

    bad = _reactant_bad_angular_branch(signβ, θs, θo, isindir, n, isvortical, isincone, TT)

    Δθ = (one(TT) - (ηtemp + λtemp^2) / a^2) / TT(2)
    up = min(Δθ + √(Δθ^2 + ηtemp / a^2), one(TT) - eps(TT))
    um = Δθ - √(Δθ^2 + ηtemp / a^2)
    m = up / um
    k = Base.ifelse(isvortical, one(TT) - m, m)
    args_v = (cos(θs)^2 - um) / (up - um)
    argo_v = (cos(θo)^2 - um) / (up - um)
    args_n = cos(θs) / √up
    argo_n = cos(θo) / √up
    invalid_v = !(((zero(TT) < argo_v) & (argo_v < one(TT))) & ((zero(TT) < args_v) & (args_v < one(TT))))
    invalid_n = !(((-one(TT) < args_n) & (args_n < one(TT))) & ((-one(TT) < argo_n) & (argo_n < one(TT))))

    safe_args_v = clamp(args_v, eps(TT), one(TT) - eps(TT))
    safe_args_n = clamp(args_n, -one(TT) + eps(TT), one(TT) - eps(TT))

    argn = (up - um) / (one(TT) - um)
    tempfac_v = inv((one(TT) - um) * √abs(um * a^2))
    tempfac_n = inv(√abs(um * a^2))
    signθs = Base.ifelse(θs > TT(Base.pi / 2), -one(TT), one(TT))
    Go_v = signθs * Go
    Gs_v = signθs * tempfac_v * ellPi(argn, asin(√safe_args_v), k)
    Gs_n = tempfac_n * ellPi(up, asin(safe_args_n), k)

    νθ = Base.ifelse(isincone, (n % 2 == 1) ⊻ (θo > θs), !isindir ⊻ (θs > TT(Base.pi / 2)))
    signs = Base.ifelse(νθ, one(TT), -one(TT))
    anscand_v = Base.ifelse(
        isindir,
        (n + 1) * Ghat - signβ * Go_v + signs * Gs_v,
        n * Ghat - signβ * Go_v + signs * Gs_v,
    )
    anscand_n = Base.ifelse(
        isindir,
        (n + 1) * Ghat - signβ * Go + signs * Gs_n,
        n * Ghat - signβ * Go + signs * Gs_n,
    )
    validi_v = Base.ifelse(invalid_v, zero(TT), one(TT))
    validi_n = Base.ifelse(invalid_n, zero(TT), one(TT))
    anscand = Base.ifelse(isvortical, anscand_v, anscand_n)
    validi = Base.ifelse(isvortical, validi_v, validi_n)
    ans = Base.ifelse(bad, zero(TT), anscand)
    validi = Base.ifelse(bad, zero(TT), validi)

    return ans, validi > zero(TT)
end

function _reactant_Gt_minotime(pix::Krang.SlowLightIntensityPixel, θs, isindir, n)
    α, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)
    signβ = sign(β)
    ηtemp = Krang.η(pix)
    λtemp = Krang.λ(pix)
    TT = typeof(met.spin)
    a = met.spin
    Go, Ghat = Krang.absGto_Gthat(pix)
    ans = zero(TT)
    isvortical = ηtemp < zero(TT)
    isincone = abs(cos(θs)) < abs(cos(θo))

    bad = _reactant_bad_angular_branch(signβ, θs, θo, isindir, n, isvortical, isincone, TT)

    Δθ = (one(TT) - (ηtemp + λtemp^2) / a^2) / TT(2)
    up = min(Δθ + √(Δθ^2 + ηtemp / a^2), one(TT) - eps(TT))
    um = Δθ - √(Δθ^2 + ηtemp / a^2)
    m = up / um
    k = Base.ifelse(isvortical, one(TT) - m, m)
    args_v = (cos(θs)^2 - um) / (up - um)
    argo_v = (cos(θo)^2 - um) / (up - um)
    args_n = cos(θs) / √up
    argo_n = cos(θo) / √up
    invalid_v = !(((zero(TT) < argo_v) & (argo_v < one(TT))) & ((zero(TT) < args_v) & (args_v < one(TT))))
    invalid_n = !(((-one(TT) < args_n) & (args_n < one(TT))) & ((-one(TT) < argo_n) & (argo_n < one(TT))))

    safe_args_v = clamp(args_v, eps(TT), one(TT) - eps(TT))
    safe_args_n = clamp(args_n, -one(TT) + eps(TT), one(TT) - eps(TT))

    tempfac_v = √abs(um / a^2)
    tempfac_n = -TT(2) * up * inv(√abs(um * a^2))
    signθs = Base.ifelse(θs > TT(Base.pi / 2), -one(TT), one(TT))
    Go_v = signθs * Go
    Gs_v = signθs * tempfac_v * _ellE(asin(√safe_args_v), k)
    Gs_n = tempfac_n * (_ellE(asin(safe_args_n), k) - ellF(asin(safe_args_n), k)) / (TT(2) * k)

    νθ = Base.ifelse(isincone, (n % 2 == 1) ⊻ (θo > θs), !isindir ⊻ (θs > TT(Base.pi / 2)))
    signs = Base.ifelse(νθ, one(TT), -one(TT))
    anscand_v = Base.ifelse(
        isindir,
        (n + 1) * Ghat - signβ * Go_v + signs * Gs_v,
        n * Ghat - signβ * Go_v + signs * Gs_v,
    )
    anscand_n = Base.ifelse(
        isindir,
        (n + 1) * Ghat - signβ * Go + signs * Gs_n,
        n * Ghat - signβ * Go + signs * Gs_n,
    )
    validi_v = Base.ifelse(invalid_v, zero(TT), one(TT))
    validi_n = Base.ifelse(invalid_n, zero(TT), one(TT))
    anscand = Base.ifelse(isvortical, anscand_v, anscand_n)
    validi = Base.ifelse(isvortical, validi_v, validi_n)
    ans = Base.ifelse(bad, zero(TT), anscand)
    validi = Base.ifelse(bad, zero(TT), validi)

    return ans, validi > zero(TT)
end

@inline function _reactant_isapprox(x, y, ::Type{T}) where {T}
    rtol = sqrt(eps(T))
    return abs(x - y) <= rtol * max(abs(x), abs(y))
end

Reactant.@reactant_overlay function Krang.radial_w_I0_terms_integrals(met::Krang.Kerr, rs, roots, τ, νr)
    numreals = sum(_isreal2, roots)
    result = (zero(real(roots[1])), zero(real(roots[1])), zero(real(roots[1])), zero(real(roots[1])))
    Reactant.@trace if numreals == 4
        result = radial_w_I0_terms_integrals_case2(met, rs, real.(roots), τ, νr)
    elseif numreals == 2
        result = radial_w_I0_terms_integrals_case3(met, rs, roots, τ)
    else
        result = radial_w_I0_terms_integrals_case4(met, rs, roots, τ)
    end
    return result
end

Reactant.@reactant_overlay function Krang.radial_w_I0_terms_integrals_case2(metric::Krang.Kerr, rs, roots, τ, νr)
    T = typeof(metric.spin)
    r1, r2, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)
    r43 = r4 - r3
    a = metric.spin
    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)
    rp3 = _regularize_zero(rp - r3, T)
    rp4 = rp - r4
    rm3 = rm - r3
    rm4 = rm - r4

    k = r32 * r41 / (r31 * r42)
    x2_s = √abs((rs - r4) / (rs - r3) * r31 / r41)
    coef = 2 / √(r31 * r42)
    n = abs(r41 / r31)
    E_s = √(r31 * r42) * _ellE(asin(clamp(x2_s, zero(T), one(T))), k)
    Π1_s = coef * ellPi(n, asin(clamp(x2_s, zero(T), one(T))), k)
    I0_total = τ

    poly_coefs = (
        r1 * r2 * r3 * r4,
        -r2 * r3 * r4 + r1 * (r2 * (-r3 - r4) - r3 * r4),
        r3 * r4 + r2 * (r3 + r4) + r1 * (r2 + r3 + r4),
        -r1 - r2 - r3 - r4,
        one(T),
    )
    signνr = Base.ifelse(νr, -one(T), one(T))

    I1_total = -r3 * I0_total - r43 * signνr * Π1_s
    I2_s = √abs(evalpoly(rs, poly_coefs)) / (rs - r3) - E_s
    I2_total = (r1 * r4 + r2 * r3) / 2 * τ - signνr * I2_s

    coef_p = 2 / √(r31 * r42) * r43 / (rp3 * rp4)
    coef_m = 2 / √(r31 * r42) * r43 / (rm3 * rm4)
    Πp_s = coef_p * ellPi(rp3 * r41 / (rp4 * r31), asin(clamp(x2_s, zero(T), one(T))), k)
    arg = _regularize_equal(rm3 * r41 / (rm4 * r31), k, T)
    Πm_s = coef_m * ellPi(arg, asin(clamp(x2_s, zero(T), one(T))), k)

    Ip_total = τ / rp3 + signνr * Πp_s
    Im_total = τ / rm3 + signνr * Πm_s

    invalid = !(x2_s < one(T))
    z = zero(T)
    return (
        Base.ifelse(invalid, z, I1_total),
        Base.ifelse(invalid, z, I2_total),
        Base.ifelse(invalid, z, Ip_total),
        Base.ifelse(invalid, z, Im_total),
    )
end

Reactant.@reactant_overlay function Krang.radial_w_I0_terms_integrals_case3(metric::Krang.Kerr, rs, roots, τ)
    T = typeof(metric.spin)
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
    A = √A2
    B = √B2
    k3 = ((A + B)^2 - r21^2) / (4 * A * B)
    temprat = B * (rs - r2) * inv(A * (rs - r1))
    x3_s = real((one(T) - temprat) * inv(one(T) + temprat))
    φ_s = acos(clamp(x3_s, -one(T), one(T)))
    αo = (B + A) / (B - A)
    αp = (B * rp2 + A * rp1) / (B * rp2 - A * rp1)
    αm = (B * rm2 + A * rm1) / (B * rm2 - A * rm1)
    coef = 2 * r21 * √(A * B) / (B2 - A2)

    Π1_s = coef * R1(αo, φ_s, k3)
    Π2_s = (coef^2) * R2(αo, φ_s, k3)

    I0_total = τ
    I1_total = -(B * r2 + A * r1) / (B + A) * I0_total + Π1_s
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

    invalid = abs(x3_s) > one(T)
    z = zero(T)
    return (
        Base.ifelse(invalid, z, I1_total),
        Base.ifelse(invalid, z, I2_total),
        Base.ifelse(invalid, z, Ip_total),
        Base.ifelse(invalid, z, Im_total),
    )
end

Reactant.@reactant_overlay function Krang.radial_w_I0_terms_integrals_case4(metric::Krang.Kerr, rs, roots, τ)
    T = typeof(metric.spin)
    _, r2, _, r4 = roots
    a = metric.spin
    a2m = a * a
    rp = one(T) + √(one(T) - a2m)
    rm = one(T) - √(one(T) - a2m)

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

    S1p_s = S1(gp, atan(x4_s) + atan(go), k4)
    S1m_s = S1(gm, atan(x4_s) + atan(go), k4)
    Π1_s = 2 / (C + D) * (a2 / go * (one(T) + go^2)) * S1(go, atan(x4_s) + atan(go), k4)
    Π2_s = 2 / (C + D) * (a2 / go * (one(T) + go^2))^2 * S2(go, atan(x4_s) + atan(go), k4)

    I1_total = -(a2 / go - b1) * τ + (-Π1_s)
    I2_total = -((a2 / go - b1)^2) * τ + 2(a2 / go - b1) * (-Π1_s) - (-Π2_s)
    Ip_total =
        -go / (a2 * (one(T) - go * x4_p)) *
        (τ - 2 / (C + D) * ((one(T) + go^2) / (go * (go + x4_p))) * (-S1p_s))
    Im_total =
        -go / (a2 * (one(T) - go * x4_m)) *
        (τ - 2 / (C + D) * ((one(T) + go^2) / (go * (go + x4_p))) * (-S1m_s))

    invalid = isnan(x4_s) | isnan(x4_p) | isnan(x4_m)
    z = zero(T)
    return (
        Base.ifelse(invalid, z, I1_total),
        Base.ifelse(invalid, z, I2_total),
        Base.ifelse(invalid, z, Ip_total),
        Base.ifelse(invalid, z, Im_total),
    )
end

Reactant.@reactant_overlay function Krang.radial_integrals(pix::Krang.AbstractPixel, rs, τ, νr)
    met = metric(pix)
    I1_o, I2_o, Ip_o, Im_o = Krang.radial_inf_integrals(met, roots(pix))
    I1_s, I2_s, Ip_s, Im_s = Krang.radial_w_I0_terms_integrals(met, rs, roots(pix), τ, νr)
    return τ, I1_o - I1_s, I2_o - I2_s, Ip_o - Ip_s, Im_o - Im_s
end

function _reactant_θs(metric::Krang.Kerr, signβ, θo, η, λ, τ)
    T = typeof(metric.spin)
    a = metric.spin
    Ghat_2 = zero(T)
    isvortical = η < zero(T)

    Δθ = (one(T) - (η + λ^2) / a^2) / T(2)
    dsc = √(Δθ^2 + η / a^2)
    up = Δθ + dsc
    um = Δθ - dsc
    m = up / um

    ans = zero(T)
    argo = zero(T)
    τs = zero(T)
    τo = zero(T)
    τhat = zero(T)
    n = zero(Int)
    isindir = false
    tempfac = inv(√abs(um * a^2))

    Reactant.@trace if isvortical
        k = one(T) - m
        argo = clamp((cos(θo)^2 - um) / (up - um), zero(T), one(T))
        τo = tempfac * ellF(asin(√argo), k)
        Ghat_2 = tempfac * ellK(k)
        τhat = 2 * Ghat_2
        Δτtemp = mod(τ, τhat) + Base.ifelse(θo > T(Base.pi / 2), -one(T), one(T)) * signβ * τo
        n = unsafe_trunc(Int, τ / τhat)
        cand1 = τhat - Δτtemp
        cand2 = Δτtemp
        absτs = abs(Base.ifelse(abs(cand1) <= abs(cand2), cand1, cand2))
        τs = Base.ifelse(θo > T(Base.pi / 2), -one(T), one(T)) * absτs
        argr = ellsn(absτs / tempfac, k)^2
        ans = acos(Base.ifelse(θo > T(Base.pi / 2), -one(T), one(T)) * √((up - um) * argr + um))
        signo = Base.ifelse(ans > T(Base.pi / 2), -one(T), one(T))
        τo = signo * τo
    else
        k = m
        argo = clamp(cos(θo) / √up, -one(T), one(T))
        τo = tempfac * ellF(asin(argo), k)
        Ghat_2 = tempfac * ellK(k)
        τhat = Ghat_2 + Ghat_2
        Δτtemp = mod(τ, τhat) + signβ * τo
        n = unsafe_trunc(Int, τ / τhat)
        signo = Base.ifelse(n % 2 == 1, -one(T), one(T))
        cand1 = signo * signβ * (τhat - Δτtemp)
        cand2 = signo * signβ * Δτtemp
        τs = Base.ifelse(abs(cand1) <= abs(cand2), cand1, cand2)
        newargs = ellsn(τs / tempfac, k)
        ans = acos(√up * newargs)
    end

    isincone = abs(cos(ans)) < abs(cos(θo))
    τ1 = zero(T)
    Reactant.@trace if isincone
        νθ = (n % 2 == 1) ⊻ (θo > ans)
        τ1 = (n + 1) * τhat - signβ * τo + Base.ifelse(νθ, one(T), -one(T)) * τs
        isindir = _reactant_isapprox(τ1, τ, T)
    elseif ans < T(Base.pi / 2)
        τ1 = (n + 1) * τhat - signβ * τo - τs
        isindir = _reactant_isapprox(τ1, τ, T)
    else
        τ1 = (n + 1) * τhat - signβ * τo + τs
        isindir = _reactant_isapprox(τ1, τ, T)
    end

    return ans, τs, τo, τhat, n, isindir
end

Reactant.@reactant_overlay function Krang.emission_inclination(pix::Krang.AbstractPixel, τ)
    met = metric(pix)
    θo = inclination(pix)
    _, β = screen_coordinate(pix)
    return _reactant_θs(met, sign(β), θo, Krang.η(pix), Krang.λ(pix), τ)
end

function _reactant_Gθ_minotime(
    pix::Krang.SlowLightIntensityPixel,
    θs,
    isindir,
    n,
)
    _, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)
    signβ = sign(β)
    ηtemp = Krang.η(pix)
    λtemp = Krang.λ(pix)
    TT = typeof(met.spin)

    a = met.spin
    a2 = a^2
    Go, Ghat = absGθo_Gθhat(pix)
    isvortical = ηtemp < zero(TT)
    isincone = abs(cos(θs)) < abs(cos(θo))
    bad = _reactant_bad_angular_branch(signβ, θs, θo, isindir, n, isvortical, isincone, TT)

    τ = zero(Go)
    valid = true
    Δθ = (one(TT) - (ηtemp + λtemp^2) / a2) / TT(2)
    desc = √(Δθ^2 + ηtemp / a2)
    up = min(Δθ + desc, one(TT) - eps(TT))
    um = Δθ - desc
    m = up / um
    τvort, validvort = _reactant_Gθ_vortical_minotime(Go, Ghat, θs, θo, signβ, n, isindir, um, up, m, a2, isincone, TT)
    τnon, validnon = _reactant_Gθ_nonvortical_minotime(Go, Ghat, θs, θo, signβ, n, isindir, um, up, m, a2, isincone, TT)
    τcandidate = Base.ifelse(isvortical, τvort, τnon)
    validi = Base.ifelse(isvortical, validvort, validnon)
    τ = Base.ifelse(bad, zero(Go), τcandidate)
    validi = Base.ifelse(bad, zero(Go), validi)
    valid = validi > zero(validi)

    return τ, valid
end

Reactant.@reactant_overlay function Krang._rs_case1_and_2(pix::Krang.AbstractPixel, rh, τ)
    radial_roots = map(real, roots(pix))
    _, _, r3, r4 = radial_roots
    T = typeof(real(r3))
    _, r31, r32, r41, r42, _ = _get_root_diffs(radial_roots...)

    k = r32 * r41 / (r31 * r42)
    rh += eps(T)
    x2_s = √abs((rh - r4) / (rh - r3) * r31 / r41)
    coef = 2 / √real(r31 * r42)
    Ir_s = zero(T) 
    @trace if (x2_s < 1)
        Ir_s = coef * ellF(asin(x2_s), k)
    end

    fo = I0_inf(pix)
    X2 = √(r31 * r42) * (fo - τ) / 2

    return_vals = (T(Inf), false, false)
    @trace if !(((r4 < rh) & (τ > (fo - Ir_s))) | (τ > 2fo))# invalid case2
        sn = r41 * ellsn(X2, k)^2
        return_vals = (r31 * r4 - r3 * sn) / (r31 - sn), X2 > zero(T), true
    end
    return return_vals
end


Reactant.@reactant_overlay function Krang._rs_case3(pix::Krang.AbstractPixel, rh, τ)
    radial_roots = roots(pix)
    r1, r2, _, _ = radial_roots
    T = typeof(real(r2))
    r21, r31, r32, r41, r42, _ = _get_root_diffs(radial_roots...)
    r1, r2, r21 = real.((r1, r2, r21))

    fo = I0_inf(pix)
    A = √abs(r32 * r42)
    B = √abs(r31 * r41)
    k = (((A + B)^2 - r21^2) / (4 * A * B))
    temprat = B * (rh - r2) / (A * (rh - r1))
    x3_s = clamp(((1 - temprat) / (1 + temprat)), -1, 1)
    coef = 1 * √inv(A * B)
    Ir_s = coef * ellF((acos(x3_s)), k)
    X3 = √(A * B) * real(fo - τ)

    return_vals = (T(Inf), false, false)
    @trace if !((X3 < zero(T)) | (τ > (fo - Ir_s)))
        cn = ellcn(X3, k)
        num = -A * r1 + B * r2 + (A * r1 + B * r2) * cn
        den = -A + B + (A + B) * cn
        return_vals = real(num / den), X3 > zero(T), true
    end

    return return_vals
end

Reactant.@reactant_overlay function Krang._rs_case4(pix::Krang.AbstractPixel, rh, τ)
    radial_roots = roots(pix)
    _, r2, _, r4 = radial_roots
    T = typeof(real(r2))
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
    Ir_s = coef * ellF(atan(x4_s) + atan(go), k4)
    fo = I0_inf(pix)

    return_vals = (T(Inf), false, false)
    @trace if !(τ > (fo - Ir_s)) 
        X4 = (C + D) / T(2) * (fo - τ)
        num = go - ellsc(X4, k4)
        den = 1 + go * ellsc(X4, k4)
        return_vals = -(a2 * num / den + b1), X4 > zero(T), true
    end

    return return_vals
end

function _reactant_emission_radius_tau(
    pix::Krang.SlowLightIntensityPixel,
    τ,
) 
    rh = horizon(metric(pix))
    numreals = sum(_isreal2, roots(pix))
    rs4, valid4 = Krang.rs_case1_and_2(pix, rh, τ)
    rs2, valid2 = Krang._rs_case3(pix, rh, τ)
    rs0, valid0 = Krang._rs_case4(pix, rh, τ)
    is4 = numreals == 4
    is2 = numreals == 2
    rs = Base.ifelse(is4, rs4, Base.ifelse(is2, rs2, rs0))
    validi = Base.ifelse(is4, valid4, Base.ifelse(is2, valid2, valid0))
    valid = validi > zero(validi)

    return rs, valid
end

function _reactant_emission_radius_theta(
    pix::Krang.SlowLightIntensityPixel,
    θs,
    isindir,
    n,
) 
    α, β = screen_coordinate(pix)
    θo = inclination(pix)
    met = metric(pix)
    TT = typeof(met.spin)
    piT = TT(Base.pi)
    isincone = ((θo <= θs) & (θs <= (piT - θo))) | (((piT - θo) <= θs) & (θs <= θo))

    rs = zero(α)
    αmin = αboundary(met, θs)
    αoutside = abs(α) >= (αmin + eps(TT))
    βbound = Base.ifelse((!isincone) & αoutside, βboundary(met, α, θo, θs), zero(α))
    angular_gate = !((!isincone) & ((abs(β) + eps(TT)) < βbound))
    result = (zero(α), zero(angular_gate))
    Reactant.@trace result = if angular_gate
        τ, τvalid = _reactant_Gθ_minotime(pix, θs, isindir, n)
        if τvalid
            _reactant_emission_radius_tau(pix, τ)
        else
            (zero(α), zero(τvalid))
        end
    else
        (zero(α), zero(angular_gate))
    end
    rs, valid = result

    return rs, valid
end

function _emission_radius(
    pix::Krang.SlowLightIntensityPixel,
    τ,
) 
    TT = typeof(metric(pix).spin)
    numreals = sum(_isreal2, roots(pix))
    rh = horizon(metric(pix))
    fo = I0_inf(pix)

    rs4, valid4 = _reactant_rs_case1_and_2(pix, rh, τ, TT)
    rs2, valid2 = _reactant_rs_case3(pix, rh, τ, TT)
    rs0, valid0 = _reactant_rs_case4(pix, rh, τ, TT)

    radial_roots_real = real.(roots(pix))
    _, _, r3, r4 = radial_roots_real
    _, r31, _, r41, r42, _ = _get_root_diffs(radial_roots_real...)
    X2 = √(r31 * r42) * (fo - τ) / 2
    νr4 = X2 > zero(TT)

    radial_roots = roots(pix)
    r1, r2, _, _ = radial_roots
    r21, r31c, r32, r41c, r42c, _ = _get_root_diffs(radial_roots...)
    r1, r2, r21 = real.((r1, r2, r21))
    A = √abs(r32 * r42c)
    B = √abs(r31c * r41c)
    X3 = √(A * B) * real(fo - τ)
    νr2 = X3 > zero(TT)

    _, r2c, _, r4c = radial_roots
    a1 = abs(imag(r4c))
    a2 = abs(imag(r2c))
    b1 = real(r4c)
    b2 = real(r2c)
    C = sqrt((a1 - a2)^2 + (b1 - b2)^2)
    D = sqrt((a1 + a2)^2 + (b1 - b2)^2)
    X4 = (C + D) / TT(2) * (fo - τ)
    νr0 = X4 > zero(TT)

    is4 = numreals == 4
    is2 = numreals == 2
    rs = Base.ifelse(is4, rs4, Base.ifelse(is2, rs2, rs0))
    νri = Base.ifelse(is4, νr4, Base.ifelse(is2, νr2, νr0))
    validi = Base.ifelse(is4, valid4, Base.ifelse(is2, valid2, valid0))
    issuccess = validi > zero(validi)

    return rs, νri, numreals, issuccess
end


#for type_list in Iterators.product((types for _ in 1:18)...)
#    if all(x-> x == Any, type_list)
#        continue
#    end
#
#    symbols = (Symbol(:T, i) for i = 1:18) 
#    for (i, symbol) in enumerate(symbols)
#        @eval $symbol = $type_list[$i]
#    end
#
#    @eval function Krang.emission_radius(
#        pix::Krang.SlowLightIntensityPixel{$T1,$T2,$T3,$T4,$T5,$T6,$T7,$T8,$T9,$T10,$T11,$T12,$T13,$T14,$T15,$T16,$T17},
#        τ::$T18,
#    ) 
#        return _reactant_get_radial_roots(met, τ)
#    end
#end

Reactant.@reactant_overlay function Krang.emission_radius(pix::Krang.AbstractPixel, τ) 
    met = metric(pix)
    a = met.spin
    T = typeof(a)
    ans = zero(T)
    νr = true

    rh = one(T) + √(one(T) - a^2)

    numreals = sum(_isreal2, roots(pix))

    @trace if numreals == 4 #case 1 & 2
        ans, νr, issuccess = _rs_case1_and_2(pix, rh, τ)
    elseif numreals == 2 #case3
        ans, νr, issuccess = _rs_case3(pix, rh, τ)
    else #case 4
        ans, νr, issuccess = _rs_case4(pix, rh, τ)
    end
    return ans, νr, numreals, issuccess
end



Reactant.@reactant_overlay function Krang.emission_radius(
    pix::Krang.SlowLightIntensityPixel,
    θs,
    isindir,
    n,
) 
    α, β = screen_coordinate(pix)
    θo = inclination(pix)
    met = metric(pix)
    TT = typeof(met.spin)
    piT = TT(Base.pi)
    isincone = ((θo <= θs) & (θs <= (piT - θo))) | (((piT - θo) <= θs) & (θs <= θo))

    err_rs = zero(TT)
    err_νr = true
    err_νθ = true
    err_numreals = 0
    err_success = false

    αmin = αboundary(met, θs)
    αoutside = abs(α) >= (αmin + eps(TT))
    βbound = Base.ifelse((!isincone) & αoutside, βboundary(met, α, θo, θs), zero(TT))
    blocked = (!isincone) & ((abs(β) + eps(TT)) < βbound)

    τ, τvalid = _reactant_Gθ_minotime(pix, θs, isindir, n)
    rs, νr, numreals, rssuccess = Krang.emission_radius(pix, τ)
    νθ = Base.ifelse(isincone, (θo > θs) ⊻ (n % 2 == 1), !isindir)
    issuccess = (!blocked) & τvalid & rssuccess

    rs = Base.ifelse(issuccess, rs, err_rs)
    νr = Base.ifelse(issuccess, νr, err_νr)
    νθ = Base.ifelse(issuccess, νθ, err_νθ)
    numreals = Base.ifelse(issuccess, numreals, err_numreals)

    return rs, νr, νθ, numreals, issuccess
end

Reactant.@reactant_overlay function Krang.Gϕ(pix::Krang.SlowLightIntensityPixel, θs, isindir, n)
    ans, valid = _reactant_Gϕ_minotime(pix, θs, isindir, n)
    TT = typeof(metric(pix).spin)
    return ans, zero(TT), zero(TT), zero(TT), false, valid
end

Reactant.@reactant_overlay function Krang.Gt(pix::Krang.SlowLightIntensityPixel, θs, isindir, n)
    ans, valid = _reactant_Gt_minotime(pix, θs, isindir, n)
    TT = typeof(metric(pix).spin)
    return ans, zero(TT), zero(TT), zero(TT), false, valid
end

Reactant.@reactant_overlay function Krang.emission_coordinates(
    pix::Krang.SlowLightIntensityPixel,
    θs,
    isindir,
    n,
)
    α, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)
    TT = typeof(met.spin)
    cosθs = cos(θs)
    cosθo = cos(θo)
    err_return = (zero(TT), zero(TT), zero(TT), false, false, false)

    αmin = αboundary(met, θs)
    βbound = Base.ifelse(abs(α) >= (αmin + eps(TT)), βboundary(met, α, θo, θs), zero(TT))
    blocked = (cosθs > abs(cosθo)) & ((abs(β) + eps(TT)) < βbound)

    τ, τvalid = _reactant_Gθ_minotime(pix, θs, isindir, n)
    rs, νr, _, rssuccess = Krang.emission_radius(pix, τ)

    a = met.spin
    λtemp = Krang.λ(pix)
    I0, I1, I2, Ip, Im = Krang.radial_integrals(pix, rs, τ, νr)

    rp = one(TT) + √(one(TT) - a^2)
    rm = one(TT) - √(one(TT) - a^2)
    It = 4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) +
         4 * I0 + 2 * I1 + I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, Gϕvalid = _reactant_Gϕ_minotime(pix, θs, isindir, n)
    Gttemp, Gtvalid = _reactant_Gt_minotime(pix, θs, isindir, n)

    emission_azimuth = Iϕ + λtemp * Gϕtemp
    emission_time_regularized = It + a^2 * Gttemp

    νθ = Base.ifelse(abs(cosθs) < abs(cosθo), (n % 2 == 1) ⊻ (θo > θs), !isindir ⊻ (θs > TT(Base.pi / 2)))
    issuccess = (!blocked) & τvalid & rssuccess & Gϕvalid & Gtvalid & !isnan(τ)

    t = Base.ifelse(issuccess, emission_time_regularized, zero(TT))
    r = Base.ifelse(issuccess, rs, zero(TT))
    ϕ = Base.ifelse(issuccess, emission_azimuth, zero(TT))
    νr = Base.ifelse(issuccess, νr, false)
    νθ = Base.ifelse(issuccess, νθ, false)

    return t, r, ϕ, νr, νθ, issuccess
end

Reactant.@reactant_overlay function Krang.emission_coordinates(
    pix::Krang.SlowLightIntensityPixel,
    τ,
)
    α, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)
    TT = typeof(met.spin)

    θs, _, _, _, n, isindir = Krang.emission_inclination(pix, τ)
    αmin = αboundary(met, θs)
    βbound = Base.ifelse(abs(α) >= (αmin + eps(TT)), βboundary(met, α, θo, θs), zero(TT))
    blocked = (cos(θs) > abs(cos(θo))) & ((abs(β) + eps(TT)) < βbound)

    a = met.spin
    λtemp = Krang.λ(pix)
    rs, νr, _, rssuccess = Krang.emission_radius(pix, τ)
    I0, I1, I2, Ip, Im = Krang.radial_integrals(pix, rs, τ, νr)

    rp = one(TT) + √(one(TT) - a^2)
    rm = one(TT) - √(one(TT) - a^2)

    It = 4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) +
         4 * I0 + 2 * I1 + I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, Gϕvalid = _reactant_Gϕ_minotime(pix, θs, isindir, n)
    Gttemp, Gtvalid = _reactant_Gt_minotime(pix, θs, isindir, n)

    emission_azimuth = Iϕ + λtemp * Gϕtemp
    emission_time_regularized = It + a^2 * Gttemp
    νθ = Base.ifelse(abs(cos(θs)) < abs(cos(θo)), (n % 2 == 1) ⊻ (θo > θs), !isindir ⊻ (θs > TT(Base.pi / 2)))
    issuccess = (!blocked) & rssuccess & Gϕvalid & Gtvalid

    t = Base.ifelse(issuccess, emission_time_regularized, zero(TT))
    r = Base.ifelse(issuccess, rs, zero(TT))
    θ = Base.ifelse(issuccess, θs, zero(TT))
    ϕ = Base.ifelse(issuccess, emission_azimuth, zero(TT))
    νr = Base.ifelse(issuccess, νr, true)
    νθ = Base.ifelse(issuccess, νθ, true)

    return t, r, θ, ϕ, νr, νθ, issuccess
end

Reactant.@reactant_overlay function Krang._raytrace(
    observation,
    pix::Krang.SlowLightIntensityPixel,
    mesh::Krang.Mesh{<:Krang.ConeGeometry,<:Krang.AbstractMaterial};
    res,
)
    geometry = mesh.geometry
    material = mesh.material
    θs = geometry.opening_angle
    subimgs = material.subimgs

    isindir = false
    for _ = 1:2
        isindir ⊻= true
        for n in subimgs
            rs = zero(typeof(metric(pix).spin))
            issuccess = false
            rs, issuccess = _reactant_emission_radius_theta(pix, θs, isindir, n)
            θs_tr = zero(rs) + θs
            intersection = Intersection(zero(rs), rs, θs_tr, zero(rs), true, true)

            cond = issuccess & (horizon(metric(pix)) < rs) & (rs < T(Inf))
            contribution = material(pix, intersection)
            observation += Base.ifelse(cond, contribution, zero(contribution))
        end
    end

    return observation
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

function _reactant_intensity_pixel(met::Krang.Kerr, α, β, θo)
    TT = typeof(met.spin)
    tempη = Krang.η(met, α, β, θo)
    tempλ = Krang.λ(met, α, θo)
    roots = zeros(Complex{typeof(tempλ)}, 4)
    roots = get_radial_roots(met, tempη, tempλ)
    r1, r2, r3, r4 = roots
    numreals = sum(_isreal2, roots)

    roots = Base.ifelse((numreals == 2) & (abs(imag(r4)) < sqrt(eps(TT))), (r1, r4, r2, r3), roots)

    I0_inf = Krang.Ir_inf(met, roots)
    τ_total = Krang.total_mino_time(met, roots)
    absGθo_Gθhat = Krang._absGθo_Gθhat(met, θo, tempη, tempλ)

    return Krang.IntensityPixel(
        met,
        (α, β),
        roots,
        I0_inf,
        τ_total,
        absGθo_Gθhat,
        θo,
        tempη,
        tempλ,
    )
end

Reactant.@reactant_overlay function Krang.SlowLightIntensityPixel(
    met::Krang.Kerr,
    α,
    β,
    θo,
) 
    return _reactant_slow_light_intensity_pixel(met, α, β, θo)
end


#for (typespin, αT, βT, θT) in Iterators.product(types, types, types, types)
#    if typespin == αT == βT == θT == Any
#        continue
#    end
#
#    @eval function Krang.SlowLightIntensityPixel(
#        met::Krang.Kerr{$typespin},
#        α::$αT,
#        β::$βT,
#        θo::$θT,
#    ) 
#        return _reactant_slow_light_intensity_pixel(met, α, β, θo)
#    end
#
#    @eval function Krang.IntensityPixel(
#        met::Krang.Kerr{$typespin},
#        α::$αT,
#        β::$βT,
#        θo::$θT,
#    )
#        return _reactant_intensity_pixel(met, α, β, θo)
#    end
#end

end
