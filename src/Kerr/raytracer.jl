include("./utils.jl")
# Follows the Formalism of Gralla & Lupsasca (https://arxiv.org/pdf/1910.12881.pdf)
export emission_radius, emission_inclination, emission_coordinates_fast_light, emission_coordinates

"""
Emission radius for point originating at inclination θs whose nth order image appears at the screen coordinate (`α`, `β`). 
Returns NaN if the emission coordinates do not exist for that screen coordinate.

# Arguments

- `metric` : Kerr{T} metric 
- `α` : Horizontal Bardeen screen coordinate
- `β` : Vertical Bardeen screen coordinate
- `θs` : Emission inclination
- `θo` : Observer inclination
- `isindir` : Is emission to observer direct or indirect
- `n` : Image index
"""
function emission_radius(metric::Kerr{T}, α, β, θs, θo, isindir, n)::Tuple{Any, Bool, Bool} where {T}
    if cos(θs) > abs(cos(θo))
        αmin = αboundary(metric, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(metric, α, θo, θs) : zero(T))
        ((abs(β) + eps(T)) < βbound) && return (T(NaN), true, true)
    end

    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    τ, _, _, _ = _Gθ(metric, sign(β), θs, θo, isindir, n, ηtemp, λtemp)
    (isnan(τ) || isinf(τ)) && return (T(NaN), true, true)

    # is θ̇s increasing or decreasing?
    νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
    # is ṙs increasing or decreasing?
    rs, νr, _ = _rs(metric, τ, get_radial_roots(metric, ηtemp, λtemp))

    return rs, νr, νθ
end

"""
Emission radius for point originating at at Mino time τ whose image appears at the screen coordinate (`α`, `β`). 
Returns NaN if the emission coordinates do not exist for that screen coordinate.

# Arguments

-`metric` : Kerr{T} metric
-`α` : Horizontal Bardeen screen coordinate
-`β` : Vertical Bardeen screen coordinate
-`θs` : Emission inclination
-`θo` : Observer inclination
-`isindir` : Normalized angular momentum
-`n` : Image index
"""
function emission_radius(metric::Kerr{T}, α, β, τ, θo) where {T}
    return _rs(metric, τ, get_radial_roots(metric, η(metric, α, β, θo), λ(metric, α, θo)))
end

"""
Emission inclination for point originating at inclination rs whose nth order image appears at screen coordinate (`α`, `β`).

# Arguments

- `metric` : Kerr{T} metric
- `α` : Horizontal Bardeen Screen Coordinates
- `β` : Vertical Bardeen Screen Coordinates
- `θo` : Observer Inclination
- `τ` : Mino Time
"""
function emission_inclination(metric::Kerr{T}, α, β, θo, rs, νr) where {T}
    τ = Ir(metric, νr, rs, η(metric, α, β, θo), λ(metric, α, θo))
    return _θs(metric, sign(β), θo, η(metric, α, β, θo), λ(metric, α, θo), τ)
end

"""
Emission inclination for point at Mino time τ whose image appears at screen coordinate (`α`, `β`).

# Arguments

- `metric` : Kerr{T} metric
- `α` : Horizontal Bardeen Screen Coordinates
- `β` : Vertical Bardeen Screen Coordinates
- `θo` : Observer Inclination
- `τ` : Mino Time
"""
function emission_inclination(metric::Kerr{T}, α, β, τ, θo) where {T}
    return _θs(metric, sign(β), θo, η(metric, α, β, θo), λ(metric, α, θo), τ)
end

"""
Emission radius and azimuthal angle for point originating at inclination θs whose nth order image appears at the screen coordinate (`α`, `β`). 
for an observer located at inclination θo.

# Arguments

- `metric` : Kerr{T} metric
- `α` : Horizontal Bardeen Screen Coordinates
- `β` : Vertical Bardeen Screen Coordinates
- `θs` : Emission Inclination
- `θo` : Observer Inclination
- `isindir` : Whether emission to observer is direct or indirect
- `n` : Image index
"""
function emission_coordinates_fast_light_inclination_parameterized(metric::Kerr{T}, α, β, θs, θo, isindir, n)::Tuple{Any, Any, Any, Bool, Bool} where {T}
    if cos(θs) > abs(cos(θo))
        αmin = αboundary(metric, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(metric, α, θo, θs) : zero(T))
        if (abs(β) + eps(T)) < βbound
            return T(NaN), T(NaN), T(NaN), false, false
        end
    end

    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    τ, _, _, _ = _Gθ(metric, sign(β), θs, θo, isindir, n, ηtemp, λtemp)
    if isnan(τ)
        return T(NaN), T(NaN), T(NaN), false, false
    end

    roots = get_radial_roots(metric, ηtemp, λtemp)
    emission_radius, νr, _ = _rs(metric, τ, roots)
    emission_azimuth = _ϕs(metric, α, β, θs, θo, emission_radius, τ, νr, roots, isindir, n)
    νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir

    return emission_radius, θs, emission_azimuth, νr, νθ
end

"""
Emission radius and azimuthal angle for point originating at inclination θs whose nth order image appears at the screen coordinate (`α`, `β`).
for an observer located at inclination θo.

# Arguments

- `metric` : Kerr{T} metric
- `α` : Horizontal Bardeen Screen Coordinates
- `β` : Vertical Bardeen Screen Coordinates
- `θs` : Emission Inclination
- `θo` : Observer Inclination
- `isindir` : Whether emission to observer is direct or indirect
- `n` : Image index
"""
function emission_coordinates(metric::Kerr{T}, α, β, θs, θo, isindir, n) where {T}
    if cos(θs) > abs(cos(θo))
        αmin = αboundary(metric, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(metric, α, θo, θs) : zero(T))

        if (abs(β) + eps(T)) < βbound
            return T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN)
        end
    end

    a = metric.spin

    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    τ, _, _, _ = _Gθ(metric, sign(β), θs, θo, isindir, n, ηtemp, λtemp)

    νθ = abs(cos(θs)) < abs(cos(θo)) ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    #isindir = !(((τ + sign(β)*τo - (νθ ? 1 : -1)*τs)/τhat) ≈ n) # TODO: Why is this here
    if (abs(cos(θs)) < abs(cos(θo)))
        isindir = ((sign(β) > 0) ⊻ (θo > π / 2))
    end
    if isnan(τ)
        return T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN)
        #else
        #println(νθ ≈ ((sign(cos(θs)) > 0)⊻(τs2< 0)))
    end

    roots = get_radial_roots(metric, ηtemp, λtemp)
    emission_radius, νr, _ = _rs(metric, τ, roots)
    numreals = sum(_isreal2.(roots))
    if numreals == 4
        roots = real.(roots)
        I0, I1, I2, Ip, Im = radial_integrals_case2(metric, emission_radius, roots, τ, νr)
    elseif numreals == 2
        if abs(imag(roots[4])) < sqrt(eps(T))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        I0, I1, I2, Ip, Im = radial_integrals_case3(metric, emission_radius, roots, τ)
    else
        I0, I1, I2, Ip, Im = radial_integrals_case4(metric, emission_radius, roots, τ)
    end

    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    It = 4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) + 4 * I0 + 2 * I1 + I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, _, _, _, _ = Gϕ(metric, α, β, θs, θo, isindir, n)
    Gttemp, _, _, _, _ = Gt(metric, α, β, θs, θo, isindir, n)

    emission_azimuth = (T(π) - Iϕ - λtemp * Gϕtemp + 4π) % T(2π)
    emission_time_regularized = (zero(T) + It + a^2 * Gttemp)

    # is θ̇s increasing or decreasing?
    νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir

    return emission_time_regularized, emission_radius, θs, emission_azimuth, νr, νθ
end

"""
Raytrace a point that appears at the screen coordinate (`α`, `β`) for an observer located at inclination θo

# Arguments

- `metric` : Kerr{T} metric
- `α` : Horizontal Bardeen Screen Coordinates
- `β` : Vertical Bardeen Screen Coordinates
- `θo` : Observer Inclination
- `τ` : Mino Time
"""
function raytrace(metric::Kerr{T}, α, β, θo, τ) where {T}
    θs, τs, _, _, n, isindir = emission_inclination(metric, α, β, θo, τ)
    νθ = τs > 0

    if cos(θs) > abs(cos(θo))
        αmin = αboundary(metric, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(metric, α, θo, θs) : zero(T))

        if (abs(β) + eps(T)) < βbound
            return T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN)
        end
    end

    a = metric.spin

    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    roots = get_radial_roots(metric, ηtemp, λtemp)
    emission_radius, νr, _ = _rs(metric, τ, roots)
    numreals = sum(_isreal2.(roots))
    if numreals == 4
        roots = real.(roots)
        I0, I1, I2, Ip, Im = radial_integrals_case2(metric, emission_radius, roots, τ, νr)
    elseif numreals == 2
        if abs(imag(roots[4])) < sqrt(eps(T))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        I0, I1, I2, Ip, Im = radial_integrals_case3(metric, emission_radius, roots, τ)
    else
        I0, I1, I2, Ip, Im = radial_integrals_case4(metric, emission_radius, roots, τ)
    end

    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    It = 4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) + 4 * I0 + 2 * I1 + I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, _, _, _, _ = Gϕ(metric, α, β, θs, θo, isindir, n)
    Gttemp, _, _, _, _ = Gt(metric, α, β, θs, θo, isindir, n)

    emission_azimuth = (-Iϕ - λtemp * Gϕtemp)# % T(2π)
    emission_time_regularized = (It + a^2 * Gttemp)


    return emission_time_regularized, emission_radius, θs, emission_azimuth, νr, νθ
end


function _rs(metric::Kerr{T}, τ, roots::NTuple{4}) where {T}
    a = metric.spin
    ans = zero(T)
    νr = true

    rh = one(T) + √(one(T) - a^2)

    numreals = sum(_isreal2.(roots))

    if numreals == 4 #case 1 & 2
        ans, νr = _rs_case1_and_2(metric, real.(roots), rh, τ)
    elseif numreals == 2 #case3
        if abs(imag(roots[4])) < sqrt(eps(T))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        ans, νr = _rs_case3(metric, roots, τ)
    else #case 4
        ans, νr = _rs_case4(metric, roots, rh, τ)
    end
    return ans, νr, numreals
end

function _rs_case1_and_2(metric::Kerr{T}, roots::NTuple{4}, rh, τ) where {T}
    _, _, r3, r4 = roots
    _, r31, r32, r41, r42, _ = _get_root_diffs(roots...)

    horizonτ, fo = Ir_case1_and_2(metric, roots, rh, true)
    r4 < rh && τ > horizonτ && return T(NaN), true# invalid case2

    k = (r32 * r41) / (r31 * r42)
    X2 = √(r31 * r42) * (fo - τ) / T(2)
    if τ > 2fo
        return T(NaN), true
    end
    sn = r41 * JacobiElliptic.sn(X2, k)^2
    return (r31 * r4 - r3 * sn) / (r31 - sn), X2 > zero(T)
end

function _rs_case3(::Kerr{T}, roots::NTuple{4,Complex}, τ) where {T}
    r1, r2, _, _ = roots
    root_diffs = _get_root_diffs(roots...)
    r21, r31, r32, r41, r42, _ = root_diffs
    r1, r2, r21 = real.((r1, r2, r21))

    A = √abs(r32 * r42)
    B = √abs(r31 * r41)
    k = (((A + B)^2 - r21^2) / (T(4) * A * B))

    fo = JacobiElliptic.F(acos(clamp((A - B) / (A + B), -1, 1)), k)
    X3 = real(fo - √(A * B) * τ)
    if X3 < zero(T)
        return T(NaN), true
    end
    cn = JacobiElliptic.cn(X3, k)
    num = -A * r1 + B * r2 + (A * r1 + B * r2) * cn
    den = -A + B + (A + B) * cn

    return real(num / den), X3 > zero(T)
end

function _rs_case4(metric::Kerr{T}, roots::NTuple{4,Complex}, rh, τ) where {T}
    r1, _, _, r4 = roots
    root_diffs = _get_root_diffs(roots...)
    _, r31, r32, r41, r42, _ = root_diffs

    τ > Ir_case4(metric, roots, root_diffs, rh)[1] && return T(NaN), true
    a2 = abs(imag(r1))
    b1 = real(r4)
    C = √abs(r31 * r42)
    D = √abs(r32 * r41)
    k4 = 4 * C * D / (C + D)^2

    go = √max((T(4) * a2^2 - (C - D)^2) / ((C + D)^2 - T(4) * a2^2), zero(T))
    fo = T(2) / (C + D) * JacobiElliptic.F(T(π / 2) + atan(go), k4)
    X4 = (C + D) / T(2) * (fo - τ)
    num = go - JacobiElliptic.sc(X4, k4)
    den = one(T) + go * JacobiElliptic.sc(X4, k4)

    return -(a2 * num / den + b1), X4 > zero(T)
end

function θs(metric::Kerr{T}, α, β, θo, τ) where {T}
    return _θs(metric, sign(β), θo, η(metric, α, β, θo), λ(metric, α, θo), τ)
end

function _θs(metric::Kerr{T}, signβ, θo, η, λ, τ) where {T}
    a = metric.spin
    Ghat_2, isvortical = zero(T), η < zero(T)

    Δθ = T(0.5) * (one(T) - (η + λ^2) / a^2)
    up = Δθ + √(Δθ^2 + η / a^2)
    um = Δθ - √(Δθ^2 + η / a^2)
    m = up / um
    k = m

    #isvortical = η < 0.
    ans, k, argo, τs, τo, τhat = zero(T), zero(T), zero(T), zero(T), zero(T), zero(T)
    tempfac = one(T) / √abs(um * a^2)
    if isvortical
        argo = (cos(θo)^2 - um) / (up - um)
        k = one(T) - m
        tempfac = inv(√abs(um * a^2))
        τo = tempfac * JacobiElliptic.F(asin(√argo), k)
        Ghat_2 = tempfac * JacobiElliptic.K(k)
        τhat = 2Ghat_2
        Δτtemp = (τ % τhat + (θo > T(π / 2) ? -one(T) : one(T)) * signβ * τo)
        n = floor(τ / τhat)
        #absτs = abs(argmin(abs, [(-one(T))^n * signβ * (Ghat - Δτtemp), (-one(T))^n * signβ * Δτtemp]))
        absτs = abs(argmin(abs, [τhat - Δτtemp, Δτtemp]))
        τs = (θo > T(π / 2) ? -one(T) : one(T)) * absτs

        argr = (JacobiElliptic.sn(absτs / tempfac, k))^2
        ans = acos((θo > T(π / 2) ? -one(T) : one(T)) * √((up - um) * argr + um))
        τo = (ans > π / 2 ? -1 : 1) * τo
    else
        argo = clamp(cos(θo) / √(up), -one(T), one(T))
        k = m
        tempfac = inv(√abs(um * a^2))
        τo = tempfac * JacobiElliptic.F(asin(argo), k)
        Ghat_2 = tempfac * JacobiElliptic.K(k)
        τhat = Ghat_2 + Ghat_2
        Δτtemp = (τ % τhat + signβ * τo)
        n = floor(τ / τhat)
        τs = argmin(abs, [(-one(T))^n * signβ * (τhat - Δτtemp), (-one(T))^n * signβ * Δτtemp])
        #τs = abs(argmin(abs, [(Ghat - Δτtemp), Δτtemp]))
        newargs = JacobiElliptic.sn(τs / tempfac, k)
        ans = acos(√up * newargs)
    end

    # TODO: Come up with a more elloquent solution for this
    isincone = abs(cos(ans)) < abs(cos(θo))
    if isincone
        νθ = (n % 2 == 1) ⊻ (θo > ans)
        τ1 = (n + 1) * τhat - signβ * τo + (νθ ? 1 : -1) * τs
        isindir = isapprox(τ1, τ, rtol=1e-5)
    elseif ans < π / 2
        τ1 = (n + 1) * τhat - signβ * τo - τs
        isindir = isapprox(τ1, τ, rtol=1e-5)
    else
        τ1 = (n + 1) * τhat - signβ * τo + τs
        isindir = isapprox(τ1, τ, rtol=1e-5)
    end

    return ans, τs, τo, τhat, n, isindir
end

function _ϕs(metric::Kerr{T}, α, β, θs, θo, rs, τ, νr, roots, isindir, n) where {T}
    λtemp = λ(metric, α, θo)
    numreals = sum(_isreal2.(roots))
    root_diffs = _get_root_diffs(roots...)
    Iϕ = zero(T)
    if numreals == 4 #case2
        roots = real.(roots)
        root_diffs = real.(root_diffs)
        Iϕ = Iϕ_case2(metric, roots, root_diffs, rs, τ, νr, λtemp)
    elseif numreals == 2 #case3
        if abs(imag(roots[4])) < eps(T)
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        Iϕ = Iϕ_case3(metric, roots, root_diffs, rs, τ, λtemp)
    else
        Iϕ = Iϕ_case4(metric, roots, root_diffs, rs, τ, λtemp)
    end
    Iϕ == T(Inf) && return zero(T)

    Gϕtemp, _, _, _ = Gϕ(metric, α, β, θs, θo, isindir, n)
    Gϕtemp == T(Inf) && return zero(T)

    if isnan(Iϕ) || Iϕ == T(Inf)
        return zero(T)
    end
    if isnan(Gϕtemp) || Gϕtemp == T(Inf)
        return zero(T)
    end
    return (T(π) - Iϕ - λtemp * Gϕtemp + 4π) % T(2π)
end

