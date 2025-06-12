# Follows the Formalism of Gralla & Lupsasca (https://arxiv.org/pdf/1910.12881.pdf).
# Functions generically return 0 when the emission coordinates do not exist for a given screen coordinate to play nice 
# with Enzyme's AD.

export emission_radius,
    emission_inclination,
    emission_coordinates_fast_light,
    emission_coordinates,
    emission_coordinates!

function emission_coordinates!(
    coordinates::Array{T,3},
    camera::AbstractCamera,
    τ::T,
) where {T}
    for I in CartesianIndices(camera.screen.pixels)
        coordinates[:, I] .= emission_coordinates(camera.screen.pixels[I], τ)[1:4]
    end
end

"""
Emission radius for point originating at inclination θs whose nth order image appears at the screen coordinate (`α`, `β`). 
Returns 0 if the emission coordinates do not exist for that screen coordinate.

# Arguments

- `pix` : Pixel information 
- `θs` : Emission inclination
- `isindir` : Is emission to observer direct or indirect
- `n` : Image index
"""
function emission_radius(
    pix::Krang.AbstractPixel,
    θs::T,
    isindir,
    n,
)::Tuple{T,Bool,Bool,Int,Bool} where {T}
    α, β = screen_coordinate(pix)
    θo = inclination(pix)
    met = metric(pix)
    isincone = θo ≤ θs ≤ (π - θo) || (π - θo) ≤ θs ≤ θo
    if !isincone
        αmin = αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(met, α, θo, θs) : zero(T))
        ((abs(β) + eps(T)) < βbound) && return zero(T), true, true, 0, false
    end

    τ, _, _, _, _, issuccess = Gθ(pix, θs, isindir, n)
    issuccess || return zero(T), true, true, 0, false

    # νθ is θ̇s increasing or decreasing?
    νθ = isincone ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
    if isincone
        νθ = (θo > θs) ⊻ (n % 2 == 1)
    end

    ans, νr, numreals, issuccess = emission_radius(pix, τ)
    return ans, νr, νθ, numreals, issuccess
end

"""
Emission radius for point originating at at Mino time τ whose image appears at the screen coordinate (`α`, `β`). 
Returns 0 if the emission coordinates do not exist for that screen coordinate.

# Arguments

-`pix` : Pixel information
- `τ` : Mino time
"""
function emission_radius(pix::AbstractPixel, τ::T)::Tuple{T,Bool,Int,Bool} where {T}
    met = metric(pix)
    a = met.spin
    ans = zero(T)
    νr = true

    rh = one(T) + √(one(T) - a^2)

    numreals = sum(Int.(_isreal2.(roots(pix))))

    if numreals == 4 #case 1 & 2
        ans, νr, issuccess = _rs_case1_and_2(pix, rh, τ)
    elseif numreals == 2 #case3
        ans, νr, issuccess = _rs_case3(pix, rh, τ)
    else #case 4
        ans, νr, issuccess = _rs_case4(pix, rh, τ)
    end
    return ans, νr, numreals, issuccess
end

"""
Emission inclination for point originating at inclination rs whose nth order image appears at screen coordinate (`α`, `β`).

# Arguments

- `pix` : Pixel information
- `rs` : Emission radius
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
"""
function emission_inclination(pix::AbstractPixel, rs, νr)
    return emission_inclination(pix, Ir(pix, νr, rs))
end

"""
Emission inclination for point at Mino time τ whose image appears at screen coordinate (`α`, `β`).

# Arguments

- `pix` : Pixel information
- `τ` : Mino Time
"""
function emission_inclination(pix::AbstractPixel, τ)
    met = metric(pix)
    θo = inclination(pix)
    _, β = pix.screen_coordinate

    return _θs(met, sign(β), θo, η(pix), λ(pix), τ)
end

"""
Emission azimuth for point at Mino time τ whose image appears at screen coordinate (`α`, `β`).

# Arguments

- `pix` : Pixel information
- `θs` : Emission Inclination
- `rs` :Emission radius
- `τ` : Mino Time
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
"""
function emission_azimuth(pix::AbstractPixel, θs, rs, τ::T, νr, isindir, n) where {T}
    met = metric(pix)
    θo = inclination(pix)

    α, _ = pix.screen_coordinate
    λtemp = λ(pix)
    Iϕ = Krang.Iϕ(pix, rs, τ, νr)
    (isnan(Iϕ) || !isfinite(Iϕ)) && return Iϕ

    Gϕtemp, _, _, _ = Gϕ(pix, θs, isindir, n)
    (isnan(Gϕtemp) || !isfinite(Gϕtemp)) && return Iϕ

    return -(Iϕ + λtemp * Gϕtemp + T(π / 2))
end

"""
Emission radius and azimuthal angle for point originating at inclination θs whose nth order image appears at the screen 
coordinate (`α`, `β`) for an observer located at inclination θo.

# Arguments

- `pix` : Pixel information
- `θs` : Emission Inclination
- `isindir` : Whether emission to observer is direct or indirect
- `n` : Image index
"""
function emission_coordinates_fast_light(pix::AbstractPixel, θs::T, isindir, n) where {T}
    # NB: I do not return θs since it is already known, and returning it might encourage bad coding errors
    α, β = screen_coordinate(pix)
    θo = inclination(pix)
    met = metric(pix)
    isincone = θo ≤ θs ≤ (π - θo) || (π - θo) ≤ θs ≤ θo
    if !isincone
        αmin = αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(met, α, θo, θs) : zero(T))
        ((abs(β) + eps(T)) < βbound) && return zero(T), zero(T), false, false, false
    end

    τ, _, _, _, _, issuccess = Gθ(pix, θs, isindir, n)
    issuccess || return (zero(T), zero(T), false, false, false)

    # is θ̇s increasing or decreasing?
    νθ = isincone ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
    if isincone
        νθ = (θo > θs) ⊻ (n % 2 == 1)
    end

    rs, νr, _, issuccess = emission_radius(pix, τ)
    issuccess || return (zero(T), zero(T), false, false, false)

    ϕs = emission_azimuth(pix, θs, rs, τ, νr, isindir, n)

    return rs, ϕs, νr, νθ, issuccess
end

"""
Ray trace a point that appears at the screen coordinate (`α`, `β`) for an observer located at inclination θo

# Arguments

- `pix` : Pixel information
- `τ` : Mino Time
"""
function emission_coordinates_fast_light(pix::AbstractPixel, τ::T) where {T}
    α, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)

    θs, τs, _, _, n, isindir = emission_inclination(pix, τ)
    νθ = τs > 0

    if cos(θs) > abs(cos(θo))
        αmin = αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(met, α, θo, θs) : zero(T))

        if (abs(β) + eps(T)) < βbound
            return zero(T), zero(T), zero(T), true, true, false
        end
    end

    a = met.spin

    λtemp = λ(pix)
    rs, νr, _, issuccess = emission_radius(pix, τ)
    issuccess || return zero(T), zero(T), zero(T), true, true, false
    _, _, _, Ip, Im = radial_integrals(pix, rs, τ, νr)

    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, _, _, _, _ = Gϕ(pix, θs, isindir, n)

    emission_azimuth = -(Iϕ + λtemp * Gϕtemp + T(π / 2))

    νθ = abs(cos(θs)) < abs(cos(θo)) ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    return rs, θs, emission_azimuth, νr, νθ, issuccess
end


"""
Emission radius and azimuthal angle for point originating at inclination θs whose nth order image appears at the screen
coordinate (`α`, `β`) for an observer located at inclination θo.

# Arguments

- `pix` : Pixel information
- `θs` : Emission Inclination
- `isindir` : Whether emission to observer is direct or indirect
- `n` : Image index
"""
function emission_coordinates(pix::AbstractPixel, θs::T, isindir, n) where {T}
    α, β = pix.screen_coordinate
    met = metric(pix)
    θo = inclination(pix)
    cosθs = cos(θs)
    cosθo = cos(θo)
    if cos(θs) > abs(cosθo)
        αmin = αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(met, α, θo, θs) : zero(T))

        if (abs(β) + eps(T)) < βbound
            return zero(T), zero(T), zero(T), false, false, false
        end
    end

    τ, _, _, _, _, issuccess = Gθ(pix, θs, isindir, n)
    issuccess || return zero(T), zero(T), zero(T), false, false, false

    rs, νr, _, issuccess = emission_radius(pix, τ)
    issuccess || return zero(T), zero(T), zero(T), false, false, false

    a = met.spin

    λtemp = λ(pix)

    νθ = abs(cosθs) < abs(cosθo) ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    #isindir = !(((τ + sign(β)*τo - (νθ ? 1 : -1)*τs)/τhat) ≈ n) # TODO: Why is this here
    if (abs(cosθs) < abs(cosθo))
        isindir = ((sign(β) > 0) ⊻ (θo > T(π / 2)))
    end
    if isnan(τ)
        return zero(T), zero(T), zero(T), false, false, false
    end

    I0, I1, I2, Ip, Im = radial_integrals(pix, rs, τ, νr)

    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    It =
        4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) +
        4 * I0 +
        2 * I1 +
        I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, _, _, _, _ = Gϕ(pix, θs, isindir, n)
    Gttemp, _, _, _, _ = Gt(pix, θs, isindir, n)

    emission_azimuth = -(Iϕ + λtemp * Gϕtemp + T(π / 2))
    emission_time_regularized = (zero(T) + It + a^2 * Gttemp)

    # is θ̇s increasing or decreasing?
    νθ = cosθs < abs(cosθo) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir

    return emission_time_regularized, rs, emission_azimuth, νr, νθ, issuccess
end

"""
Ray trace a point that appears at the screen coordinate (`α`, `β`) for an observer located at inclination θo

# Arguments

- `pix` : Pixel information
- `τ` : Mino Time
"""
function emission_coordinates(pix::AbstractPixel, τ::T) where {T}
    α, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)

    θs, _, _, _, n, isindir = emission_inclination(pix, τ)

    if cos(θs) > abs(cos(θo))
        αmin = αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(met, α, θo, θs) : zero(T))

        if (abs(β) + eps(T)) < βbound
            return zero(T), zero(T), zero(T), zero(T), true, true, false
        end
    end

    a = met.spin

    λtemp = λ(pix)
    rs, νr, _, issuccess = emission_radius(pix, τ)
    issuccess || return zero(T), zero(T), zero(T), zero(T), true, true, false
    I0, I1, I2, Ip, Im = radial_integrals(pix, rs, τ, νr)

    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    It =
        4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) +
        4 * I0 +
        2 * I1 +
        I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, _, _, _, _ = Gϕ(pix, θs, isindir, n)
    Gttemp, _, _, _, _ = Gt(pix, θs, isindir, n)

    emission_azimuth = -(Iϕ + λtemp * Gϕtemp + T(π / 2))
    emission_time_regularized = (It + a^2 * Gttemp)

    νθ = abs(cos(θs)) < abs(cos(θo)) ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    return emission_time_regularized, rs, θs, emission_azimuth, νr, νθ, issuccess
end

function θs(metric::Kerr{T}, α, β, θo, τ) where {T}
    return _θs(metric, sign(β), θo, η(metric, α, β, θo), λ(metric, α, θo), τ)
end

function _θs(metric::Kerr{T}, signβ, θo, η, λ, τ) where {T}
    a = metric.spin
    Ghat_2, isvortical = zero(T), η < zero(T)

    Δθ = (1 - (η + λ^2) / a^2) / 2
    dsc = √(Δθ^2 + η / a^2)
    up = Δθ + dsc
    um = Δθ - dsc
    m = up / um
    k = m

    #isvortical = η < 0.
    ans, k, argo, τs, τo, τhat = zero(T), zero(T), zero(T), zero(T), zero(T), zero(T)
    tempfac = one(T) / √abs(um * a^2)
    if isvortical
        argo = clamp((cos(θo)^2 - um) / (up - um), zero(T), one(T))
        k = one(T) - m
        tempfac = inv(√abs(um * a^2))
        τo = tempfac * JacobiElliptic.F(asin(√argo), k)
        Ghat_2 = tempfac * JacobiElliptic.K(k)
        τhat = 2Ghat_2
        Δτtemp = (τ % τhat + (θo > T(π / 2) ? -one(T) : one(T)) * signβ * τo)
        n = floor(τ / τhat)
        absτs = abs(argmin(abs, (τhat - Δτtemp, Δτtemp)))
        τs = (θo > T(π / 2) ? -one(T) : one(T)) * absτs

        argr = (JacobiElliptic.sn(absτs / tempfac, k))^2
        ans = acos((θo > T(π / 2) ? -one(T) : one(T)) * √((up - um) * argr + um))
        sign = ans > T(π / 2) ? -one(T) : one(T)
        τo = sign * τo
    else
        argo = clamp(cos(θo) / √(up), -one(T), one(T))
        k = m
        tempfac = inv(√abs(um * a^2))
        τo = tempfac * JacobiElliptic.F(asin(argo), k)
        Ghat_2 = tempfac * JacobiElliptic.K(k)
        τhat = Ghat_2 + Ghat_2
        Δτtemp = (τ % τhat + signβ * τo)
        n = floor(τ / τhat)
        sign = (n % 2 == 1) ? -one(T) : one(T)
        τs = argmin(abs, (sign * signβ * (τhat - Δτtemp), sign * signβ * Δτtemp))
        newargs = JacobiElliptic.sn(τs / tempfac, k)
        ans = acos(√up * newargs)
    end

    # TODO: Come up with a more elloquent solution for this
    isincone = abs(cos(ans)) < abs(cos(θo))
    if isincone
        νθ = (n % 2 == 1) ⊻ (θo > ans)
        τ1 = (n + 1) * τhat - signβ * τo + (νθ ? 1 : -1) * τs
        isindir = isapprox(τ1, τ, rtol = sqrt(eps(T)))
    elseif ans < T(π / 2)
        τ1 = (n + 1) * τhat - signβ * τo - τs
        isindir = isapprox(τ1, τ, rtol = sqrt(eps(T)))
    else
        τ1 = (n + 1) * τhat - signβ * τo + τs
        isindir = isapprox(τ1, τ, rtol = sqrt(eps(T)))
    end

    return ans, τs, τo, τhat, n, isindir
end
