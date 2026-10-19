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
@inline function emission_radius(pix::Krang.AbstractPixel, θs::T, isindir, n) where {T}
    α, β = screen_coordinate(pix)
    θo = inclination(pix)
    met = @inline metric(pix)
    isincone = θo ≤ θs ≤ (π - θo) || (π - θo) ≤ θs ≤ θo
    err_return = (zero(T), true, true, 0, false)
    if !isincone
        αmin = @inline αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? @inline(βboundary(met, α, θo, θs)) : zero(T))
        ((abs(β) + eps(T)) < βbound) && return err_return
    end

    τ, _, _, _, _, issuccess = Gθ(pix, θs, isindir, n)
    issuccess || return err_return

    # νθ is θ̇s increasing or decreasing?
    νθ = isincone ? (θo > θs) ⊻ (n % 2 == 1) : !isindir

    ans, νr, numreals, issuccess = @inline emission_radius(pix, τ)
    issuccess || return err_return
    return ans, νr, νθ, numreals, issuccess
end

"""
Emission radius for point originating at at Mino time τ whose image appears at the screen coordinate (`α`, `β`). 
Returns 0 if the emission coordinates do not exist for that screen coordinate.

# Arguments

-`pix` : Pixel information
- `τ` : Mino time
"""
@inline function emission_radius(pix::AbstractPixel, τ::T)::Tuple{T,Bool,Int,Bool} where {T}
    met = metric(pix)
    a = met.spin

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
Emission inclination for a point originating at radius `rs` whose image appears at screen coordinate (`α`, `β`).

# Arguments

- `pix` : Pixel information
- `rs` : Emission radius
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
"""
@inline function emission_inclination(pix::AbstractPixel, rs, νr)
    return emission_inclination(pix, Ir(pix, νr, rs))
end

"""
Emission inclination for point at Mino time τ whose image appears at screen coordinate (`α`, `β`).

# Arguments

- `pix` : Pixel information
- `τ` : Mino Time
"""
@inline function emission_inclination(pix::AbstractPixel, τ)
    met = metric(pix)
    θo = inclination(pix)
    _, β = pix.screen_coordinate

    # `absGθo_Gθhat(pix)` is exactly the `(τo, τhat)` `_θs` would recompute, so hand the
    # cached pair over instead of paying for another `F` and `K`.
    Go, Ghat = absGθo_Gθhat(pix)
    return _θs(met, sign(β), θo, η(pix), λ(pix), τ, Go, Ghat)
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
@inline function emission_azimuth(
    pix::AbstractPixel,
    θs,
    rs,
    τ::T,
    νr,
    isindir,
    n,
) where {T}
    Gϕtemp, _, _, _ = @inline Gϕ(pix, θs, isindir, n)
    return @inline emission_azimuth(pix, rs, τ, νr, Gϕtemp)
end

"""
Emission azimuth from an already-computed `Gϕ`.

For callers that already hold `Gϕ` and would otherwise pay for a second pass through
[`Gϕ`](@ref), which redoes the whole angular setup.

# Arguments

- `pix` : Pixel information
- `rs` : Emission radius
- `τ` : Mino Time
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
- `Gϕtemp` : The ϕ angular antiderivative for this emission inclination
"""
@inline function emission_azimuth(pix::AbstractPixel, rs, τ, νr, Gϕtemp)
    Iϕ = Krang.Iϕ(pix, rs, τ, νr)
    (isnan(Iϕ) || !isfinite(Iϕ)) && return Iϕ
    (isnan(Gϕtemp) || !isfinite(Gϕtemp)) && return Iϕ

    return (Iϕ + λ(pix) * Gϕtemp)
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
@inline function emission_coordinates_fast_light(
    pix::AbstractPixel,
    θs::T,
    isindir,
    n,
) where {T}
    # NB: I do not return θs since it is already known, and returning it might encourage bad coding errors
    α, β = screen_coordinate(pix)
    θo = inclination(pix)
    met = metric(pix)
    isincone = θo ≤ θs ≤ (π - θo) || (π - θo) ≤ θs ≤ θo
    err_return = (zero(T), zero(T), false, false, false)
    if !isincone
        αmin = αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? βboundary(met, α, θo, θs) : zero(T))
        ((abs(β) + eps(T)) < βbound) && return err_return
    end

    τ, _, _, _, _, issuccess = @inline Gθ(pix, θs, isindir, n)
    issuccess || return err_return

    # is θ̇s increasing or decreasing?
    νθ = isincone ? (θo > θs) ⊻ (n % 2 == 1) : !isindir

    rs, νr, _, issuccess = emission_radius(pix, τ)
    issuccess || return err_return

    # `Gϕ` is deliberately left until after the `emission_radius` gate: on a pixel that does
    # not cache `absGϕo_Gϕhat` it is the single most expensive term here, and hoisting it up
    # to share the angular setup with `Gθ` costs far more on rays that never get this far.
    ϕs = @inline emission_azimuth(pix, θs, rs, τ, νr, isindir, n)

    return rs, ϕs, νr, νθ, issuccess
end

"""
Ray trace a point that appears at the screen coordinate (`α`, `β`) for an observer located at inclination θo

# Arguments

- `pix` : Pixel information
- `τ` : Mino Time
"""
@inline function emission_coordinates_fast_light(pix::AbstractPixel, τ::T) where {T}
    α, β = @inline screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)

    θs, _, _, _, n, isindir = @inline emission_inclination(pix, τ)
    cosθs = cos(θs)
    cosθo = cos(θo)

    err_return = (zero(T), zero(T), zero(T), true, true, false)
    if cosθs > abs(cosθo)
        αmin = @inline αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? @inline(βboundary(met, α, θo, θs)) : zero(T))

        if (abs(β) + eps(T)) < βbound
            return err_return
        end
    end

    a = met.spin

    λtemp = λ(pix)
    rs, νr, _, issuccess = @inline emission_radius(pix, τ)
    issuccess || return err_return
    _, _, _, Ip, Im = @inline radial_integrals(pix, rs, τ, νr)

    rtΔ = √(one(T) - a^2)
    rp = one(T) + rtΔ
    rm = one(T) - rtΔ

    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, _, _, _, _ = @inline Gϕ(pix, θs, isindir, n)

    emission_azimuth = (Iϕ + λtemp * Gϕtemp)

    νθ = abs(cosθs) < abs(cosθo) ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
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
@inline function emission_coordinates(pix::AbstractPixel, θs::T, isindir, n) where {T}
    α, β = pix.screen_coordinate
    met = metric(pix)
    θo = inclination(pix)
    cosθs = cos(θs)
    cosθo = cos(θo)
    isincone = abs(cosθs) < abs(cosθo)
    err_return = (zero(T), zero(T), zero(T), false, false, false)
    if cosθs > abs(cosθo)
        αmin = @inline αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? @inline(βboundary(met, α, θo, θs)) : zero(T))

        if (abs(β) + eps(T)) < βbound
            return err_return
        end
    end

    τ, _, _, _, _, issuccess = @inline Gθ(pix, θs, isindir, n)
    issuccess || return err_return

    rs, νr, _, issuccess = @inline emission_radius(pix, τ)
    issuccess || return err_return

    a = met.spin

    λtemp = λ(pix)

    # `isindir` is only overwritten on the `isincone` branch, which is exactly the branch
    # whose `νθ` does not read it, so this value still holds after the reassignment below.
    νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    if isincone
        isindir = ((sign(β) > 0) ⊻ (θo > T(π / 2)))
    end
    if isnan(τ)
        return err_return
    end

    I0, I1, I2, Ip, Im = @inline radial_integrals(pix, rs, τ, νr)

    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    It =
        4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) +
        4 * I0 +
        2 * I1 +
        I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    # One angular setup for both; bit-identical to calling `Gϕ` and `Gt` separately.
    Gϕtemp, Gttemp, _ = @inline _Gϕt(pix, θs, isindir, n)

    emission_azimuth = (Iϕ + λtemp * Gϕtemp)
    emission_time_regularized = (zero(T) + It + a^2 * Gttemp)

    return emission_time_regularized, rs, emission_azimuth, νr, νθ, issuccess
end

"""
Ray trace a point that appears at the screen coordinate (`α`, `β`) for an observer located at inclination θo

# Arguments

- `pix` : Pixel information
- `τ` : Mino Time
"""
@inline function emission_coordinates(pix::AbstractPixel, τ::T) where {T}
    α, β = screen_coordinate(pix)
    met = metric(pix)
    θo = inclination(pix)

    θs, _, _, _, n, isindir = @inline emission_inclination(pix, τ)
    err_return = (zero(T), zero(T), zero(T), zero(T), true, true, false)
    if cos(θs) > abs(cos(θo))
        αmin = @inline αboundary(met, θs)
        βbound = (abs(α) >= (αmin + eps(T)) ? @inline(βboundary(met, α, θo, θs)) : zero(T))

        if (abs(β) + eps(T)) < βbound
            return err_return
        end
    end

    a = met.spin

    λtemp = λ(pix)
    rs, νr, _, issuccess = @inline emission_radius(pix, τ)
    issuccess || return err_return
    I0, I1, I2, Ip, Im = @inline radial_integrals(pix, rs, τ, νr)

    rp = one(T) + √(one(T) - a^2)
    rm = one(T) - √(one(T) - a^2)

    It =
        4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) +
        4 * I0 +
        2 * I1 +
        I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    # One angular setup for both; bit-identical to calling `Gϕ` and `Gt` separately.
    Gϕtemp, Gttemp, _ = @inline _Gϕt(pix, θs, isindir, n)

    emission_azimuth = (Iϕ + λtemp * Gϕtemp)
    emission_time_regularized = (It + a^2 * Gttemp)

    νθ = abs(cos(θs)) < abs(cos(θo)) ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > T(π / 2))
    return emission_time_regularized, rs, θs, emission_azimuth, νr, νθ, issuccess
end

function θs(metric::Kerr{T}, α, β, θo, τ) where {T}
    return _θs(metric, sign(β), θo, η(metric, α, β, θo), λ(metric, α, θo), τ)
end

@inline _θs(metric::Kerr{T}, signβ, θo, η, λ, τ) where {T} =
    _θs(metric, signβ, θo, η, λ, τ, _absGθo_Gθhat(metric, θo, η, λ)...)

@doc """
    _θs(metric, signβ, θo, η, λ, τ, Go, Ghat)

Emission inclination at Mino time `τ`, with the angular antiderivative pair
`(Go, Ghat) = `[`_absGθo_Gθhat`](@ref)` supplied by the caller.
"""
@inline function _θs(metric::Kerr{T}, signβ, θo, η, λ, τ, Go, Ghat) where {T}
    a = metric.spin
    isvortical = η < zero(T)

    Δθ = (1 - (η + λ^2) / a^2) / 2
    # `_absGθo_Gθhat` hands back `(0, 0)` for rays whose angular turning points are
    # degenerate (`a → 0`, or `argo` outside its branch). Reusing that would put a zero in
    # `τ / τhat`, so those fall back to computing the pair here, with the unguarded turning
    # points the cache-free path has always used.
    usable = isfinite(Ghat) & (Ghat > zero(T))
    dsc = usable ? √max(Δθ^2 + η / a^2, zero(T)) : √(Δθ^2 + η / a^2)
    up = usable ? min(Δθ + dsc, one(T) - eps(T)) : Δθ + dsc
    um = Δθ - dsc
    m = up / um

    τo = Go
    τhat = Ghat
    if isvortical
        k = one(T) - m
        tempfac = inv(√abs(um * a^2))
        if !usable
            τo =
                tempfac * JacobiElliptic.F(
                    asin(√clamp((cos(θo)^2 - um) / (up - um), zero(T), one(T))),
                    k,
                )
            τhat = 2 * tempfac * JacobiElliptic.K(k)
        end
        Δτtemp = (τ % τhat + (θo > T(π / 2) ? -one(T) : one(T)) * signβ * τo)
        n = unsafe_trunc(Int, τ / τhat)
        absτs = abs(argmin(abs, (τhat - Δτtemp, Δτtemp)))
        τs = (θo > T(π / 2) ? -one(T) : one(T)) * absτs

        argr = (JacobiElliptic.sn(absτs / tempfac, k))^2
        ans = acos((θo > T(π / 2) ? -one(T) : one(T)) * √((up - um) * argr + um))
        sign = ans > T(π / 2) ? -one(T) : one(T)
        τo = sign * τo
    else
        k = m
        tempfac = inv(√abs(um * a^2))
        if !usable
            τo =
                tempfac * JacobiElliptic.F(asin(clamp(cos(θo) / √(up), -one(T), one(T))), k)
            τhat = 2 * tempfac * JacobiElliptic.K(k)
        end
        Δτtemp = (τ % τhat + signβ * τo)
        n = unsafe_trunc(Int, τ / τhat)
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
