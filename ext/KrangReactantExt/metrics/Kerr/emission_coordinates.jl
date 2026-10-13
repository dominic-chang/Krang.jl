# metrics/Kerr/emission_coordinates.jl
# Reactant specialisations mirroring the corresponding file in `src/`.

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
        Δτtemp =
            mod(τ, τhat) + Base.ifelse(θo > T(Base.pi / 2), -one(T), one(T)) * signβ * τo
        n = unsafe_trunc(Int, τ / τhat)
        cand1 = τhat - Δτtemp
        cand2 = Δτtemp
        absτs = abs(Base.ifelse(abs(cand1) <= abs(cand2), cand1, cand2))
        τs = Base.ifelse(θo > T(Base.pi / 2), -one(T), one(T)) * absτs
        argr = ellsn(absτs / tempfac, k)^2
        ans = acos(
            Base.ifelse(θo > T(Base.pi / 2), -one(T), one(T)) * √((up - um) * argr + um),
        )
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

function _reactant_emission_radius_tau(pix::Krang.SlowLightIntensityPixel, τ)
    rh = horizon(metric(pix))
    numreals = sum(x->Base.ifelse(_isreal2(x), 1, 0), roots(pix))
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

function _reactant_emission_radius_theta(pix::Krang.SlowLightIntensityPixel, θs, isindir, n)
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

Reactant.@reactant_overlay function Krang.emission_radius(pix::Krang.AbstractPixel, τ)
    met = metric(pix)
    a = met.spin
    T = typeof(a)
    ans = zero(T)
    νr = true

    rh = one(T) + √(one(T) - a^2)

    numreals = sum(x->Base.ifelse(_isreal2(x), 1, 0), roots(pix))

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
    It =
        4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) +
        4 * I0 +
        2 * I1 +
        I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, Gϕvalid = _reactant_Gϕ_minotime(pix, θs, isindir, n)
    Gttemp, Gtvalid = _reactant_Gt_minotime(pix, θs, isindir, n)

    emission_azimuth = Iϕ + λtemp * Gϕtemp
    emission_time_regularized = It + a^2 * Gttemp

    νθ = Base.ifelse(
        abs(cosθs) < abs(cosθo),
        (n % 2 == 1) ⊻ (θo > θs),
        !isindir ⊻ (θs > TT(Base.pi / 2)),
    )
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

    It =
        4 / (rp - rm) * (rp * (rp - a * λtemp / 2) * Ip - rm * (rm - a * λtemp / 2) * Im) +
        4 * I0 +
        2 * I1 +
        I2
    Iϕ = 2a / (rp - rm) * ((rp - a * λtemp / 2) * Ip - (rm - a * λtemp / 2) * Im)

    Gϕtemp, Gϕvalid = _reactant_Gϕ_minotime(pix, θs, isindir, n)
    Gttemp, Gtvalid = _reactant_Gt_minotime(pix, θs, isindir, n)

    emission_azimuth = Iϕ + λtemp * Gϕtemp
    emission_time_regularized = It + a^2 * Gttemp
    νθ = Base.ifelse(
        abs(cos(θs)) < abs(cos(θo)),
        (n % 2 == 1) ⊻ (θo > θs),
        !isindir ⊻ (θs > TT(Base.pi / 2)),
    )
    issuccess = (!blocked) & rssuccess & Gϕvalid & Gtvalid

    t = Base.ifelse(issuccess, emission_time_regularized, zero(TT))
    r = Base.ifelse(issuccess, rs, zero(TT))
    θ = Base.ifelse(issuccess, θs, zero(TT))
    ϕ = Base.ifelse(issuccess, emission_azimuth, zero(TT))
    νr = Base.ifelse(issuccess, νr, true)
    νθ = Base.ifelse(issuccess, νθ, true)

    return t, r, θ, ϕ, νr, νθ, issuccess
end
