# metrics/Kerr/emission_coordinates.jl
# Reactant specialisations mirroring the corresponding file in `src/`.

"""
Branch-free three-way select over the radial-root cases, shared by both `emission_radius`
entry points. All three `_rs_case*` helpers are evaluated and `ifelse`-selected on the real
root count rather than dispatched through a `@trace if`: each of them is well defined for
any root layout (case 1&2 takes real parts, case 3 works in the complex plane, case 4 reads
imaginary parts), and each already masks its own failure to the documented zero sentinel, so
the values the select discards can be garbage but never a throw.
"""
@inline function _reactant_rs_cases(pix::Krang.AbstractPixel, rh, τ)
    numreals = sum(x -> Base.ifelse(_isreal2(x), 1, 0), roots(pix))

    rs4, νr4, success4 = _rs_case1_and_2(pix, rh, τ)
    rs2, νr2, success2 = _rs_case3(pix, rh, τ)
    rs0, νr0, success0 = _rs_case4(pix, rh, τ)

    is4 = numreals == 4
    is2 = numreals == 2
    rs = Base.ifelse(is4, rs4, Base.ifelse(is2, rs2, rs0))
    νr = Base.ifelse(is4, νr4, Base.ifelse(is2, νr2, νr0))
    issuccess = Base.ifelse(is4, success4, Base.ifelse(is2, success2, success0))

    return rs, νr, numreals, issuccess
end

Reactant.@reactant_overlay function Krang.emission_radius(pix::Krang.AbstractPixel, τ)
    met = metric(pix)
    a = met.spin
    T = typeof(a)
    rh = one(T) + √(one(T) - a^2)

    return _reactant_rs_cases(pix, rh, τ)
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
