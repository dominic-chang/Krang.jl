"""
Returns the screen polarization associated with a killing spinor κ as seen seen by an assymptotic observer.
"""
function screen_polarisation(metric::Kerr{T}, κ::Complex, θ, α, β) where {T}# Eq 31 10.1103/PhysRevD.104.044060
    a = metric.spin
    κ1 = real(κ)
    κ2 = imag(κ)

    μ = -(α + a * sin(θ))
    norm = sqrt(μ^2 + β^2)
    fα = (β * κ2 - μ * κ1) / norm
    fβ = (β * κ1 + μ * κ2) / norm 

    return fα, fβ
end

evpa(fα, fβ) = atan(-fα, fβ)

"""
Calculates the polarization of a photon emitted from a fluid particle with momentum f_u and observed by an asymptotic observer.
"""
function synchrotronPolarization(metric::Kerr{T}, α, β, ri, θs, θo, magfield::SVector{3,T}, βfluid::SVector{3,T}, νr::Bool, θsign::Bool) where {T}


    a = metric.spin
    βv = βfluid[1]
    θz = βfluid[2]
    ϕz = βfluid[3]

    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)


    curr_p_bl_d = p_bl_d(metric, ri, θs, ηtemp, λtemp, νr, θsign)

    curr_p_bl_u = metric_uu(metric, ri, θs) * curr_p_bl_d
    p_zamo_u = jac_zamo_u_bl_d(metric, ri, θs) * curr_p_bl_u
    p_fluid_u = jac_fluid_u_zamo_d(metric, βv, θz, ϕz) * p_zamo_u
    magfieldx, magfieldy, magfieldz = magfield
    _, p_fluid_ux, p_fluid_uy, p_fluid_uz = p_fluid_u ./ p_fluid_u[1]
    vec = @SVector[p_fluid_uy * magfieldz - p_fluid_uz * magfieldy, p_fluid_uz * magfieldx - p_fluid_ux * magfieldz, p_fluid_ux * magfieldy - p_fluid_uy * magfieldx]
    norm = √sum(vec .* vec) + eps(T)
    f_fluid_u = SVector(zero(eltype(vec)), vec[1], vec[2], vec[3])
    f_zamo_u = jac_fluid_u_zamo_d(metric, -βv, θz, ϕz) * f_fluid_u
    f_bl_u = jac_bl_u_zamo_d(metric, ri, θs) * f_zamo_u
    z = zero(T)
    o = one(T)
    sinθs = sin(θs)
    A = @SMatrix [
        z o z z
        -o z z a*sinθs^2
        z z z z
        z -a*sinθs^2 z z
    ]
    B = @SMatrix [
        z z -a*sinθs z
        z z z z
        a*sinθs z z -(ri^2 + a^2)
        z z (ri^2+a^2) z
    ]
    f_temp_d = ((A - B * im) * (ri - a * cos(θs) * im)) * (f_bl_u)
    κ = sum(curr_p_bl_u .* f_temp_d)
    κ = κ * Krang._pow(conj(κ) * κ, -T(0.5))

    eα, eβ = screen_polarisation(metric, κ, θo, α, β) .* norm

    return eα, eβ, inv(p_fluid_u[1]), abs(p_fluid_u[1] / p_fluid_u[4])
end

"""
    $TYPEDEF

   Linear polarization material from https://doi.org/10.3847/1538-4357/abf117
"""
struct ElectronSynchrotronPowerLawPolarization <: AbstractMaterial end

function nan2zero(x)
    return isnan(x) ? zero(eltype(x)) : x
end

"""
    Functor for the NarayanPolarization material
"""
function (linpol::ElectronSynchrotronPowerLawPolarization)(pix::AbstractPixel, geometry::ConeGeometry{T,A}) where {T, A}
    magfield, fluid_velocity, subimgs, profile, σ = geometry.attriributes

    θs = geometry.opening_angle
    θo = inclination(pix)
    met = metric(pix)
    α, β = screen_coordinate(pix)

    observation = StokesParams(zero(T), zero(T), zero(T), zero(T))

    for n in subimgs
        for isindir in (true, false)
            νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
            rs, νr, _ =  emission_radius(pix, geometry.opening_angle, isindir, n)
            eα, eβ, redshift, lp = synchrotronPolarization(met, α, β, rs, θs, θo, magfield, fluid_velocity, νr, νθ)

            prof = profile(rs)*max(redshift , eps(T))^(T(3)+σ)
            q = T(-(eα^2 - eβ^2) + eps(T))
            u = T(-2*eα*eβ + eps(T))
            i = hypot(q, u)^(1+σ)*lp*prof
            observation += StokesParams(nan2zero(i), nan2zero(q), nan2zero(u), zero(T))
        end
    end
    return observation
end

function (linpol::ElectronSynchrotronPowerLawPolarization)(pix::AbstractPixel, geometry::UnionGeometry)
    #return @SVector[0f0, 0f0, 0f0, 0f0]
    return linpol(pix, geometry.geometry1) + linpol(pix, geometry.geometry2)
end