"""
Calculates the intensity of a photon emitted from a fluid particle with momentum f_u and observed by an asymptotic observer.
"""
function synchrotronIntensity(metric::Kerr{K}, α, β, ri, θs, θo, magfield::SVector{3,A}, βfluid::SVector{3,B}, νr::Bool, θsign::Bool) where {K, A, B}
    T = promote_type(K, A, B)
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

    return norm, inv(p_fluid_u[1]), abs(p_fluid_u[1] / p_fluid_u[4])
end
struct ElectronSynchrotronPowerLawIntensity <: AbstractMaterial end

function (prof::ElectronSynchrotronPowerLawIntensity)(pix::AbstractPixel, geometry::ConeGeometry{T,A}) where {T,A}
    magfield, fluid_velocity, subimgs, profile, σ = geometry.attributes

    θs = geometry.opening_angle
    θo = inclination(pix)
    met = metric(pix)
    α, β = screen_coordinate(pix)

    observation = zero(T)

    isindir = false
    for _ in 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        for n in subimgs
            νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
            rs, νr, _ = emission_radius(pix, geometry.opening_angle, isindir, n)
            if rs ≤ horizon(met) || isnan(rs)
                continue
            end
            norm, redshift, lp = synchrotronIntensity(met, α, β, rs, θs, θo, magfield, fluid_velocity, νr, νθ)

            prof = profile(rs) * max(redshift, eps(T))^(T(3) + σ)
            i = norm^(1 + σ) * lp * prof
            observation += nan2zero(i)
        end
    end
    return observation
end

function (linpol::ElectronSynchrotronPowerLawIntensity)(pix::AbstractPixel, geometry::UnionGeometry)
    return linpol(pix, geometry.geometry1) + linpol(pix, geometry.geometry2)
end