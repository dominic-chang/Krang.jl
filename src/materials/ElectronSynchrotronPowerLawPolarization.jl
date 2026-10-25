"""
Returns the screen polarization associated with a killing spinor κ as seen seen by an asymptotic observer.
"""
@inline function screen_polarization(metric::Kerr{T}, κ::Complex, θ, α, β) where {T}# Eq 31 10.1103/PhysRevD.104.044060
    a = metric.spin
    κ1 = real(κ)
    κ2 = imag(κ)

    μ = -(α + a * sin(θ))
    norm = sqrt(μ^2 + β^2)
    fα = (β * κ2 - μ * κ1) / norm
    fβ = (β * κ1 + μ * κ2) / norm

    return fα, fβ
end

function evpa(fα, fβ)
    atan(-fα, fβ)
end

"""
Calculates the polarization of a photon emitted from a fluid particle with momentum f_u and observed by an asymptotic observer.
"""
@inline function synchrotronPolarization(
    metric::Kerr{T},
    α,
    β,
    ri,
    θs,
    θo,
    magnetic_field::SVector{3,T},
    βfluid::SVector{3,T},
    νr::Bool,
    θsign::Bool,
) where {T}

    βv = βfluid[1]
    θz = βfluid[2]
    ϕz = βfluid[3]

    p_zamo_u, p_fluid_u, vec, norm = synchrotron_emission_frame(
        metric,
        α,
        β,
        ri,
        θs,
        θo,
        magnetic_field,
        βfluid,
        νr,
        θsign,
    )

    # Undo the boost and the ZAMO transform
    f_fluid_u = SVector(zero(eltype(vec)), vec[1], vec[2], vec[3])
    f_zamo_u = jac_fluid_u_zamo_d(metric, -βv, θz, ϕz, f_fluid_u)
    f_bl_u = jac_bl_u_zamo_d(metric, ri, θs, f_zamo_u)
    curr_p_bl_u = jac_bl_u_zamo_d(metric, ri, θs, p_zamo_u)

    κ1, κ2 = penrose_walker(
        metric,
        ri,
        θs,
        @SVector([curr_p_bl_u[1], curr_p_bl_u[2], curr_p_bl_u[4], curr_p_bl_u[3]]),
        @SVector([f_bl_u[1], f_bl_u[2], f_bl_u[4], f_bl_u[3]])
    )
    # normalising `κ` by `√(conj(κ)κ)` is just dividing by its modulus
    κ = complex(κ1, κ2) * inv(hypot(κ1, κ2))

    eα, eβ = screen_polarization(metric, κ, θo, α, β) .* norm

    return eα, eβ, inv(p_fluid_u[1]), abs(p_fluid_u[1] / p_fluid_u[4])
end

"""
    $TYPEDEF

   A struct representing the linear polarization of synchrotron radiation from electrons following a power-law energy distribution. This model is based on the material described in https://doi.org/10.3847/1538-4357/abf117.

   # Fields
    - `magnetic_field::SVector{3, T}`: The magnetic field vector components (x, y, z).
    - `fluid_velocity::SVector{3, T}`: The fluid velocity vector components (speed, inclination angle, azimuthal angle).
    - `spectral_index::T`: The spectral index of the electron energy distribution.
    - `R::T`: The characteristic radius of the emissivity profile.
    - `p1::T`: The first power-law index.
    - `p2::T`: The second power-law index.
    - `subimgs::NTuple{N, Int}`: The sub-images to ray trace.
"""
struct ElectronSynchrotronPowerLawPolarization{N,T} <: AbstractMaterial
    magnetic_field::SVector{3,T}
    fluid_velocity::SVector{3,T}
    spectral_index::T
    R::T
    p1::T
    p2::T
    subimgs::NTuple{N,Int}

    @doc """
        ElectronSynchrotronPowerLawPolarization(
        magnetic_fieldx::T, 
        magnetic_fieldy::T, 
        magnetic_fieldz::T, 
        fluid_speed::T, 
        fluid_inclination_angle::T, 
        fluid_azimuthal_angle::T, 
        spectral_index::T,
        R::T,
        p1::T,
        p2::T,
        subimgs::NTuple{N, Int},
    ) -> ElectronSynchrotronPowerLawPolarization{N, T}

    Constructs an `ElectronSynchrotronPowerLawPolarization` object.

    # Arguments
    - `magnetic_fieldx::T`: The x-component of the magnetic field vector.
    - `magnetic_fieldy::T`: The y-component of the magnetic field vector.
    - `magnetic_fieldz::T`: The z-component of the magnetic field vector.
    - `fluid_speed::T`: The speed component of the fluid velocity vector.
    - `fluid_inclination_angle::T`: The inclination angle component of the fluid velocity vector.
    - `fluid_azimuthal_angle::T`: The azimuthal angle component of the fluid velocity vector.
    - `spectral_index::T`: The spectral index of the electron energy distribution.
    - `R::T`: The characteristic radius of emissivity profile.
    - `p1::T`: The first power-law index.
    - `p2::T`: The second power-law index.
    - `subimgs::NTuple{N, Int}`: The sub-images to ray trace.

    # Returns
    - `ElectronSynchrotronPowerLawPolarization{N, T}`: A new instance of `ElectronSynchrotronPowerLawPolarization`.
"""
    function ElectronSynchrotronPowerLawPolarization(
        magnetic_fieldx::T,
        magnetic_fieldy::T,
        magnetic_fieldz::T,
        fluid_speed::T,
        fluid_inclination_angle::T,
        fluid_azimuthal_angle::T,
        spectral_index::T,
        R::T,
        p1::T,
        p2::T,
        subimgs::NTuple{N,Int},
    ) where {N,T}
        new{N,T}(
            SVector(magnetic_fieldx, magnetic_fieldy, magnetic_fieldz),
            SVector(fluid_speed, fluid_inclination_angle, fluid_azimuthal_angle),
            spectral_index,
            R,
            p1,
            p2,
            subimgs,
        )
    end
end

function returnTrait(::ElectronSynchrotronPowerLawPolarization{N,T}) where {N,T}
    return SimplePolarizationTrait()
end

function nan2zero(x)
    return isnan(x) ? zero(eltype(x)) : x
end

"""
    Functor for the ElectronSynchrotronPowerLawPolarization material
"""
@inline function (linpol::ElectronSynchrotronPowerLawPolarization{N,T})(
    pix::AbstractPixel,
    intersection,
) where {N,T}
    (; magnetic_field, fluid_velocity, spectral_index, R, p1, p2) = linpol
    (; rs, θs, νr, νθ) = intersection


    θo = inclination(pix)
    met = metric(pix)
    α, β = screen_coordinate(pix)

    eα, eβ, redshift, lp = synchrotronPolarization(
        met,
        α,
        β,
        rs,
        θs,
        θo,
        magnetic_field,
        fluid_velocity,
        νr,
        νθ,
    )

    rat = (rs / R)
    prof = rat^p1 / (1 + rat^(p1 + p2)) * redshift^(T(3) + spectral_index)
    q_temp = T(-(eα^2 - eβ^2) + eps(T))
    u_temp = T(-2 * eα * eβ + eps(T))
    mag = hypot(q_temp, u_temp)

    # Add a clamp to lp to help remove hot pixels
    i = mag^(one(T) + spectral_index) * min(lp, T(1e2)) * prof
    q = i * q_temp / mag
    u = i * u_temp / mag

    return StokesParams(i, q, u, zero(T))
end

lightTrait(::ElectronSynchrotronPowerLawPolarization) = FastLightTrait()
symmetryTrait(::ElectronSynchrotronPowerLawPolarization) = AxisymmetricTrait()
