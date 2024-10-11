"""
Calculates the intensity of a photon emitted from a fluid particle with momentum f_u and observed by an asymptotic observer.
"""
function synchrotronIntensity(metric::Kerr{T}, α, β, ri, θs, θo, magnetic_field::SVector{3,T}, βfluid::SVector{3,T}, νr::Bool, θsign::Bool) where {T}

    βv = βfluid[1]
    θz = βfluid[2]
    ϕz = βfluid[3]

    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)

    curr_p_bl_d = p_bl_d(metric, ri, θs, ηtemp, λtemp, νr, θsign)

    curr_p_bl_u = metric_uu(metric, ri, θs) * curr_p_bl_d
    p_zamo_u = jac_zamo_u_bl_d(metric, ri, θs) * curr_p_bl_u
    p_fluid_u = jac_fluid_u_zamo_d(metric, βv, θz, ϕz) * p_zamo_u
    magnetic_fieldx, magnetic_fieldy, magnetic_fieldz = magnetic_field
    _, p_fluid_ux, p_fluid_uy, p_fluid_uz = p_fluid_u ./ p_fluid_u[1]
    vec = @SVector[p_fluid_uy * magnetic_fieldz - p_fluid_uz * magnetic_fieldy, p_fluid_uz * magnetic_fieldx - p_fluid_ux * magnetic_fieldz, p_fluid_ux * magnetic_fieldy - p_fluid_uy * magnetic_fieldx]
    norm = √sum(vec .* vec) + eps(T)

    return norm, inv(p_fluid_u[1]), abs(p_fluid_u[1] / p_fluid_u[4])
end


"""
    $TYPEDEF

   A struct representing the linear polarization intensity of synchrotron radiation from electrons following a power-law energy distribution. This model is based on the material described in https://doi.org/10.3847/1538-4357/abf117.

   # Fields
    - `magnetic_field::SVector{3, T}`: The magnetic field vector components (x, y, z).
    - `fluid_velocity::SVector{3, T}`: The fluid velocity vector components (speed, inclination angle, azimuthal angle).
    - `spectral_index::T`: The spectral index of the electron energy distribution.
    - `R::T`: The characteristic radius of the emissivity profile.
    - `p1::T`: The first power-law index.
    - `p2::T`: The second power-law index.
    - `subimgs::NTuple{N, Int}`: The sub-images to raytrace.
"""
struct ElectronSynchrotronPowerLawIntensity{N, T} <: AbstractMaterial 
    magnetic_field::SVector{3, T}
    fluid_velocity::SVector{3, T}
    spectral_index::T
    R::T
    p1::T
    p2::T
    subimgs::NTuple{N, Int}

    @doc """
        ElectronSynchrotronPowerLawIntensity(
        magnetic_fieldx::T, 
        magnetic_fieldy::T, 
        magnetic_fieldz::T, 
        fluid_speed::T, 
        fluid_inclination_angle::T, 
        fluid_azimuthal_angle::T, 
        spectral_index::T
        R::T,
        p1::T,
        p2::T,
        subimgs::NTuple{N, Int},
    ) -> ElectronSynchrotronPowerLawIntensity{N, T}

    Constructs an `ElectronSynchrotronPowerLawIntensity` object.

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
    - `subimgs::NTuple{N, Int}`: The sub-images to raytrace.

    # Returns
    - `ElectronSynchrotronPowerLawIntensity{N, T}`: A new instance of `ElectronSynchrotronPowerLawPolarization`.
"""   
    function ElectronSynchrotronPowerLawIntensity(
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
    ) where {N,T}
        new{N, T}(SVector(magnetic_fieldx, magnetic_fieldy, magnetic_fieldz), SVector(fluid_speed, fluid_inclination_angle, fluid_azimuthal_angle), spectral_index, R, p1, p2, subimgs)
    end
end

"""
    Functor for the the ElectronSynchrotronPowerLawIntensity material
"""
function (linpol::ElectronSynchrotronPowerLawIntensity{N,T})(pix::AbstractPixel, geometry::ConeGeometry{T}) where {N,T}
    (;magnetic_field, fluid_velocity, subimgs, R, p1, p2, spectral_index) = linpol

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
            if !(horizon(met) ≤ rs < T(Inf))
                continue
            end
            norm, redshift, lp = synchrotronIntensity(met, α, β, rs, θs, θo, magnetic_field, fluid_velocity, νr, νθ)

            rat = (rs/R)
            prof = rat^p1/(1+rat^(p1+p2)) * max(redshift, eps(T))^(T(3) + spectral_index)
            i = norm^(1 + spectral_index) * lp * prof
            observation += isnan(i) ? zero(T) : i
            #observation += nan2zero(rs)
        end
    end
    return observation
end
