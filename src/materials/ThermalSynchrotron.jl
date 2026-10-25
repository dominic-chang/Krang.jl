# Physical constants in SI units used by the thermal synchrotron material.
const C_SI = 2.99792458e8      # speed of light, m s⁻¹
const ME_SI = 9.1093837015e-31 # electron mass, kg
const QE_SI = 1.602176634e-19  # elementary charge, C
const KB_SI = 1.380649e-23     # Boltzmann constant, J K⁻¹
const H_SI = 6.62607015e-34    # Planck constant, J s
const KE_SI = 8.9875517923e9   # Coulomb constant, N m² C⁻²
const GAUSS_TO_TESLA = 1.0e-4

"""
Analytic approximation to the thermal synchrotron integral, eq. A18 of
[Dexter (2016)](https://doi.org/10.1093/mnras/stw1526); the fit is originally that of
[Mahadevan, Narayan & Yi (1996)](https://doi.org/10.1086/177422), who quote a maximum
error of 0.0036. It stands in for \$R(x) = \\int_x^\\infty K_{5/3}(ξ)\\,dξ\$ and holds for
`x > 1`.

```math
I(x) = 2.5651\\left(1 + 1.92x^{-1/3} + 0.9977x^{-2/3}\\right)\\exp\\left(-1.8899x^{1/3}\\right)
```

# Arguments

- `x` : Frequency in units of the critical frequency, `ν/ν_c`
"""
@inline function thermal_synchrotron_shape(x::T) where {T}
    cbrtx = cbrt(x)
    invcbrtx = inv(cbrtx)
    return T(2.5651) *
           (one(T) + T(1.92) * invcbrtx + T(0.9977) * invcbrtx * invcbrtx) *
           exp(-T(1.8899) * cbrtx)
end

"""
Planck function \$B_ν\$ in SI units (W m⁻² Hz⁻¹ sr⁻¹).

# Arguments

- `ν` : Frequency in Hz
- `temperature` : Temperature in K
"""
@inline function blackbody(ν::T, temperature) where {T}
    return 2 * T(H_SI) * ν^3 / T(C_SI)^2 / expm1(T(H_SI) * ν / (T(KB_SI) * temperature))
end

"""
    $TYPEDEF

   A synchrotron emitting disk of thermal (Maxwell-Jüttner) electrons whose density,
   temperature and magnetic field each follow a power law in radius. Emission is absorbed
   on its way out, so the material returns a [`Krang.AbsorptiveIntensity`](@ref) and
   `render` performs the radiative transfer across sub-images.

   Follows the model of
   [Desire, Cárdenas-Avendaño & Chael (2025)](https://doi.org/10.3847/1538-4357/adac4d)
   (their eqs. 13-14), with the thermal synchrotron integral replaced by the
   [`Krang.thermal_synchrotron_shape`](@ref) approximation,

```math
j_ν = \\frac{k_e n_{th} q_e^2 ν}{2\\sqrt3\\,Θ_e^2} I(ν/ν_c),
\\qquad ν_c = \\frac{3}{4π}\\frac{q_e B \\sin θ_b}{m_e}Θ_e^2,
\\qquad Θ_e = \\frac{k_B T_e}{m_e c^2},
```

   with absorption set by Kirchhoff's law, \$α_ν = j_ν/B_ν\$, and an optical depth
   \$τ = α_ν\\,ℓ/c\$. The ray crosses a shell of angular thickness \$Δθ = 1/2\$ about the
   emitting 2D surface, so \$ℓ = (r_s r_g/2)\\,\\ell_p\$ with
   \$\\ell_p = |p^{\\hat t}/p^{\\hat θ}|\$, the obliquity of the ray to the surface normal
   \$\\hat θ\$. This quantity diverges for a ray traveling parallel to the emitting surface.
   Note the proper thickness is \$\\sqrt{g_{θθ}}\\,Δθ = \\sqrt{r_s^2 + a^2\\cos^2θ_s}\\,Δθ\$ to
   first order in \$Δθ\$, so taking it to be \$r_s r_g Δθ\$ is exact only for the equatorial
   disk; the depth through a cone is underestimated by
   \$\\sqrt{1 + a^2\\cos^2θ_s/r_s^2}\$, which reaches a few per cent close to the horizon.

   The emissivity is evaluated at the observed frequency and the result is redshifted by
   \$g^3\$.

   Radii are in gravitational radii and `rg` is specified in metres. Every other
   quantity is SI except for the magnetic field, which is in Gauss.

   # Fields
    - `magnetic_field::SVector{3, T}`: Magnetic field vector `(x, y, z)` at `rb`, in Gauss.
    - `fluid_velocity::SVector{3, T}`: Fluid velocity `(speed, inclination, azimuth)`.
    - `p_dens::T`: Power-law index of the number density.
    - `p_temp::T`: Power-law index of the electron temperature.
    - `p_b::T`: Power-law index of the magnetic field magnitude.
    - `n_th0::T`: Number density at `rb`, in m⁻³.
    - `t_e0::T`: Electron temperature at `rb`, in K.
    - `rb::T`: Reference radius, in gravitational radii.
    - `ν::T`: Observed frequency, in Hz.
    - `rg::T`: Gravitational radius `GM/c²`, in metres.
    - `subimgs::NTuple{N, Int}`: The sub-images to ray trace.
"""
struct ThermalSynchrotron{N,T} <: AbstractMaterial
    magnetic_field::SVector{3,T}
    fluid_velocity::SVector{3,T}
    p_dens::T
    p_temp::T
    p_b::T
    n_th0::T
    t_e0::T
    rb::T
    ν::T
    rg::T
    subimgs::NTuple{N,Int}

    @doc """
        ThermalSynchrotron(
        magnetic_fieldx::T,
        magnetic_fieldy::T,
        magnetic_fieldz::T,
        fluid_speed::T,
        fluid_inclination_angle::T,
        fluid_azimuthal_angle::T,
        p_dens::T,
        p_temp::T,
        p_b::T,
        n_th0::T,
        t_e0::T,
        rb::T,
        ν::T,
        rg::T,
        subimgs::NTuple{N, Int},
    ) -> ThermalSynchrotron{N, T}

    Constructs a `ThermalSynchrotron` object.

    # Arguments
    - `magnetic_fieldx::T`: The x-component of the magnetic field at `rb`, in Gauss.
    - `magnetic_fieldy::T`: The y-component of the magnetic field at `rb`, in Gauss.
    - `magnetic_fieldz::T`: The z-component of the magnetic field at `rb`, in Gauss.
    - `fluid_speed::T`: The speed component of the fluid velocity vector.
    - `fluid_inclination_angle::T`: The inclination angle component of the fluid velocity vector.
    - `fluid_azimuthal_angle::T`: The azimuthal angle component of the fluid velocity vector.
    - `p_dens::T`: The power-law index of the number density.
    - `p_temp::T`: The power-law index of the electron temperature.
    - `p_b::T`: The power-law index of the magnetic field magnitude.
    - `n_th0::T`: The number density at `rb`, in m⁻³.
    - `t_e0::T`: The electron temperature at `rb`, in K.
    - `rb::T`: The reference radius, in gravitational radii.
    - `ν::T`: The observed frequency, in Hz.
    - `rg::T`: The gravitational radius `GM/c²`, in metres.
    - `subimgs::NTuple{N, Int}`: The sub-images to ray trace.

    # Returns
    - `ThermalSynchrotron{N, T}`: A new instance of `ThermalSynchrotron`.
"""
    function ThermalSynchrotron(
        magnetic_fieldx::T,
        magnetic_fieldy::T,
        magnetic_fieldz::T,
        fluid_speed::T,
        fluid_inclination_angle::T,
        fluid_azimuthal_angle::T,
        p_dens::T,
        p_temp::T,
        p_b::T,
        n_th0::T,
        t_e0::T,
        rb::T,
        ν::T,
        rg::T,
        subimgs::NTuple{N,Int},
    ) where {N,T}
        new{N,T}(
            SVector(magnetic_fieldx, magnetic_fieldy, magnetic_fieldz),
            SVector(fluid_speed, fluid_inclination_angle, fluid_azimuthal_angle),
            p_dens,
            p_temp,
            p_b,
            n_th0,
            t_e0,
            rb,
            ν,
            rg,
            subimgs,
        )
    end
end

"""
Emissivity, source function and optical depth of a [`Krang.ThermalSynchrotron`](@ref) slab
at the emission point of a ray.

Returns `(g³Bν, τ)` in SI units: the source function redshifted to the observer, and the
optical depth of the slab. The emissivity carries an extra factor of `c`, which is what
makes `τ = jν/Bν * ℓ/c` dimensionless.

# Arguments

- `metric` : Kerr metric
- `α`, `β` : Bardeen screen coordinates
- `ri`, `θs` : Emission radius and inclination
- `θo` : Observer inclination
- `magnetic_field` : Magnetic field at the emission radius, in Gauss
- `βfluid` : Fluid velocity `(speed, inclination, azimuth)`
- `νr`, `θsign` : Signs of the radial and polar momenta at emission
- `ν` : Observed frequency, in Hz
- `n_th` : Number density at the emission radius, in m⁻³
- `temperature` : Electron temperature at the emission radius, in K
- `rg` : Gravitational radius `GM/c²`, in metres
"""
@inline function thermalSynchrotronIntensity(
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
    ν,
    n_th,
    temperature,
    rg,
) where {T}
    _, p_fluid_u, _, bsinθb = synchrotron_emission_frame(
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

    # `p_fluid_u` is ordered {t, r, ϕ, θ}
    redshift = inv(p_fluid_u[1])
    lp = abs(p_fluid_u[1] / p_fluid_u[4])

    θe = T(KB_SI) * temperature / (T(ME_SI) * T(C_SI)^2)
    νc = 3 / (4 * T(π)) * T(QE_SI) * (bsinθb * T(GAUSS_TO_TESLA)) * θe^2 / T(ME_SI)
    jν =
        T(KE_SI) * n_th * T(QE_SI)^2 * ν / (2 * √T(3) * θe^2) *
        thermal_synchrotron_shape(ν / νc)
    source = blackbody(ν, temperature)

    # Kirchhoff's law for the absorption coefficient, over the path length the ray spends
    # in the cone: a shell of angular thickness 1/2, crossed at obliquity `lp`. Using
    # `ri` for the proper thickness per radian is exact for an equatorial disk and
    # understates a cone by √(1 + a²cos²θs/ri²). `jν` carries an extra factor of `c`,
    # which the `/C_SI` here removes.
    τ = (jν / source) * (ri * rg * lp / 2) / T(C_SI)

    return redshift^3 * source, τ
end

"""
    Functor for the ThermalSynchrotron material
"""
@inline function (material::ThermalSynchrotron{N,T})(
    pix::AbstractPixel,
    intersection,
) where {N,T}
    (; magnetic_field, fluid_velocity, n_th0, t_e0, rb, p_dens, p_temp, p_b) = material
    (; rs, θs, νr, νθ) = intersection

    θo = inclination(pix)
    met = metric(pix)
    α, β = screen_coordinate(pix)

    rat = rs / rb
    n_th = n_th0 * rat^p_dens
    temperature = t_e0 * rat^p_temp
    # the power law scales the magnitude and leaves the direction alone
    b_field = magnetic_field * rat^p_b

    source, τ = thermalSynchrotronIntensity(
        met,
        α,
        β,
        rs,
        θs,
        θo,
        b_field,
        fluid_velocity,
        νr,
        νθ,
        material.ν,
        n_th,
        temperature,
        material.rg,
    )

    return AbsorptiveIntensity(source * -expm1(-τ), τ)
end

lightTrait(::ThermalSynchrotron) = FastLightTrait()
symmetryTrait(::ThermalSynchrotron) = AxisymmetricTrait()
returnTrait(::ThermalSynchrotron) = AbsorptiveIntensityTrait()
