"""
    $TYPEDEF

Abstract Material
"""
abstract type AbstractMaterial end

"""
    $TYPEDEF

Trait that selects which emission-coordinate call a material needs. Materials whose
emissivity does not depend on the emission azimuth `ϕs` are [`AxisymmetricTrait`](@ref)
and only need the (cheaper) emission radius; everything else is
[`NonAxisymmetricTrait`](@ref).
"""
abstract type AbstractSymmetryTrait end

"""
    $TYPEDEF

Material is axisymmetric: `ϕs` is not needed, so only the emission radius is computed.
"""
struct AxisymmetricTrait <: AbstractSymmetryTrait end

"""
    $TYPEDEF

Material depends on the emission azimuth, so related quantities are computed.
"""
struct NonAxisymmetricTrait <: AbstractSymmetryTrait end

"""
    symmetryTrait(material::AbstractMaterial)

Return the [`AbstractSymmetryTrait`](@ref) of `material`. Defaults to
[`NonAxisymmetricTrait`](@ref); axisymmetric materials should define

    Krang.symmetryTrait(::MyMaterial) = Krang.AxisymmetricTrait()
"""
symmetryTrait(::AbstractMaterial) = NonAxisymmetricTrait()

"""
    $TYPEDEF

Trait that selects whether a material needs the coordinate time of emission.
[`SlowLightTrait`](@ref) materials get the full `ts` from `emission_coordinates`;
[`FastLightTrait`](@ref) materials skip it, and `ts` is reported as zero.
"""
abstract type AbstractLightTrait end

"""
    $TYPEDEF

Material does not use the emission time, so the time integrals are skipped.
"""
struct FastLightTrait <: AbstractLightTrait end

"""
    $TYPEDEF

Material depends on the emission time, so `ts` is computed.
"""
struct SlowLightTrait <: AbstractLightTrait end

"""
    lightTrait(material::AbstractMaterial)

Return the [`AbstractLightTrait`](@ref) of `material`. Defaults to
[`SlowLightTrait`](@ref); materials that do not need the emission time should define

    Krang.lightTrait(::MyMaterial) = Krang.FastLightTrait()
"""
lightTrait(::AbstractMaterial) = SlowLightTrait()

abstract type AbstractReturnTrait end
abstract type AbstractPolarizationTrait <: AbstractReturnTrait end
struct SimplePolarizationTrait <: AbstractPolarizationTrait end
struct SimpleIntensityTrait <: AbstractReturnTrait end

function returnTrait(mat::AbstractMaterial)
    return SimpleIntensityTrait()
end
