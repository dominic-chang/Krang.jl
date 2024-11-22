"""
    $TYPEDEF

Abstract Material
"""
abstract type AbstractMaterial end

#TODO: replace with trait
isOcclusive(material::AbstractMaterial) = false
#TODO: replace with trait
isAxisymmetric(material::AbstractMaterial) = false
#TODO: replace with trait
isFastLight(material::AbstractMaterial) = false

yield(material::AbstractMaterial) = 0

abstract type AbstractReturnTrait end
abstract type AbstractPolarizationTrait <: AbstractReturnTrait end
struct SimplePolarizationTrait <: AbstractPolarizationTrait end
struct SimpleIntensityTrait <: AbstractReturnTrait end 

function returnTrait(mat::AbstractMaterial)
    return SimpleIntensityTrait() 
end
