module Krang
using DocStringExtensions
using StaticArrays
using JacobiElliptic
using PolarizedTypes
@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(TYPEDSIGNATURES)
    $(DOCSTRING)
    """

include("metrics/AbstractMetric.jl")
include("metrics/Kerr/Kerr.jl")
include("geometries/geometry_types.jl")
include("cameras/camera_types.jl")
include("cameras/SlowLightIntensityCamera.jl")
include("cameras/IntensityCamera.jl")
include("metrics/Kerr/misc.jl")
include("metrics/Kerr/emission_coordinates.jl")
include("metrics/Kerr/api.jl")
include("materials/physicsUtils.jl")
include("materials/material_types.jl")
include("materials/CoordinateRadius.jl")
include("materials/CoordinatePoint.jl")
include("materials/ElectronSynchrotronPowerLawPolarization.jl")
include("materials/ElectronSynchrotronPowerLawIntensity.jl")
include("schemes/schemes.jl")
include("schemes/RayTrace.jl")
include("raytracer_api.jl")
end
