module Krang
using DocStringExtensions
using StaticArrays
using JacobiElliptic
@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(TYPEDSIGNATURES)
    $(DOCSTRING)
    """

# Write your package code here.
include("metrics/AbstractMetric.jl")
include("metrics/Kerr/Kerr.jl")
include("geometries/geometry_types.jl")
include("cameras/camera_types.jl")
include("cameras/SlowLightIntensityCamera.jl")
include("cameras/IntensityCamera.jl")
include("metrics/Kerr/raytracer.jl")
include("metrics/Kerr/api.jl")
include("materials/material_types.jl")
include("materials/PowerLawPolarization.jl")
include("materials/observations.jl")

end
