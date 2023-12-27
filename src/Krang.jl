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
include("metrics/Kerr.jl")
include("Geometries/geometry_types.jl")
include("Geometries/misc.jl")
include("Cameras/camera_types.jl")
include("Cameras/SlowLightIntensityCamera.jl")
include("Cameras/IntensityCamera.jl")
include("Kerr/raytracer.jl")
include("Kerr/api.jl")
include("Materials/material_types.jl")
include("Materials/PowerLawPolarization.jl")
include("Materials/observations.jl")

end
