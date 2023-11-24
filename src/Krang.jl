module Krang
using DocStringExtensions
using StaticArrays
using JacobiElliptic
using UnPack
@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(TYPEDSIGNATURES)
    $(DOCSTRING)
    """

# Write your package code here.
include("metrics/AbstractMetric.jl")
include("metrics/Kerr.jl")
include("Cameras/camera_types.jl")
include("Cameras/SlowLightCamera.jl")
include("Cameras/BasicCamera.jl")
include("Kerr/raytracer.jl")
include("Kerr/api.jl")
include("Kerr/polarization.jl")
end
