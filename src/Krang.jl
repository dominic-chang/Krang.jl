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
include("Cameras/BasicGPUCamera.jl")
include("Kerr/raytracer.jl")
include("Kerr/api.jl")
include("Observables/observable_types.jl")
include("Observables/PowerLawPolarization.jl")
#include("Observables/LensingBand.jl")
include("Observables/observations.jl")

end
