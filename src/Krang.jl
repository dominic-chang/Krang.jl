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
include("Kerr/raytracer.jl")
include("Kerr/api.jl")
include("Kerr/polarization.jl")
end
