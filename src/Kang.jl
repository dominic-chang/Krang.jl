module Kang
using DocStringExtensions
using StaticArrays
using FastElliptic
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

end
