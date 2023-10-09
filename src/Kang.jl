module Kang
using DocStringExtensions
using StaticArrays
using FastElliptic

# Write your package code here.
include("metrics/AbstractMetric.jl")
include("metrics/Kerr.jl")
include("Kerr/raytracer.jl")
include("Kerr/api.jl")

end
