module Krang
using DocStringExtensions
using Adapt
using StaticArrays
using JacobiElliptic
using PolarizedTypes
using Roots # For Level set geometry
import GeometryBasics
using Rotations
@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(TYPEDSIGNATURES)
                                         $(DOCSTRING)
                                         """

# True iff KrangMetalExt is loaded; gate CPU-only code (e.g. `@warn`) so the rest of the function can compile to Metal.
has_metal(_...) = false

include("metrics/AbstractMetric.jl")
include("metrics/Kerr/Kerr.jl")
include("cameras/camera_types.jl")
include("cameras/SlowLightIntensityCamera.jl")
include("cameras/IntensityCamera.jl")
include("materials/AbstractMaterialTypes.jl")
include("geometries/geometry_types.jl")
include("geometries/mesh_geometry.jl")
include("geometries/level_set_geometry.jl")
include("metrics/Kerr/misc.jl")
include("metrics/Kerr/emission_coordinates.jl")
include("materials/physicsUtils.jl")
include("materials/ElectronSynchrotronPowerLawPolarization.jl")
include("materials/ElectronSynchrotronPowerLawIntensity.jl")
include("materials/ThermalSynchrotron.jl")
include("schemes/schemes.jl")
include("schemes/RayTrace.jl")

end
