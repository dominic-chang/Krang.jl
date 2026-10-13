module KrangReactantExt
#include("SlowLightIntensityPixel.jl")

using Krang
import Krang.JacobiElliptic as JacobiElliptic
using Reactant
import Reactant: TracedRNumber


types = (TracedRNumber, Any)
Ts = Union{TracedRNumber,Any}
_isreal2 = Krang._isreal2
Ir_inf = Krang.Ir_inf
horizon = Krang.horizon
Ir_s = Krang.Ir_s
Iϕ_inf_case2 = Krang.Iϕ_inf_case2
Iϕ_inf_case3 = Krang.Iϕ_inf_case3
Iϕ_inf_case4 = Krang.Iϕ_inf_case4
It_inf_case2 = Krang.It_inf_case2
It_inf_case3 = Krang.It_inf_case3
It_inf_case4 = Krang.It_inf_case4
radial_inf_integrals_case2 = Krang.radial_inf_integrals_case2
radial_inf_integrals_case3 = Krang.radial_inf_integrals_case3
radial_inf_integrals_case4 = Krang.radial_inf_integrals_case4
radial_w_I0_terms_integrals_case2 = Krang.radial_w_I0_terms_integrals_case2
radial_w_I0_terms_integrals_case3 = Krang.radial_w_I0_terms_integrals_case3
radial_w_I0_terms_integrals_case4 = Krang.radial_w_I0_terms_integrals_case4
_get_root_diffs = Krang._get_root_diffs
screen_coordinate = Krang.screen_coordinate
metric = Krang.metric
inclination = Krang.inclination
roots = Krang.roots
I0_inf = Krang.I0_inf
absGθo_Gθhat = Krang.absGθo_Gθhat
αboundary = Krang.αboundary
βboundary = Krang.βboundary
Intersection = Krang.Intersection
R1 = Krang.R1
R2 = Krang.R2
S1 = Krang.S1
S2 = Krang.S2
regularized_Pi = Krang.regularized_Pi
regularized_R1 = Krang.regularized_R1
regularized_R2 = Krang.regularized_R2
regularizedS1 = Krang.regularizedS1
regularizedS2 = Krang.regularizedS2
ellF = JacobiElliptic.F
ellK = JacobiElliptic.K
ellPi = JacobiElliptic.Pi
ellsn = JacobiElliptic.sn
ellcn = JacobiElliptic.cn
ellsc = JacobiElliptic.sc
_rs_case1_and_2 = Krang._rs_case1_and_2
_rs_case3 = Krang._rs_case3
_rs_case4 = Krang._rs_case4

@inline _ellE(m) = JacobiElliptic.E(m)
@inline _ellE(φ, m) = JacobiElliptic.E(φ, m)

include("utils.jl")
include("metrics/Kerr/misc.jl")
include("metrics/Kerr/emission_coordinates.jl")
include("cameras/IntensityCamera.jl")
include("cameras/SlowLightIntensityCamera.jl")
include("geometries/geometry_types.jl")

end
