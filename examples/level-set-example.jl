# # Raytracing a Level set geometry
# A level set geoemtry is defined by a constraint equations $f(x,y,z)=0$.
# We will ray trace an example parabaloid geometry in this example as a simple geometric jet model.
using Krang
import CairoMakie as GLMk
using GeometryBasics
using FileIO

# Lets create a camera with a screen of 20Mx20M at a resolution of 200x200 pixels for a high spin black hole.
metric = Krang.Kerr(0.9) # Kerr metric with a spin of 0.99
θo = 89 / 180 * π # Inclination angle of the observer
ρmax = 20.0 # Horizontal and Vertical extent of the screen
sze = 200 # Resolution of the screen is sze x sze

camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze)

# We will define a parabaloid geometry to the set of all points satisfied by the equation $\cos(θ)=1-(r/r_h)^n$.
struct Parabaloid{T} <: Krang.AbstractLevelSetGeometry{T}
    rh::T
    index::T
end
function (geometry::Parabaloid)(x,y,z)
    r = sqrt(x^2+y^2+z^2)
    return 1-(r/geometry.rh)^geometry.index*(1-z/r)
end

# The jet will be emit a constant intensity whose physics we define in the `XMaterial`.
# [!NOTE] We are ignoring relativistic effects in this example. 
struct XMaterial <: Krang.AbstractMaterial end
function (mat::XMaterial)(pix, intersection) where T
    return 1.0
end

# We will ray trace the geometry and plot the image.
parabaloid = Parabaloid(Krang.horizon(metric), 1.0)

fig = GLMk.Figure();
ax = GLMk.Axis(fig[1, 1], aspect = 1)


intersections = raytrace(camera, Krang.Mesh(parabaloid, XMaterial()), res = 1_00)

# And plot the image with GLMakie, 

GLMk.heatmap!(ax, intersections, colormap = :afmhot);
fig

GLMk.save("geometric_jet.png", fig)

# ![Simple geometric jet model](geometric_jet.png)
