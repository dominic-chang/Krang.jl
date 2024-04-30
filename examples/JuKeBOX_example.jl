# # Creating a Custom Dual Cone Model

# In this example, we will raytrace the region around a Kerr blackhole as seen by an observer stationed at infinity.
# We will show the emission coordinates of the n=0 (direct) and n=1 (indirect) photons as they are emitted from the 
# source, at a fixed inclination angle from the blackhole's spin axis.
#
# First, let's import Krang and CairoMakie for plotting.
#using Revise
using Krang
using WGLMakie 
using Metal

curr_theme = Theme(
    Axis = (
        xticksvisible = false, 
        xticklabelsvisible = false,
        xgridvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        ygridvisible = false,
        gridvisible = false,
        topspinevisible = false,
        bottomspinevisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
        ),
    Text = (
        align = (:left,:center),
        ),
)
set_theme!(merge!(curr_theme, theme_latexfonts()))

#
# We will use a 0.99 spin Kerr blackhole viewed by an assymptotic observer at an inclination angle of θo=π/4. 
# A region spanned by radii between the horizon and 20M at varying inclinations will be raytraced onto the 20Mx20M 
# screen of the observer.
sze = 150;
metric = Krang.Kerr(-0.94);
rmin = Krang.horizon(metric)
ρmax = 10f0;
η1 = 2.6444786738735804f0

# Create the material

fig = Figure(resolution = (500, 1000));

ax = Axis(fig[1, 1:2], aspect=1)
textTag(loc, text) = begin
    ax = Axis(loc, height = 30, width=160)
    xlims!(ax,0,1)
    WGLMakie.text!(ax, 0,0; text=text)
end


ivals = Matrix{Float32}(undef, sze, sze)
sl_θ = WGLMakie.Makie.Slider(fig[2, 1], range = 0f0:0.1f0:1f0π, startvalue = 17f0/180*π)
textTag(fig[2, 2], "Spin Axis Inclination")
sl_θs = WGLMakie.Makie.Slider(fig[3, 1], range = 0f0:0.1f0:1f0π/2, startvalue = 40f0/180*π)
textTag(fig[3, 2], "Jet Opening Angle")
sl_a = WGLMakie.Makie.Slider(fig[4, 1], range = 0.001f0:0.1f0:1f0, startvalue = 0.5f0)
textTag(fig[4, 2], "Spin")
sl_β = WGLMakie.Makie.Slider(fig[5, 1], range = 0f0:0.1f0:0.999f0, startvalue = 0.5f0)
textTag(fig[5, 2], "Fluid Speed")
sl_χ = WGLMakie.Makie.Slider(fig[6, 1], range = -1f0π:0.1f0:1f0π, startvalue = 0f0)
textTag(fig[6, 2], "Fluid Direction")
sl_R = WGLMakie.Makie.Slider(fig[7, 1], range = 1f0:0.1f0:10f0, startvalue = 3f0)
textTag(fig[7, 2], "Number Density\nCharateristic Radius")
sl_p1 = WGLMakie.Makie.Slider(fig[8, 1], range = 1f0:0.1f0:10f0, startvalue = 3f0)
textTag(fig[8, 2], "Number Density Inner\nRadial Exponent")
sl_p2 = WGLMakie.Makie.Slider(fig[9, 1], range = 1f0:0.1f0:10f0, startvalue = 3f0)
textTag(fig[9, 2], "Number Density Outer\nRadial Exponent")
sl_ι = WGLMakie.Makie.Slider(fig[10, 1], range = 0f0:0.1f0:1f0π, startvalue = π/2f0)
textTag(fig[10, 2], "Magnetic Field Inclination Angle")
sl_η = WGLMakie.Makie.Slider(fig[11, 1], range = -1f0π:0.1f0:1f0π, startvalue = π/2f0)
textTag(fig[11, 2], "Magnetic Field Azimuthal Angle")
sl_σ = WGLMakie.Makie.Slider(fig[12, 1], range = -1f0:0.1f0:10f0, startvalue = 0f0)
textTag(fig[12, 2], "Energy Spectral Index")

lift(
    sl_θ.value,
    sl_θs.value,
    sl_a.value,
    sl_β.value,
    sl_χ.value,
    sl_R.value,
    sl_p1.value,
    sl_p2.value,
    sl_ι.value,
    sl_η.value,
    sl_σ.value
) do θo, θs, a, βv, χ, R, p1, p2, ι, η1, σ
metric = Krang.Kerr(a);
camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
η2 = 1f0π-η1
magfield1 = Krang.SVector(sin(ι)*cos(η1), sin(ι)*sin(η1), cos(ι));
magfield2 = Krang.SVector(sin(ι)*cos(η2), sin(ι)*sin(η2), cos(ι));
vel = Krang.SVector(βv, (π/2f0), χ);
material = Krang.ElectronSynchrotronPowerLawIntensity();

profile(r::Float32) = let R=R, p1=p1, p2=p2
    return (r/R)^p1/(1f0+(r/R)^(p1+p2))
end

geometry1 = Krang.ConeGeometry((θs), (magfield1, vel, (0,1,2), profile, σ))
geometry2 = Krang.ConeGeometry((1f0π-θs), (magfield2, vel, (0,1,2), profile, σ))
geometry = geometry1 ⊕ geometry2

mesh = Krang.Mesh(geometry, material)
ivals .= Array(mesh.material.(MtlArray(camera.screen.pixels), Ref(mesh.geometry)))


#@Threads.threads for I in CartesianIndices(camera.screen.pixels)
#    ivals[I] = mesh.material(camera.screen.pixels[I], mesh.geometry)
#end

heatmap!(ax, ivals)
end



fig