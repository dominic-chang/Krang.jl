# # Raytracing with inclination

# In this example, we will raytrace the region around a Kerr blackhole as seen by an observer stationed at infinity.
# We will show the emission coordinates of the n=0 (direct) and n=1 (indirect) photons as they are emitted from the 
# source, at a fixed inclination angle from the blackhole's spin axis.
#
# First, let's import Krang and CairoMakie for plotting.
using CairoMakie
using Krang
using ImageFiltering

curr_theme = Theme(
    Axis = (
        xticksvisible = false, 
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        ),
)
set_theme!(merge!(curr_theme, theme_latexfonts()))
figure = Figure(resolution=(600, 600), fontfamily="Computer Modern", fontface="bold");
ax = Axis(figure[1, 1], aspect=1)

#
# We will use a 0.99 spin Kerr blackhole viewed by an assymptotic observer at an inclination angle of θo=π/4. 
# A region spanned by radii between the horizon and 20M at varying inclinations will be raytraced onto the 20Mx20M 
# screen of the observer.
metric = Krang.Kerr(0.99);
θo = 20 * π / 180;
θs = 45 * π / 180
sze = 1000;
rmin = Krang.horizon(metric)
rmax = 10;
ρmax = 10;

camera = Krang.SlowLightIntensityIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
magfield = Krang.SVector(0.0, 0.0, 1.0);
vel = Krang.SVector(0.0, 0.1, 0.0);
material = Krang.PowerLawPolarization(magfield, vel);

θs = 45 * π / 180
geometry1 = Krang.ConeGeometry(θs)
geometry2 = Krang.ConeGeometry(π-θs)
R = 4
p1 = 10
p2 = 3

profile(r) = (r/R)^p1/(1+(r/R)^(p1+p2))

iquv = Krang.observe(geometry1, material, camera; subimgs=[0], profile = profile)
i_vals = [i[1] for i in iquv]
iquv = Krang.observe(geometry2, material, camera; subimgs=[0], profile = profile)
i_vals += [i[1] for i in iquv]

hm = heatmap!(ax,i_vals, colormap=:afmhot)

iquv1 = Krang.observe(geometry1, material, camera; subimgs=[1,], profile = profile)
iquv2 = Krang.observe(geometry1, material, camera; subimgs=[1,], profile = profile)
i1_vals = Vector{Float64}()
for (i,j,k) in zip(iquv1, iquv2, i_vals)
    if i[1] + j[1] > k
        push!(i1_vals, i[1]+j[1])
    else
        push!(i1_vals, NaN)
    end
end
#i1_vals = imfilter(reshape(i1_vals, (sze, sze)), Kernel.gaussian(3))
i1_vals = reshape(i1_vals, (sze, sze))
heatmap!(ax,i1_vals, colormap=:grays)

display(figure)