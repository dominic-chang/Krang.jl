# # Raytracing a Cone

# In this example, we will raytrace the region around a Kerr blackhole as seen by an observer stationed at infinity.
# We will show the emission coordinates of the n=0 (direct) and n=1 (indirect) photons as they are emitted from the 
# source, at a fixed inclination angle from the blackhole's spin axis.
#
# First, let's import Krang and CairoMakie for plotting.
using CairoMakie
using Krang
#
# We will use a 0.99 spin Kerr blackhole viewed by an assymptotic observer at an inclination angle of θo=π/4. 
# A region spanned by radii between the horizon and 20M at varying inclinations will be raytraced onto the 20Mx20M 
# screen of the observer.
metric = Krang.Kerr(1.0);
sze = 50;
rmin = Krang.horizon(metric)
rmax = 10;
ρmax = 10;
#observer = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);

# Let us now create a figure to plot the emission coordinates on,

fig = Figure(resolution=(600, 600), fontfamily="Computer Modern", fontface="bold");
# and use this figure make an animation by looping over the inclination angle θs.
# This loop will plot the emission coordinates for each θs.
#αvals = range(-ρmax, ρmax, length=sze);
#βvals = range(-ρmax, ρmax, length=sze);

θs = π/2
θo = 85π/180
n = 0
camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
CairoMakie.heatmap(
    [isnan(x[1]) || x[1] >10.0 ? 0.0 : x[1] for x in Krang.emission_radius.(camera.screen.pixels, θs, true, n)] .+ 
    [isnan(x[1]) || x[1] > 10.0 ? 0.0 : x[1] for x in Krang.emission_radius.(camera.screen.pixels, θs, false, n)]
)


