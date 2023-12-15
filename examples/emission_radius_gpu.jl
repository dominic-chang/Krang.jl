using CairoMakie
using Krang
using Metal

metric = Krang.Kerr(Float32(0.99));
θo = Float32(85 * π / 180);
θs = Float32(90 * π / 180);
sze = 500;
rmin = Krang.horizon(metric)
rmax = Float32(15);
ρmax = Float32(15);
n = 0

#camera = Krang.SlowLightCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
gpucamera = Krang.BasicGPUCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
camera = Krang.BasicCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);

mtlarr = gpucamera.screen.pixels
arr = camera.screen.pixels
CairoMakie.heatmap(
    [isnan(x[1]) || x[1] > 10.0 ? 0.0 : x[1] for x in collect(Krang.emission_radius.(mtlarr, θs, true, n))] .+ 
    [isnan(x[1]) || x[1] > 10.0 ? 0.0 : x[1] for x in collect(Krang.emission_radius.(mtlarr, θs, false, n))]
)
