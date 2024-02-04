using Krang
#using WGLMakie

metric = Krang.Kerr(0.9f0);
θo = 18f0 * π / 180;
sze = 400;
rmin = Krang.horizon(metric)
ρmax = 15f0;

χ = 2.7242920822576653f0
ι = 1.0113763707982746f0
βv = 0.3096525555378556f0
σ = -0.34430247937749187f0
η1 = 0.01817958773540323f0

camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
η2 = Float32(π-η1)
magfield1 = Krang.SVector(sin(ι)*cos(η1), sin(ι)*sin(η1), cos(ι));
magfield2 = Krang.SVector(sin(ι)*cos(η2), sin(ι)*sin(η2), cos(ι));
vel = Krang.SVector(βv, (π/2f0), χ);
material = Krang.ElectronSynchrotronPowerLawIntensity();

θs = (40f0 * π / 180)
function profile(r)
    R = 5.59261856012756f0
    p1 = 0.158968560099307f0
    p2 = 4.192568509151882f0

    return (r/R)^p1/(1+(r/R)^(p1+p2))
end
geometry1 = Krang.ConeGeometry((θs), (magfield1, vel, (0,1,2), profile, σ))
geometry2 = Krang.ConeGeometry((π-θs), (magfield2, vel, (0,1,2), profile, σ))
geometry = geometry1 ⊕ geometry2# ⊕ geometry1

mesh = Krang.Mesh(geometry, material)
using Metal
pixels = MtlArray(camera.screen.pixels)
using CairoMakie
arr = MtlArray(zeros(Float32, sze, sze))

@Metal.time arr = mesh.material.(pixels, Ref(mesh.geometry)) 
fig = Figure(resolution=(600, 600), fontfamily="Computer Modern", fontface="bold");
ax = Axis(fig[1,1], aspect=1)
heatmap!(ax, arr |> Array, colormap=:afmhot)
display(fig)
#mesh.material.(metal(camera.screen.pixels), Ref(mesh.geometry))

