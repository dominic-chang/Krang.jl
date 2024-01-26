using Krang
#using WGLMakie

metric = Krang.Kerr(0.9);
θo = 45 * π / 180;
sze = 400;
rmin = Krang.horizon(metric)
rmax = 10.0;
ρmax = 10.0;
σ = 1.0
χ = π/2
ι = π/2
η1 = π/4
βv = 0.9
σζ = 3.0

camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
η2 = π-η1
magfield1 = Krang.SVector(sin(ι)*cos(η1), sin(ι)*sin(η1), cos(ι));
magfield2 = Krang.SVector(sin(ι)*cos(η2), sin(ι)*sin(η2), cos(ι));
vel = Krang.SVector(βv, (π/2), χ);
material = Krang.PowerLawPolarization();

θs = (40 * π / 180)
function profile(r)
    R = 4.0
    p1 = 3.0
    p2 = 2.0
    return (r/R)^p1/(1+(r/R)^(p1+p2))
end
geometry1 = Krang.ConeGeometry((θs), (magfield1, vel, (0,1,2), profile, σ, σζ))
geometry2 = Krang.ConeGeometry((π-θs), (magfield2, vel, (0,1,2), profile, σ, σζ))
geometry = geometry1 ⊕ geometry2

mesh = Krang.Mesh(geometry, material)
mesh.material.(camera.screen.pixels, Ref(mesh.geometry))
