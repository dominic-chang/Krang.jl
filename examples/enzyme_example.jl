using Revise 
using Krang
using Enzyme
Enzyme.API.printall!(false)
Enzyme.Compiler.RunAttributor[] = false
Enzyme.Compiler.CheckNan[] = true

function test(α, β)
    a = 0.94
    η1 = π/4
    ι = π/2
    βv = 0.5
    χ = π-η1
    θs = π/2
    R = 5.0
    p1 = 4.0
    σ = 3.0
    p2 = 3.0
    θo = π/4
    metric = Krang.Kerr(a);
    η2 = 1f0π-η1
    magfield1 = Krang.SVector(sin(ι)*cos(η1), sin(ι)*sin(η1), cos(ι));
    magfield2 = Krang.SVector(sin(ι)*cos(η2), sin(ι)*sin(η2), cos(ι));
    vel = Krang.SVector(βv, (π/2f0), χ);
    material = Krang.ElectronSynchrotronPowerLawIntensity();

    @inline profile(r) = let R=R, p1=p1, p2=p2
        return (r/R)^p1/(1f0+(r/R)^(p1+p2))
    end

    geometry1 = Krang.ConeGeometry((θs), (magfield1, vel, (0,1,2), profile, σ))
    geometry2 = Krang.ConeGeometry((1f0π-θs), (magfield2, vel, (0,1,2), profile, σ))
    geometry = geometry1 ⊕ geometry2

    mesh = Krang.Mesh(geometry, material)
    pix = Krang.IntensityPixel(metric, α, β, θo)
    return mesh.material(pix, mesh.geometry)
end


test(1.0, 5.0)

autodiff(ReverseWithPrimal, test, Active, Active(1.0), Active(5.0))