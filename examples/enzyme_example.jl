using Revise
using Krang
using Enzyme
using ForwardDiff
Enzyme.API.printall!(false)
Enzyme.Compiler.RunAttributor[] = false
Enzyme.Compiler.CheckNan[] = true

struct TestMaterial <: Krang.AbstractMaterial end
function (mat::TestMaterial)(pix::Krang.AbstractPixel, geometry::Krang.ConeGeometry{T,A}) where {T,A}
    isindir = false
    ans = zero(T)
    for _ in 1:2
        isindir ⊻= true
        rs, _, _ = emission_radius(pix, geometry.opening_angle, isindir, 0)
        if rs ≤ horizon(Krang.metric(pix)) || isinf(rs)
            continue
        end
        ans += rs
    end
    return ans
end

function test(α, β)
    a = 0.94
    θo = π/4
    θs = π/2
    metric = Krang.Kerr(a);
    material = TestMaterial();
    geometry = Krang.ConeGeometry((θs),)


    mesh = Krang.Mesh(geometry, material)
    pix = Krang.IntensityPixel(metric, α, β, θo)
    #return Krang.emission_radius(pix , θs, false, 0)[1]

    return mesh.material(pix, mesh.geometry)
end



function test(α, β)
    a = 0.94
    η1 = π/4
    ι = π/2
    βv = 0.5
    χ = π-η1
    θs = π/2.1
    R = 3.0
    p1 = 4.0
    σ = 3.0
    p2 = 3.0
    θo = π/4.0
    metric = Krang.Kerr(a);
    η2 = 1f0π-η1
    magfield1 = Krang.SVector(sin(ι)*cos(η1), sin(ι)*sin(η1), cos(ι));
    magfield2 = Krang.SVector(sin(ι)*cos(η2), sin(ι)*sin(η2), cos(ι));
    vel = Krang.SVector(βv, (π/2f0), χ);
    material = Krang.ElectronSynchrotronPowerLawIntensity();
#
    @inline profile(r) = let R=R, p1=p1, p2=p2
        return (r/R)^p1/(1f0+(r/R)^(p1+p2))
    end
#
    geometry1 = Krang.ConeGeometry((θs), (magfield1, vel, (0,1), profile, σ))
    geometry2 = Krang.ConeGeometry((1f0π-θs), (magfield2, vel, (0,), profile, σ))
    geometry = geometry1 ⊕ geometry2
#
    mesh = Krang.Mesh(geometry, material)
    pix = Krang.IntensityPixel(metric, α, β, θo)
    #return Krang.emission_radius(pix , θs, false, 0)[1]
    return mesh.material(pix, geometry)
end


test(1.0, 5.0)

using CairoMakie
curr_theme = Theme(
    Axis = (
        xticksvisible = false, 
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        ),
    Heatmap = (
        rasterize=true,
    )
)
set_theme!(merge!(curr_theme, theme_latexfonts()))
fig = Figure()

Ivals = []
ax1 = Axis(fig[1,1], aspect=1, title="Intensity")
for β in -10:0.11:10
    for α in -10:0.11:10
        append!(Ivals, test(α,β))
    end
end
heatmap!(ax1, reshape(Ivals, 182, 182))
display(fig)

Ivals = []
ax2 = Axis(fig[1,2], aspect=1, title="Enzyme (Reverse)")
for β in -10:0.11:10
    for α in -10:0.11:10
        append!(Ivals, hypot(autodiff(Enzyme.Reverse, test, Active, Active(α), Active(β))[1]...))
    end
end
heatmap!(ax2, reshape(Ivals, 182, 182))
display(fig)

Ivals2 = []
ax3 = Axis(fig[2,1], aspect=1, title="ForwardDiff")
for β in -10:0.11:10
    for α in -10:0.11:10
        append!(Ivals2, hypot(ForwardDiff.gradient((a)->test(a[1],a[2]),[α,β])...))
    end
end
heatmap!(ax3, reshape(Ivals2, 182, 182))
display(fig)


ax4 = Axis(fig[2,2], aspect=1, title="Difference")
heatmap!(ax4, reshape((Ivals .- Ivals2) ./ Ivals, 182, 182))
display(fig)
