# # Raytracing with inclination

# In this example, we will raytrace the region around a Kerr blackhole as seen by an observer stationed at infinity.
# We will show the emission coordinates of the n=0 (direct) and n=1 (indirect) photons as they are emitted from the 
# source, at a fixed inclination angle from the blackhole's spin axis.
#
# First, let's import Krang and CairoMakie for plotting.
using Krang
using CairoMakie

curr_theme = Theme(
    Axis = (
        xticksvisible = false, 
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        aspect=1
        ),
    Heatmap = (
        rasterize=true,
    )
)
set_theme!(merge!(curr_theme, theme_latexfonts()))

metric = Krang.Kerr(0.5);
θo = 1 * π / 180;
sze = 4000;
rmin = 2.5#Krang.horizon(metric)
rmax = 3.0;
ρmax = 7.0;

using Optimization

#
# We will use a 0.99 spin Kerr blackhole viewed by an assymptotic observer at an inclination angle of θo=π/4. 
# A region spanned by radii between the horizon and 20M at varying inclinations will be raytraced onto the 20Mx20M 
# screen of the observer.

# Create Figure
fig = Figure(resolution=(1400, 1400));
axes_list = [
    [
        Axis(fig[i, 1], title=(i==1 ? "Regularized Time" : ""), titlesize=20, ylabel=(i==1 ? L"n=0" : i==2 ? L"n=1" : L"n=2"), ylabelsize=20),
        Axis(fig[i, 2], title=(i==1 ? "Radius" : ""), titlesize=20),
        Axis(fig[i, 3], title=(i==1 ? "Azimuth" : ""), titlesize=20),
    ] for i in 1:3
]


# Initialize Camera and Pre-Allocate Memory for data to be plotted
coordinates = (zeros(sze, sze) for _ in 1:3)
camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
material = Krang.CoordinatePoint();
colormaps = (:afmhot, :afmhot, :hsv)
colorrange = ((-20, 20), (0, rmax), (0, 2π))

# Draw Function
function draw!(axes_list, camera, material, coordinates, rmin, rmax, θs)
    times, radii, azimuths = coordinates 
    map(axes -> empty!.(axes), axes_list)

    meshes = [Krang.Mesh(Krang.ConeGeometry(θs, (i, rmin, rmax)), material) for i in 0:2]
    
    for i in 1:3
        @Threads.threads for I in CartesianIndices(camera.screen.pixels)
            t, r, _, az = meshes[i].material(camera.screen.pixels[I], meshes[i].geometry)
            
            azimuths[I] = az == 0 ? NaN : az
            times[I] = t == 0 ? NaN : t
            radii[I] = r == 0 ? NaN : r
            #times[I], radii[I], _, azimuths[I] = meshes[i].material(camera.screen.pixels[I], meshes[i].geometry)
        end
        coordinates = (times, radii, azimuths)
        for j in 1:3
            heatmap!(axes_list[i][j], coordinates[j], colormap = colormaps[j], colorrange=colorrange[j], nan_color=:black)
        end
    end
end

# Create the animation of Cone of Emission Coordinates
#recording = CairoMakie.record(fig, "coordinate.gif", range(0.0, π, length=180), framerate=12) do θs
θs = π/2
draw!(axes_list, camera, material, coordinates, rmin, rmax, θs)
#end
display(fig)

# ![image](coordinate.gif\sq

nan2zero(x) = isnan(x) ? 0.0 : x
using Roots
function get_isoradial_curve(rs, φ, θs, θo, met, n)
    function f(ρ)
        α = ρ*cos(φ)
        β = ρ*sin(φ)
        pix = Krang.SlowLightIntensityPixel(met, α, β, θo)
        rref = nan2zero(Krang.emission_coordinates(pix, θs, true, n)[2]) + nan2zero(Krang.emission_coordinates(pix, θs, false, n)[2])
        (rref-rs)
    end
    find_zero(f, [0.95√27,1.5√27])
end

rs = 7
met = Krang.Kerr(0.001)
θo = 1/180*π
φvals = range(2π*0.001, 2π*0.999, length = 1000)
ρvals = get_isoradial_curve.(rs, φvals, π/2, θo, Ref(met), 0)
xvals = ρvals .* cos.(φvals)
yvals = ρvals .* sin.(φvals)

fig = Figure()
ax = Axis(fig[1,1], aspect=1)
xlims!(ax, -20,20)
ylims!(ax, -20,20)
lines!(ax, xvals, yvals)
display(fig)