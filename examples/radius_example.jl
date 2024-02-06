# # Raytracing with inclination

# In this example, we will raytrace the region around a Kerr blackhole as seen by an observer stationed at infinity.
# We will show the emission coordinates of the n=0 (direct) and n=1 (indirect) photons as they are emitted from the 
# source, at a fixed inclination angle from the blackhole's spin axis.
#
# First, let's import Krang and CairoMakie for plotting.
using Krang
using CairoMakie

metric = Krang.Kerr(0.99);
θo = π-45 * π / 180;
sze = 400;
rmin = Krang.horizon(metric)
rmax = 10.0;
ρmax = 15.0;

curr_theme = Theme(
    Axis = (
        xticksvisible = false, 
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        ),
    Heatmap = (
        rasterize=true,
        colorrange = (rmin, rmax)
    )
)
set_theme!(merge!(curr_theme, theme_latexfonts()))

#
# We will use a 0.99 spin Kerr blackhole viewed by an assymptotic observer at an inclination angle of θo=π/4. 
# A region spanned by radii between the horizon and 20M at varying inclinations will be raytraced onto the 20Mx20M 
# screen of the observer.

# Create the material
fig = Figure(resolution=(700, 300), fontfamily="Computer Modern", fontface="bold");
axes = (
    Axis(fig[1, 1], aspect=1, title=L"n=0",titlesize=40),
    Axis(fig[1, 2], aspect=1, title=L"n=1",titlesize=40),
    Axis(fig[1, 3], aspect=1, title=L"n=2",titlesize=40)
)
observation = zeros(sze, sze)
camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
flag = true

function draw!(axes, camera, observation, rmin, rmax, θs, flag)
    empty!.(axes)

    material = Krang.CoordinateRadius();

    meshes = [Krang.Mesh(Krang.ConeGeometry(θs, (i, rmin, rmax)), material) for i in 0:2]
    
    hm = nothing
    for (mesh,ax) in zip(meshes, axes)
        @Threads.threads for I in CartesianIndices(camera.screen.pixels)
            observation[I] = mesh.material(camera.screen.pixels[I], mesh.geometry)
        end
        hm = heatmap!(ax, observation)
    end

    if flag
        cb = Colorbar(fig[1, 4], hm; labelsize=30, ticklabelsize=20)
        flag = false
    end
end
recording = CairoMakie.record(fig, "radius.gif", range(0.0, π, length=180), framerate=12) do θs
    draw!(axes, camera, observation, rmin, rmax, θs, flag)
end

# ![image](radius.gif)
