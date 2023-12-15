# # Raytracing with inclination

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
metric = Krang.Kerr(0.99);
θo = 85 * π / 180;
sze = 500;
rmin = Krang.horizon(metric)
rmax = 10;
ρmax = 10;

camera = Krang.SlowLightCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);

curr_theme = Theme(
    fontsize=20,
    Axis=(
        xticksvisible=false,
        xticklabelsvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false,
        leftspinevisible=false,
        rightspinevisible=false,
        topspinevisible=false,
        bottomspinevisible=false,
        titlefontsize=30,
    ),
)

set_theme!(curr_theme)

# Let us now create a figure to plot the emission coordinates on,
fig = Figure(resolution=(800, 300));
# and use this figure make an animation by looping over the inclination angle θs.
# This loop will plot the emission coordinates for each θs.
recording = record(fig, "emission_coordinates.gif", range(0, π, length=180), framerate=15) do θs
    empty!(fig)

    time, radius, inclination, azimuth = [size(camera.screen.pixels) |> zeros for i in 1:4]

    for n in [1,0]
        Threads.@threads for I in CartesianIndices(time)
            coordinates = Krang.emission_coordinates(camera.screen.pixels[I], θs, true, n)
            if !any(isnan.(coordinates[1:4])) && coordinates[2] < rmax
                time[I], radius[I], inclination[I], azimuth[I] = coordinates
            end

            coordinates = Krang.emission_coordinates(camera.screen.pixels[I], θs, false, n)
            if !any(isnan.(coordinates[1:4])) && coordinates[2] < rmax
                time[I], radius[I], inclination[I], azimuth[I] = coordinates
            end
        end

        Threads.@threads for I in CartesianIndices(time)
            coordinates = Krang.emission_coordinates(camera.screen.pixels[I], θs, true, n)
            if !any(isnan.(coordinates[1:4])) && coordinates[2] < rmax
                time[I], radius[I], inclination[I], azimuth[I] = coordinates
            end

            coordinates = Krang.emission_coordinates(camera.screen.pixels[I], θs, false, n)
            if !any(isnan.(coordinates[1:4])) && coordinates[2] < rmax
                time[I], radius[I], inclination[I], azimuth[I] = coordinates
            end
        end
    end

    data = (time, radius, azimuth)
    titles = (L"\text{Regularized Time }(t_s)", L"\text{Radius }(r_s)", L"\text{Azimuth } (\phi_s)")
    colormaps = (:afmhot, :afmhot, :hsv)
    colorrange = ((-20, 20), (0, rmax), (0, 2π))

    for i in 1:3
        hm = heatmap!(
            Axis(getindex(fig, 1, (2i)-1); aspect=1, title=titles[i]),
            data[i],
            colormap=colormaps[i],
            colorrange=colorrange[i]
        )
        cb = Colorbar(fig[1, 2i], hm; ticklabelfont="Computer Modern", labelsize=30, ticklabelsize=20)
    end

    ax = Axis(fig[2, 3], height=60)
    hidedecorations!(ax)
    CairoMakie.text!(ax,0,100; text=L"θ_s=%$(Int(floor(θs*180/π)))^\circ")

    display(fig)
end

# ![image](emission_coordinates.gif)
