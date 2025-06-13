# # Coordinates with inclination (θs)

# This example shows how to access coordinate information from the raytracing process with inclination based functions.
# We will ray trace the coordinate points on a sequence of cones in the region around a Kerr black hole as seen by an observer stationed at infinity.
# We will show the emission coordinates of the n=0 (direct) and n=1 (indirect) photons that are emitted from the 
# source, at a fixed inclination angles with respect to the black hole's spin axis.
# ## Setup
# First, let's import Krang and CairoMakie for plotting.
using Krang
using CairoMakie

#set the theme for the plots
curr_theme = Theme(
    Axis = (
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        aspect = 1,
    ),
    Heatmap = (rasterize = true,),
)
set_theme!(merge!(curr_theme, theme_latexfonts()))

# We will use a 0.99 spin Kerr black hole viewed by an asymptotic observer at an inclination angle of θo=π/4. 
# A region spanned by radii between the horizon and 10M at varying inclinations will be ray traced onto the 15x15 
# screen of the observer.

metric = Krang.Kerr(0.01); # Kerr metric with a spin of 0.99
θo = 45 * π / 180; # inclination angle of the observer. θo ∈ (0, π)
ro  = 20.0
sze = 40; # resolution of the screen is sze x sze
rmin = Krang.horizon(metric); # minimum radius to be ray traced
rmax = 6.0; # maximum radius to be ray traced
ρmax = 10/ro; # horizontal and vertical limits of the screen


# Initialize Camera and pre-allocate memory for data to be plotted
coordinates = (zeros(sze, sze) for _ = 1:3)
camera = Krang.SlowLightIntensityCamera(metric, θo, ro, -ρmax, ρmax, -ρmax, ρmax, sze);
camera2 = Krang.SlowLightIntensityCamera(metric, θo, -10, 10, -10, 10, sze);
colormaps = (:afmhot, :afmhot, :hsv)
colorrange = ((-20, 20), (0, rmax), (0, 2π))


[i.screen_coordinate[1] for i in camera2.screen.pixels]  .≈ [i.screen_coordinate[1] .* ro for i in camera.screen.pixels]

maximum([abs.(camera2.screen.pixels[i].η - camera.screen.pixels[i].η) for i in range(1, length(camera.screen.pixels))])
maximum([abs.(camera2.screen.pixels[i].λ - camera.screen.pixels[i].λ) for i in range(1, length(camera.screen.pixels))])


## Defining a function to get the coordinates of the geometry
# Let's define a function that will return the coordinates of a ray when it intersects with a cone of opening angle $\theta_s$.
# Coordinate information can be accessed using the `emission_coordinates(pixel, θs, isindir, n)` function, which returns the 
# coordinate information for a sub-image `n` associated with the ray at the pixel `pixel` which intersects a cone with opening angle `θs`.
# We will define a function which introduces some basic occlusion effects by checking if the ray has intersected with the cone on 
# the 'far-side' or the 'near-side' from the observer.
function coordinate_point(
    pix::Krang.AbstractPixel,
    geometry::Krang.ConeGeometry{T,A},
) where {T,A}
    n, rmin, rmax = geometry.attributes
    θs = geometry.opening_angle

    coords = ntuple(j -> zero(T), Val(4))

    isindir = false
    for _ = 1:2 # Looping over isindir this way is needed to get Metal.jl to work
        isindir ⊻= true
        ts, rs, ϕs = emission_coordinates(pix, geometry.opening_angle, isindir, n)
        if rs ≤ rmin || rs ≥ rmax
            continue
        end
        coords = isnan(rs) ? observation : (ts, rs, θs, ϕs)
    end
    return coords
end

# ## Drawing coordinate points
# This function plots the coordinates associated with the n=0,1,2 sub-images of a cone with opening angle θs.
function draw!(axes_list, camera, coordinates, rmin, rmax, θs)
    times, radii, azimuths = coordinates
    map(axes -> empty!.(axes), axes_list)

    geometries = (Krang.ConeGeometry(θs, (i, rmin, rmax)) for i = 0:2)

    for (i, geometry) in enumerate(geometries)
        rendered_scene = coordinate_point.(camera.screen.pixels, Ref(geometry))
        for I in CartesianIndices(rendered_scene)
            temp = rendered_scene[I][1]
            times[I] =  temp 
            radii[I] = rendered_scene[I][2]
            azimuths[I] = mod2pi(rendered_scene[I][4]) # azimuths are in radians
        end
        coordinates = (times, radii, azimuths)
        for j = 1:3
            heatmap!(
                axes_list[i][j],
                coordinates[j],
                colormap = colormaps[j],
                colorrange = colorrange[j],
            )
        end
    end
end

function draw2!(axes_list, camera, coordinates, rmin, rmax, θs)
    times, radii, azimuths = coordinates
    map(axes -> empty!.(axes), axes_list)

    geometries = (Krang.ConeGeometry(θs, (i, rmin, rmax)) for i = 0:2)

    for (i, geometry) in enumerate(geometries)
        rendered_scene = coordinate_point.(camera.screen.pixels, Ref(geometry))
        for I in CartesianIndices(rendered_scene)
            temp = rendered_scene[I][1]
            times[I] =  temp  == 0 ? 0 : temp - 2log(ro) -ro
            radii[I] = rendered_scene[I][2]
            azimuths[I] = mod2pi(rendered_scene[I][4]) # azimuths are in radians
        end
        coordinates = (times, radii, azimuths)
        for j = 1:3
            heatmap!(
                axes_list[i][j],
                coordinates[j],
                colormap = colormaps[j],
                colorrange = colorrange[j],
            )
        end
    end
end


# Create Figure
fig = Figure(resolution = (700, 700));
axes_list = [
    [
        Axis(
            fig[i, 1],
            title = (i == 1 ? "Regularized Time" : ""),
            titlesize = 20,
            ylabel = (i == 1 ? L"n=0" : i == 2 ? L"n=1" : L"n=2"),
            ylabelsize = 20,
        ),
        Axis(fig[i, 2], title = (i == 1 ? "Radius" : ""), titlesize = 20),
        Axis(fig[i, 3], title = (i == 1 ? "Azimuth" : ""), titlesize = 20),
    ] for i = 1:3
]

draw2!(axes_list, camera, coordinates, rmin, rmax, π/2)
display(fig)

fig = Figure(resolution = (700, 700));
axes_list = [
    [
        Axis(
            fig[i, 1],
            title = (i == 1 ? "Regularized Time" : ""),
            titlesize = 20,
            ylabel = (i == 1 ? L"n=0" : i == 2 ? L"n=1" : L"n=2"),
            ylabelsize = 20,
        ),
        Axis(fig[i, 2], title = (i == 1 ? "Radius" : ""), titlesize = 20),
        Axis(fig[i, 3], title = (i == 1 ? "Azimuth" : ""), titlesize = 20),
    ] for i = 1:3
]

draw!(axes_list, camera2, coordinates, rmin, rmax, π/2)
display(fig)



# Create the animation of Cone of Emission Coordinates
#recording = CairoMakie.record(
#    fig,
#    "coordinate.gif",
#    range(0.0, π, length = 180),
#    framerate = 12,
#) do θs
#    draw!(axes_list, camera, coordinates, rmin, rmax, θs)
#end

# ![Emission coordinates of cones](coordinate.gif)

# > [!IMPORTANT]
# > The GPU can be used in this example with an appropriate broadcast.
# 
# ```julia
# using CUDA
# 
# rendered_scene = Array(coordinate_point.(CuArray(camera.screen.pixels), Ref(geometry)))
# ```
