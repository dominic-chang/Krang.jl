using Krang
using CairoMakie
using Metal

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

metric = Krang.Kerr(0.99f0); # Kerr metric with a spin of 0.99
θo = 45f0 * π / 180; # inclination angle of the observer. θo ∈ (0, π)
sze = 400; # resolution of the screen is sze x sze
rmin = Krang.horizon(metric); # minimum radius to be ray traced
rmax = 10f0; # maximum radius to be ray traced
ρmax = 15f0; # horizontal and vertical limits of the screen

# Create Figure
fig = Figure(resolution=(700, 700));
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
colormaps = (:afmhot, :afmhot, :hsv)
colorrange = ((-20, 20), (0, rmax), (0, 2π))

function coordinate_point(pix::Krang.AbstractPixel, geometry::Krang.ConeGeometry{T,A}) where {T, A}
    n, rmin, rmax = geometry.attributes
    θs = geometry.opening_angle

    coords = ntuple(j -> zero(T), Val(4))

    isindir = false 
    for _ in 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        ts, rs, θs, ϕs =  emission_coordinates(pix, geometry.opening_angle, isindir, n)
        if rs ≤ rmin || rs ≥ rmax
            continue
        end
        coords = (ts, rs, θs, ϕs)
    end
    return coords 
end

# Draw Function
# This function draws the coordinates associated with the n=0,1,2 subimages of a cone with opening angle θs.
function draw!(axes_list, camera, coordinates, rmin, rmax, θs)
    times, radii, azimuths = coordinates 
    map(axes -> empty!.(axes), axes_list)

    geometries = (Krang.ConeGeometry(θs, (i, rmin, rmax)) for i in 0:2)
    
    for (i, geometry) in enumerate(geometries)
        rendered_scene = Array(coordinate_point.(MtlArray(camera.screen.pixels), Ref(geometry)))
        for I in CartesianIndices(rendered_scene)
            times[I] = rendered_scene[I][1]
            radii[I] = rendered_scene[I][2]
            azimuths[I] = rendered_scene[I][4]
        end
        coordinates = (times, radii, mod2pi.(azimuths ))
        for j in 1:3
            heatmap!(axes_list[i][j], coordinates[j], colormap = colormaps[j], colorrange=colorrange[j])
        end
    end
end

# Create the animation of Cone of Emission Coordinates
recording = CairoMakie.record(fig, "coordinate.gif", range(0f0, 1f0π, length=180), framerate=12) do θs
    draw!(axes_list, camera, coordinates, rmin, rmax, θs)
end

