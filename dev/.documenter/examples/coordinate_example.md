


# Raytracing with inclination {#Raytracing-with-inclination}

In this example, we will raytrace the region around a Kerr blackhole as seen by an observer stationed at infinity. We will show the emission coordinates of the n=0 (direct) and n=1 (indirect) photons as they are emitted from the source, at a fixed inclination angle from the blackhole&#39;s spin axis.

First, let&#39;s import Krang and CairoMakie for plotting.

```julia
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

metric = Krang.Kerr(0.99);
θo = 45 * π / 180;
sze = 400;
rmin = Krang.horizon(metric)
rmax = 10.0;
ρmax = 15.0;
```


We will use a 0.99 spin Kerr blackhole viewed by an assymptotic observer at an inclination angle of θo=π/4. A region spanned by radii between the horizon and 20M at varying inclinations will be raytraced onto the 20Mx20M screen of the observer.

Create Figure

```julia
fig = Figure(resolution=(700, 700));
axes_list = [
    [
        Axis(fig[i, 1], title=(i==1 ? "Regularized Time" : ""), titlesize=20, ylabel=(i==1 ? L"n=0" : i==2 ? L"n=1" : L"n=2"), ylabelsize=20),
        Axis(fig[i, 2], title=(i==1 ? "Radius" : ""), titlesize=20),
        Axis(fig[i, 3], title=(i==1 ? "Azimuth" : ""), titlesize=20),
    ] for i in 1:3
]
```


```
3-element Vector{Vector{Makie.Axis}}:
 [Axis (0 plots), Axis (0 plots), Axis (0 plots)]
 [Axis (0 plots), Axis (0 plots), Axis (0 plots)]
 [Axis (0 plots), Axis (0 plots), Axis (0 plots)]
```


Initialize Camera and Pre-Allocate Memory for data to be plotted

```julia
coordinates = (zeros(sze, sze) for _ in 1:3)
camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
material = Krang.CoordinatePoint();
colormaps = (:afmhot, :afmhot, :hsv)
colorrange = ((-20, 20), (0, rmax), (0, 2π))
```


```
((-20, 20), (0, 10.0), (0, 6.283185307179586))
```


Draw Function

```julia
function draw!(axes_list, camera, material, coordinates, rmin, rmax, θs)
    times, radii, azimuths = coordinates
    map(axes -> empty!.(axes), axes_list)

    meshes = [Krang.Mesh(Krang.ConeGeometry(θs, (i, rmin, rmax)), material) for i in 0:2]

    for i in 1:3
        @Threads.threads for I in CartesianIndices(camera.screen.pixels)
            times[I], radii[I], _, azimuths[I] = meshes[i].material(camera.screen.pixels[I], meshes[i].geometry)
        end
        coordinates = (times, radii, azimuths)
        for j in 1:3
            heatmap!(axes_list[i][j], coordinates[j], colormap = colormaps[j], colorrange=colorrange[j])
        end
    end
end
```


```
draw! (generic function with 1 method)
```


Create the animation of Cone of Emission Coordinates

```julia
recording = CairoMakie.record(fig, "coordinate.gif", range(0.0, π, length=180), framerate=12) do θs
    draw!(axes_list, camera, material, coordinates, rmin, rmax, θs)
end
```


```
"coordinate.gif"
```



![](coordinate.gif)



---


_This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl)._
