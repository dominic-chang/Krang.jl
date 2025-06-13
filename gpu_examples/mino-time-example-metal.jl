using Krang
import GLMakie as GLMk
using Metal
GLMk.Makie.inline!(true)

metric = Krang.Kerr(0.99f0); # Kerr spacetime with 0.99 spin
θo = 85.0f0 * π / 180; # Observer inclination angle with respect to spin axis
sze = 200; # Number of pixels along each axis of the screen
ρmax = 5.0f0; # Size of the screen

camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);

curr_theme = GLMk.Theme(# Makie theme
    fontsize = 20,
    Axis = (
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
        topspinevisible = false,
        bottomspinevisible = false,
        titlefontsize = 30,
    ),
)

GLMk.set_theme!(GLMk.merge(curr_theme, GLMk.theme_latexfonts()))

fig = GLMk.Figure(resolution = (500, 600));

recording =
    GLMk.record(fig, "raytrace.gif", range(0.1f0, 3.0f0, length = 290), framerate = 15) do τ
        GLMk.empty!(fig)

        coordinates = Array(
            ((x, τ) -> emission_coordinates(x, τ)[1:4]).(MtlArray(camera.screen.pixels), τ),
        )
        time = [i[1] for i in coordinates]
        radius = [i[2] for i in coordinates]
        inclination = [i[3] for i in coordinates]
        azimuth = mod2pi.([i[4] for i in coordinates])

        data = (time, radius, inclination, azimuth)
        titles = (
            GLMk.L"\text{Regularized Time }(t_s)",
            GLMk.L"\text{Radius }(r_s)",
            GLMk.L"\text{Inclination }(\theta_s)",
            GLMk.L"\text{Azimuth } (\phi_s)",
        )
        colormaps = (:afmhot, :afmhot, :afmhot, :hsv)
        colorrange = ((-20, 20), (0, 10), (0, π), (0, 2π))
        indices = ((1, 1), ())

        for i = 1:4
            hm = GLMk.heatmap!(
                GLMk.Axis(
                    getindex(fig, (i > 2 ? 2 : 1), (iszero(i % 2) ? 3 : 1));
                    aspect = 1,
                    title = titles[i],
                ),
                data[i],
                colormap = colormaps[i],
                colorrange = colorrange[i],
            )
            cb = GLMk.Colorbar(
                fig[(i > 2 ? 2 : 1), (iszero(i % 2) ? 3 : 1)+1],
                hm;
                labelsize = 30,
                ticklabelsize = 20,
            )
        end

        ax = GLMk.Axis(fig[3, 1:3], height = 60)
        GLMk.hidedecorations!(ax)
        GLMk.text!(ax, 0, 100; text = GLMk.L"θ_o=%$(Int(floor(θo*180/π)))^\circ")
        GLMk.rowgap!(fig.layout, 1, GLMk.Fixed(0))
    end

# ![image](raytrace.gif)

camera = Krang.SlowLightIntensityCamera(metric, θo, -3, 3, -3, 3, 4);

fig = GLMk.Figure()
ax = GLMk.Axis3(fig[1, 1], aspect = (1, 1, 1))
GLMk.xlims!(ax, (-3, 3))
GLMk.ylims!(ax, (-3, 3))
GLMk.zlims!(ax, (-3, 3))
lines_to_plot =
    Krang.generate_rays(MtlArray(camera.screen.pixels), 5_000; A = MtlArray) |> Array

sphere = GLMk.Sphere(GLMk.Point(0.0, 0.0, 0.0), horizon(metric))
GLMk.mesh!(ax, sphere, color = :black) # Sphere to represent black hole

for i = 1:4
    for j = 1:4
        GLMk.lines!(ax, lines_to_plot[i, j, :, :])
    end
end
fig

GLMk.save("mino_time_rays.png", fig)

# ![Photons trajectories around Kerr black hole in Boyer-Lindquist Coordinates](mino_time_rays.png)
