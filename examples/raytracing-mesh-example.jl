# # Raytracing a polygon mesh
using Krang
import GLMakie as GLMk
using GeometryBasics
using FileIO

metric = Krang.Kerr(0.99) # Kerr metric with a spin of 0.99
θo = 90/180*π # Inclination angle of the observer
ρmax = 12.0 # Horizontal and Vertical extent of the screen
sze = 100 # Resolution of the screen is sze x sze

camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze)

# Import external mesh
bunny_mesh = translate(
    rotate(
        scale(
            load(download("https://graphics.stanford.edu/~mdfisher/Data/Meshes/bunny.obj", "bunny.obj")), 100
        ),π/2, 1.0, 0.0, 0.0
    ), 2.0,7.0,-10.0
);

# Raytrace the mesh embedded in the Kerr space time. The emission picked up by a ray will be the sum of all the intersections of the ray with the mesh
# The mesh is embedded in the space-time by assuming the world coordinates of the mesh is the same as the Boyer-Lindquist coordinates of the space time.
# The ray tracing is done by scanning over the image one line at a time

# Let us now raytrace the image to see what the mesh looks like from the observer's perspective.

intersections = raytrace(camera, bunny_mesh);

# And plot the image with GLMakie, 

fig = GLMk.Figure()
ax = GLMk.Axis(fig[1,1], aspect=1)  
GLMk.heatmap!(ax, intersections, colormap=:afmhot);

save("mesh_geometry_example.png", fig);

# ![image](mesh_geometry_example.png)

# Below is a short plotting routine that generates a video showing the scene from three different perspectives, and live renders the image.

fig = GLMk.Figure()
ax1 = GLMk.Axis3(fig[1,1], aspect=(1,1,1), title="Top view of scene", elevation = π/2, azimuth = π)
ax2 = GLMk.Axis3(fig[1,2], aspect=(1,1,1), title="Side view of scene", elevation = 0.0, azimuth = 3π/2)
ax3 = GLMk.Axis3(fig[2,1], aspect=(1,1,1), title="Isometric view of scene")

ax = GLMk.Axis(fig[2,2], aspect=1, title="Heatmap of ray intersections with mesh")

for a in [ax1, ax2, ax3]
    GLMk.xlims!(a, (-10, 10)) 
    GLMk.ylims!(a, (-10, 10)) 
    GLMk.zlims!(a, (-10, 10)) 
    GLMk.hidedecorations!(a)
end

GLMk.hidedecorations!(ax)
sphere = GLMk.Sphere(GLMk.Point(0.0,0.0,0.0), horizon(metric)) # Sphere to represent black hole
lines_to_plot = Krang.generate_ray.(camera.screen.pixels, 100) # 100 is the number of steps to take along the ray

img = zeros(sze, sze)
recording = GLMk.record(fig, "mesh.mp4", 1:sze*sze, framerate=120) do i
    line = lines_to_plot[i] 

    img[i] = intersections[i]

    GLMk.empty!(ax)
    GLMk.empty!(ax1)
    GLMk.empty!(ax2)
    GLMk.empty!(ax3)

    for a in [ax1, ax2, ax3]
        GLMk.mesh!(a, bunny_mesh, color=[parse(GLMk.Colorant, "rgba(0%, 50%, 70%,1.0)") for tri in bunny_mesh.position], colormap = :blues, transparency=true)
        GLMk.mesh!(a, sphere, color=:black) 
    end

    GLMk.lines!(ax3, line, color=:red)
    GLMk.heatmap!(ax, img, colormap=:afmhot, colorrange=(0, 8))
end
# ```@raw html 
# <video loop muted playsinline controls src="./mesh.mp4" />
# ```

