# # Raytracing a polygon mesh
using Krang
import GLMakie as GLMk
using GeometryBasics
using FileIO

metric = Krang.Kerr(0.99) # Kerr metric with a spin of 0.99
θo = 90/180*π # Inclination angle of the observer
ρmax = 12.0 # Horizontal and Vertical extent of the screen
sze = 400 # Resolution of the screen is sze x sze

camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze)

# Import external mesh
bunny_mesh = translate(
    rotate(
        scale(
            load(download("https://graphics.stanford.edu/~mdfisher/Data/Meshes/bunny.obj", "bunny.obj")), 100
        ),π/2, 1.0, 0.0, 0.0
    ), 2.0,6.0,-10.0
)

# Raytrace the mesh embedded in the Kerr space time. The emission picked up by a ray will be the sum of all the intersections of the ray with the mesh
# The mesh is embedded in the space-time by assuming the world coordinates of the mesh is the same as the Boyer-Lindquist coordinates of the space time.
# Below is some code that generates a gif showing a global view of what the intersection between the mesh and each ray looks like.

fig = GLMk.Figure()
ax = GLMk.Axis3(fig[1,1], aspect=(1,1,1))#, elevation = π/2, azimuth = π)
sphere = GLMk.Sphere(GLMk.Point(0.0,0.0,0.0), horizon(metric))
lines_to_plot = Krang.generate_ray.(camera.screen.pixels, 10)

recording = GLMk.record(fig, "mesh.gif", range(1,100), framerate=15) do i
    line = lines_to_plot[i] |> collect
    GLMk.empty!(fig)
    ax = GLMk.Axis3(fig[1,1], aspect=(1,1,1))#, elevation = π/2, azimuth = π)
    GLMk.xlims!(ax, (-10, 10)) 
    GLMk.ylims!(ax, (-10, 10)) 
    GLMk.zlims!(ax, (-10, 10)) 

    GLMk.mesh!(ax, bunny_mesh, color=[parse(GLMk.Colorant, "rgba(0%, 50%, 70%,1.0)") for tri in bunny_mesh.position], colormap = :blues, transparency=true, shading=true)
    GLMk.mesh!(ax, sphere, color=:black) # Sphere to represent black hole

    GLMk.lines!(ax, line, color=:red)
end
    
fig

# ![image](mesh.gif)


# Let us now raytrace the image to see what the mesh looks like from the observer's perspective.

intersections = raytrace(camera, bunny_mesh)


# And plot the image with GLMakie, 

fig = GLMk.Figure()
ax = GLMk.Axis(fig[1,1], aspect=1, xreversed=true)  
GLMk.heatmap!(ax, intersections, colormap=:afmhot)
fig

save("mesh_geometry_example.png", fig)

# ![image](mesh_geometry_example.png)