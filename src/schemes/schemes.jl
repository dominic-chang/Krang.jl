export render, render!, render_cpu_threaded!
abstract type AbstractScheme end

function render(camera::AbstractCamera, scene::Scene)  

    mapreduce(
        mesh -> mesh.material.(camera.screen.pixels, Ref(mesh.geometry)),
        +,
        scene
    )
end

@kernel function _render!(store, pixels, mesh::Mesh)
    I = @index(Global)
    for I in CartesianIndices(pixels)
        @inbounds store[I] = mesh.material(pixels[I], mesh.geometry)
    end
end

function render!(store, camera::AbstractCamera, scene::Scene)  
    backend = get_backend(store)
    @assert backend == get_backend(camera.screen.pixels)
    @assert size(store) == size(camera.screen.pixels)
    mapreduce(
        mesh -> begin _render!(backend)(store, camera.screen.pixels, mesh, ndrange=size(store)); store end,
        +,
        scene
    )
end

function render_cpu_threaded!(store, camera::AbstractCamera, scene::Scene)  
    @assert size(store) == size(camera.screen.pixels)
    mapreduce(
        mesh -> begin
            @Threads.threads for I in CartesianIndices(camera.screen.pixels)
                store[I] = mesh.material(camera.screen.pixels[I], mesh.geometry)
            end
            store
        end,
        +,
        scene
    )
end