export render, render!
abstract type AbstractScheme end

function render(camera::AbstractCamera, scene::Scene)  

    mapreduce(
        mesh -> mesh.material.(camera.screen.pixels, Ref(mesh.geometry)),
        +,
        scene
    )
end

function render!(store, camera::AbstractCamera, scene::Scene)  
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