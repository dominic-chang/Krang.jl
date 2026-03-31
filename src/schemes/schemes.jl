export render, render!, render_cpu_threaded!
abstract type AbstractScheme end

function render(camera::AbstractCamera, scene::Scene)
    return render.(camera.screen.pixels, Ref(scene))
end

function render(pixel::AbstractPixel, scene::Scene)
    render(returnTrait(scene[1].material), pixel, scene)
end

function _render_init(pixel::AbstractPixel, scene::Scene)
    return zero(raytrace(pixel, scene[1]))
end

function render(::AbstractReturnTrait, pixel::AbstractPixel, scene::Scene)
    return _render(_render_init(pixel, scene), pixel, scene)
end

function render(::AbstractPolarizationTrait, pixel::AbstractPixel, scene::Scene)
    return _render(_render_init(pixel, scene), pixel, scene)
end

function _render(observation, pixel::AbstractPixel, scene::Scene) where {T}
    mesh = scene[1]

    for itr = 1:length(scene)
        mesh = scene[itr]
        observation += raytrace(pixel, mesh)
    end
    return observation
end

function render_cpu_threaded!(store, camera::AbstractCamera, scene::Scene)
    @assert size(store) == size(camera.screen.pixels)
    mapreduce(
        mesh -> begin
            Threads.@threads for I in CartesianIndices(camera.screen.pixels)
                store[I] = mesh.material(camera.screen.pixels[I], mesh.geometry)
            end
            store
        end,
        +,
        scene,
    )
end
