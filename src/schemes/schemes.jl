export render, render!, render_cpu_threaded!
abstract type AbstractScheme end

function render(camera::AbstractCamera, scene::Scene)
    return render.(camera.screen.pixels, Ref(scene))
end

function render(pixel::AbstractPixel, scene::Scene)
    render(returnTrait(scene[1].material), pixel, scene)
end

function render(::AbstractReturnTrait, pixel::AbstractPixel{T}, scene::Scene) where {T}
    return _render(zero(T), pixel, scene)
end

function render(
    ::AbstractPolarizationTrait,
    pixel::AbstractPixel{T},
    scene::Scene,
) where {T}
    return _render(StokesParams(zero(T), zero(T), zero(T), zero(T)), pixel, scene)
end

function _render(observation, pixel::AbstractPixel{T}, scene::Scene) where {T}
    mesh = scene[1]

    for itr = 1:length(scene)
        mesh = scene[itr]
        observation += raytrace(pixel, mesh)
    end
    return observation
end

@kernel function _render!(store, pixels, mesh::Mesh)
    I = @index(Global)
    @inbounds store[I] = mesh.material(pixels[I], mesh.geometry)
end

function render!(store, camera::AbstractCamera, scene::Scene)
    backend = get_backend(store)
    @assert backend == get_backend(camera.screen.pixels)
    @assert size(store) == size(camera.screen.pixels)
    mapreduce(
        mesh -> begin
            _render!(backend)(store, camera.screen.pixels, mesh, ndrange = size(store))
            store
        end,
        +,
        scene,
    )
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
