export render, render!, render_cpu_threaded!
abstract type AbstractScheme end

"""
    render!(store, camera::AbstractCamera, scene::Scene)

Render `scene` through `camera` into `store`, dispatching on the backend `store` lives on.

Provided by `KrangKernelAbstractionsExt`: load `KernelAbstractions` (plus the package owning
the array type, e.g. `Metal` or `CUDA`) to get a method. A `Matrix` store runs on CPU threads,
a device array runs on that device.
"""
function render! end

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

# The optical depth has done its job once every mesh has been accumulated, so only the
# intensity is handed back.
function render(::AbsorptiveIntensityTrait, pixel::AbstractPixel, scene::Scene)
    return _render(_render_init(pixel, scene), pixel, scene).I
end

function _render(observation, pixel::AbstractPixel, scene::Scene)
    mesh = scene[1]

    for itr = 1:length(scene)
        mesh = scene[itr]
        observation = _accumulate(observation, raytrace(pixel, mesh))
    end
    return observation
end

function render_cpu_threaded!(store, camera::AbstractCamera, scene::Scene)
    @assert size(store) == size(camera.screen.pixels)
    Threads.@threads for I in CartesianIndices(camera.screen.pixels)
        pixel = camera.screen.pixels[I]
        store[I] = sum(mesh.material(pixel, mesh.geometry) for mesh in scene)
    end
    return store
end
