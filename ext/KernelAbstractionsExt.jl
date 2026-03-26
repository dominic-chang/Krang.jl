module KernelAbstractions
using Krang
using KernelAbstractions

@kernel function _generate_screen_intensity_pixel!(
    screen,
    met::Kerr{T},
    αmin,
    αmax,
    βmin,
    βmax,
    θo,
    res,
) where {T}
    I, J = @index(Global, NTuple)
    α = αmin + (αmax - αmin) * (T(I) - 1) / (res - 1)
    β = βmin + (βmax - βmin) * (T(J) - 1) / (res - 1)
    screen[I, J] = IntensityPixel(met, α, β, θo)
end

@doc """
    IntensityScreen(met::Kerr{T}, αmin::T, αmax::T, βmin::T, βmax::T, θo::T, res::Int; A=Matrix) where {T}

Creates an intensity screen for the given Kerr metric. 
This camera caches information for fast light computations.

# Arguments
- `met::Kerr{T}`: The Kerr metric.
- `αmin::T`: Minimum value of α.
- `αmax::T`: Maximum value of α.
- `βmin::T`: Minimum value of β.
- `βmax::T`: Maximum value of β.
- `θo::T`: Observer's inclination angle. θo ∈ (0, π).
- `res::Int`: Resolution of the screen.
- `A=Matrix`: Optional argument to specify the type of matrix to use. A GPUMatrix can be used for GPU computations.

# Returns
- `IntensityScreen{T, A}`: An intensity screen object.
"""
function IntensityScreen(
    met::Kerr{T},
    αmin::T,
    αmax::T,
    βmin::T,
    βmax::T,
    θo::T,
    res,
    A,
) where {T}
    screen = A(Matrix{IntensityPixel{T}}(undef, res, res))

    backend = get_backend(screen)

    _generate_screen_intensity_pixel!(backend)(
        screen,
        met,
        αmin,
        αmax,
        βmin,
        βmax,
        θo,
        res,
        ndrange = (res, res),
    )

    IntensityScreen((αmin, αmax), (βmin, βmax), screen)
end

@kernel function _generate_screen_slow_light!(
    screen,
    met::Kerr{T},
    αmin,
    αmax,
    βmin,
    βmax,
    θo,
    res,
) where {T}
    I, J = @index(Global, NTuple)
    α = αmin + (αmax - αmin) * (T(I) - 1) / (res - 1)
    β = βmin + (βmax - βmin) * (T(J) - 1) / (res - 1)
    screen[I, J] = SlowLightIntensityPixel(met, α, β, θo)
end

function SlowLightIntensityScreen(
    met::Kerr{T},
    αmin,
    αmax,
    βmin,
    βmax,
    θo,
    res,
    A
) where {T}
    screen = A(Matrix{SlowLightIntensityPixel{T}}(undef, res, res))
    backend = get_backend(screen)
    _generate_screen_slow_light!(backend)(
        screen,
        met,
        αmin,
        αmax,
        βmin,
        βmax,
        θo,
        res,
        ndrange = (res, res),
    )
    new{T,typeof(screen)}((αmin, αmax), (βmin, βmax), screen)
end

#TODO: debug this
@kernel function _generate_rays!(
    rays::AbstractArray{Intersection{T}},
    pixels::AbstractMatrix{P},
    res::Int,
) where {T,P<:AbstractPixel}
    (I, J, K) = @index(Global, NTuple)
    k = Int(K)
    pixel = pixels[I, J]
    actual_res =
        unsafe_trunc(Int, sum((Krang._isreal2.(Krang.roots(pixel))))) == 4 ? res + 1 : res
    τf = total_mino_time(pixel)

    Δτ = τf / actual_res

    ts, rs, θs, ϕs, νr, νθ = (2.0f0, 0.0f0, 0.0f0, 0.0f0, true, true)
    rays[I, J, K] = Intersection(ts, rs, θs, ϕs, νr, νθ)
    emission_coordinates(pixel, Δτ * k)


end

function generate_rays(
    pixels::AbstractMatrix{P},
    res::Int;
    A = Array,
) where {P<:AbstractPixel{T}} where {T}
    dims = (size(pixels)..., res)
    rays = A{Intersection{T}}(undef, dims...)
    backend = get_backend(rays)
    _generate_rays!(backend)(rays, pixels, res, ndrange = dims)
    return rays
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




end


