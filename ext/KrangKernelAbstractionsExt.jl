module KrangKernelAbstractionsExt

using Krang
using KernelAbstractions

@kernel function _generate_screen_intensity_pixel!(
    screen,
    met::Krang.Kerr{T},
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
    screen[I, J] = Krang.IntensityPixel(met, α, β, θo)
end

function Krang.IntensityScreen(
    met::Krang.Kerr{T},
    αmin::T,
    αmax::T,
    βmin::T,
    βmax::T,
    θo::T,
    res,
    A,
) where {T}
    screen = A(Matrix{Krang.IntensityPixel{T}}(undef, res, res))
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

    return Krang.IntensityScreen{typeof(screen)}((αmin, αmax), (βmin, βmax), screen)
end

@kernel function _generate_screen_slow_light!(
    screen,
    met::Krang.Kerr{T},
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
    screen[I, J] = Krang.SlowLightIntensityPixel(met, α, β, θo)
end

function Krang.SlowLightIntensityScreen(
    met::Krang.Kerr{T},
    αmin::T,
    αmax::T,
    βmin::T,
    βmax::T,
    θo::T,
    res,
    A,
) where {T}
    screen = A(Matrix{Krang.SlowLightIntensityPixel{T}}(undef, res, res))
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

    return Krang.SlowLightIntensityScreen{typeof(screen)}(
        (αmin, αmax),
        (βmin, βmax),
        screen,
    )
end

@kernel function _generate_rays!(
    rays::AbstractArray{Krang.Intersection{T}},
    pixels::AbstractMatrix{P},
    res::Int,
) where {T,P<:Krang.AbstractPixel}
    I, J, K = @index(Global, NTuple)
    k = Int(K)
    pixel = pixels[I, J]
    actual_res =
        unsafe_trunc(Int, sum(Krang._isreal2.(Krang.roots(pixel)))) == 4 ? res + 1 : res
    τf = Krang.total_mino_time(pixel)
    Δτ = τf / actual_res

    ts, rs, θs, ϕs, νr, νθ, _ = Krang.emission_coordinates(pixel, Δτ * k)
    ϕs = ϕs % T(2π)
    rays[I, J, K] = Krang.Intersection(T(ts), T(rs), T(θs), T(ϕs), νr, νθ)
end

function Krang.generate_rays(
    pixels::AbstractMatrix{<:Krang.AbstractPixel},
    res::Int;
    A = Array,
) 
    T = typeof(Krang.metric(first(pixels)).spin)
    dims = (size(pixels)..., res)
    rays = A{Krang.Intersection{T}}(undef, dims...)
    backend = get_backend(rays)
    _generate_rays!(backend)(rays, pixels, res, ndrange = dims)
    return rays
end

@kernel function _render!(store, pixels, mesh::Krang.Mesh)
    I = @index(Global)
    @inbounds store[I] = mesh.material(pixels[I], mesh.geometry)
end

function Krang.render!(store, camera::Krang.AbstractCamera, scene::Krang.Scene)
    backend = get_backend(store)
    @assert backend == get_backend(camera.screen.pixels)
    @assert size(store) == size(camera.screen.pixels)
    fill!(store, zero(eltype(store)))
    temp = similar(store)

    for mesh in scene
        _render!(backend)(temp, camera.screen.pixels, mesh, ndrange = size(store))
        store .+= temp
    end

    return store
end

end
