module KrangKernelAbstractionsExt

# Note: Fusing pixel construction and ray-tracing one kernel does not finish compiling on Metal 
using Krang
using KernelAbstractions
using Adapt: adapt

const KA = KernelAbstractions

function _coord_vector(num, lo, hi, res)
    arr = similar(num, eltype(num), res)
    copyto!(arr, collect(range(lo, hi, length = res)))
    return arr
end

##----------------------------------------------------------------------------------------
# Screen construction
##----------------------------------------------------------------------------------------

@kernel function _fill_intensity_screen!(
    α,
    β,
    r1,
    r2,
    r3,
    r4,
    I0_inf,
    total_mino_time,
    Gθo,
    Gθhat,
    ηm,
    λm,
    @Const(αvals),
    @Const(βvals),
    met,
    θo,
)
    i, j = @index(Global, NTuple)
    @inbounds a, b = αvals[i], βvals[j]
    pix = Krang.IntensityPixel(met, a, b, θo)
    @inbounds begin
        α[i, j] = a
        β[i, j] = b
        r1[i, j], r2[i, j], r3[i, j], r4[i, j] = Krang.roots(pix)
        I0_inf[i, j] = Krang.I0_inf(pix)
        total_mino_time[i, j] = Krang.total_mino_time(pix)
        Gθo[i, j], Gθhat[i, j] = Krang.absGθo_Gθhat(pix)
        ηm[i, j] = Krang.η(pix)
        λm[i, j] = Krang.λ(pix)
    end
end

"""
    IntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res, A) where {T}

Build a fast-light screen with a `KernelAbstractions` kernel, storing each cached invariant in
an array of type `A` (e.g. `Matrix`, `MtlMatrix`, `CuMatrix`). The backend is taken from `A`,
and the result is numerically identical to the serial `src/` constructor.
"""
function Krang.IntensityScreen(
    met::Krang.Kerr{T},
    αmin,
    αmax,
    βmin,
    βmax,
    θo,
    res,
    A,
) where {T}
    real_field() = A(Matrix{T}(undef, res, res))
    cplx_field() = A(Matrix{Complex{T}}(undef, res, res))

    α, β = real_field(), real_field()
    I0_inf, total_mino_time = real_field(), real_field()
    Gθo, Gθhat = real_field(), real_field()
    ηm, λm = real_field(), real_field()
    r1, r2, r3, r4 = cplx_field(), cplx_field(), cplx_field(), cplx_field()

    αvals = _coord_vector(α, T(αmin), T(αmax), res)
    βvals = _coord_vector(α, T(βmin), T(βmax), res)

    backend = KA.get_backend(α)
    _fill_intensity_screen!(backend)(
        α,
        β,
        r1,
        r2,
        r3,
        r4,
        I0_inf,
        total_mino_time,
        Gθo,
        Gθhat,
        ηm,
        λm,
        αvals,
        βvals,
        met,
        T(θo);
        ndrange = (res, res),
    )
    KA.synchronize(backend)

    return Krang.IntensityScreen(
        (T(αmin), T(αmax)),
        (T(βmin), T(βmax)),
        met,
        T(θo),
        α,
        β,
        r1,
        r2,
        r3,
        r4,
        I0_inf,
        total_mino_time,
        Gθo,
        Gθhat,
        ηm,
        λm,
    )
end

@kernel function _fill_slow_light_screen!(
    α,
    β,
    r1,
    r2,
    r3,
    r4,
    I0_inf,
    total_mino_time,
    Iϕ_inf,
    It_inf,
    I1,
    I2,
    Ip,
    Im,
    Gθo,
    Gθhat,
    Gϕo,
    Gϕhat,
    Gto,
    Gthat,
    ηm,
    λm,
    @Const(αvals),
    @Const(βvals),
    met,
    θo,
)
    i, j = @index(Global, NTuple)
    @inbounds a, b = αvals[i], βvals[j]
    pix = Krang.SlowLightIntensityPixel(met, a, b, θo)
    @inbounds begin
        α[i, j] = a
        β[i, j] = b
        r1[i, j], r2[i, j], r3[i, j], r4[i, j] = Krang.roots(pix)
        I0_inf[i, j] = Krang.I0_inf(pix)
        total_mino_time[i, j] = Krang.total_mino_time(pix)
        Iϕ_inf[i, j] = pix.Iϕ_inf
        It_inf[i, j] = pix.It_inf
        I1[i, j] = pix.I1_inf_m_I0_terms
        I2[i, j] = pix.I2_inf_m_I0_terms
        Ip[i, j] = pix.Ip_inf_m_I0_terms
        Im[i, j] = pix.Im_inf_m_I0_terms
        Gθo[i, j], Gθhat[i, j] = pix.absGθo_Gθhat
        Gϕo[i, j], Gϕhat[i, j] = pix.absGϕo_Gϕhat
        Gto[i, j], Gthat[i, j] = pix.absGto_Gthat
        ηm[i, j] = Krang.η(pix)
        λm[i, j] = Krang.λ(pix)
    end
end

"""
    SlowLightIntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res, A) where {T}

Slow-light counterpart of the [`IntensityScreen`](@ref) kernel constructor.

!!! warning "Metal"
    This does not finish compilation on Metal.
"""
function Krang.SlowLightIntensityScreen(
    met::Krang.Kerr{T},
    αmin,
    αmax,
    βmin,
    βmax,
    θo,
    res,
    A,
) where {T}
    real_field() = A(Matrix{T}(undef, res, res))
    cplx_field() = A(Matrix{Complex{T}}(undef, res, res))

    α, β = real_field(), real_field()
    I0_inf, total_mino_time = real_field(), real_field()
    Iϕ_inf, It_inf = real_field(), real_field()
    I1, I2, Ip, Im = real_field(), real_field(), real_field(), real_field()
    Gθo, Gθhat = real_field(), real_field()
    Gϕo, Gϕhat = real_field(), real_field()
    Gto, Gthat = real_field(), real_field()
    ηm, λm = real_field(), real_field()
    r1, r2, r3, r4 = cplx_field(), cplx_field(), cplx_field(), cplx_field()

    αvals = _coord_vector(α, T(αmin), T(αmax), res)
    βvals = _coord_vector(α, T(βmin), T(βmax), res)

    backend = KA.get_backend(α)
    _fill_slow_light_screen!(backend)(
        α,
        β,
        r1,
        r2,
        r3,
        r4,
        I0_inf,
        total_mino_time,
        Iϕ_inf,
        It_inf,
        I1,
        I2,
        Ip,
        Im,
        Gθo,
        Gθhat,
        Gϕo,
        Gϕhat,
        Gto,
        Gthat,
        ηm,
        λm,
        αvals,
        βvals,
        met,
        T(θo);
        ndrange = (res, res),
    )
    KA.synchronize(backend)

    return Krang.SlowLightIntensityScreen(
        (T(αmin), T(αmax)),
        (T(βmin), T(βmax)),
        met,
        T(θo),
        α,
        β,
        r1,
        r2,
        r3,
        r4,
        I0_inf,
        total_mino_time,
        Iϕ_inf,
        It_inf,
        I1,
        I2,
        Ip,
        Im,
        Gθo,
        Gθhat,
        Gϕo,
        Gϕhat,
        Gto,
        Gthat,
        ηm,
        λm,
    )
end

##----------------------------------------------------------------------------------------
# Rendering
##----------------------------------------------------------------------------------------

# `pixels` is the lazy `*PixelArray` view, so each thread materialises its own pixel from the
# screen's field arrays. The view is not a device array and so cannot have Metal annotation.
@kernel function _render_screen!(store, pixels, scene)
    I = @index(Global, Cartesian)
    @inbounds store[I] = Krang.render(pixels[I], scene)
end

"""
    render!(store, camera::AbstractCamera, scene::Scene)

Render `scene` through `camera` into `store` with a `KernelAbstractions` kernel. The backend is
taken from `store`, so a `Matrix` runs on CPU threads and a device array runs on that device.
The camera's screen and the scene are `adapt`ed to the backend, which means an on-device screen
(see the kernel [`IntensityScreen`](@ref) constructor) is used in place rather than copied.

Each thread owns one pixel and loops over the whole scene, so meshes stay in registers.
"""
function Krang.render!(
    store::AbstractMatrix,
    camera::Krang.AbstractCamera,
    scene::Krang.Scene,
)
    pixels = camera.screen.pixels
    size(store) == size(pixels) || throw(
        DimensionMismatch("store is $(size(store)) but the screen is $(size(pixels))"),
    )
    backend = KA.get_backend(store)
    _render_screen!(backend)(
        store,
        adapt(backend, pixels),
        adapt(backend, scene);
        ndrange = size(store),
    )
    KA.synchronize(backend)
    return store
end

end # module
