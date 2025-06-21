export SlowLightIntensityCamera

"""
    $TYPEDEF

Intensity Pixel Type.
"""
struct SlowLightIntensityPixel{T} <: AbstractPixel{T}
    metric::Kerr{T}
    "Pixel screen_coordinate"
    screen_coordinate::NTuple{2,T}
    "Radial roots"
    roots::NTuple{4,Complex{T}}
    "Radial anti-derivative at observer"
    I0_o::T
    "Radial anti-derivative at infinity"
    I0_inf::T
    "Total possible Mino time"
    total_mino_time::T
    "Radial phi anti-derivative"
    Iϕ_inf::T
    "Radial time anti-derivative"
    It_inf::T
    I1_o_m_I0_terms::T
    I2_o_m_I0_terms::T
    Ip_o_m_I0_terms::T
    Im_o_m_I0_terms::T
    "Angular anti-derivative"
    absGθo_Gθhat::NTuple{2,T}
    "Angular ϕ anti-derivative"
    absGϕo_Gϕhat::NTuple{2,T}
    "Angular t anti-derivative"
    absGto_Gthat::NTuple{2,T}
    "Half orbit of angular t anti-derivative"
    θo::T
    ro::T
    η::T
    λ::T
    @doc """
        SlowLightIntensityPixel(met::Kerr{T}, α::T, β::T, θo::T) where {T}

    Construct a `SlowLightIntensityPixel` object with the given Kerr metric, screen coordinates, and inclination.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric.
    - `α::T`: The Bardeen α value (screen coordinate).
    - `β::T`: The Bardeen β value (screen coordinate).
    - `θo::T`: The inclination angle.

    # Returns
    - A `SlowLightIntensityPixel` object initialized with the given parameters.

    # Details
    This function calculates the η and λ values using the provided Kerr metric and screen coordinates. 
    It then computes the radial roots and adjusts them if necessary. 
    It also calculates the radial and angular anti-derivatives. 
    Finally, it initializes a `SlowLightIntensityPixel` object with the calculated values and the provided parameters.
    """
    function SlowLightIntensityPixel(met::Kerr{T}, α::T, β::T, θo) where {T}
        tempη = Krang.η(met, α, β, θo)
        tempλ = Krang.λ(met, α, θo)
        roots = Krang.get_radial_roots(met, tempη, tempλ)
        numreals = sum(_isreal2.(roots))
        if (numreals == 2) && (abs(imag(roots[4]) / real(roots[4])) < eps(T))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        I1, I2, Ip, Im = radial_inf_integrals(met, roots)

        I0_inf = Krang.Ir_inf(met, roots)
        new{T}(
            met,
            (α, β),
            roots,
            I0_inf,
            I0_inf,
            total_mino_time(met, roots),
            Krang.Iϕ_inf(met, roots, tempλ),
            Krang.It_inf(met, roots, tempλ),
            I1,
            I2,
            Ip,
            Im,
            Krang._absGθo_Gθhat(met, θo, tempη, tempλ),
            Krang._absGϕo_Gϕhat(met, θo, tempη, tempλ),
            Krang._absGto_Gthat(met, θo, tempη, tempλ),
            θo,
            T(Inf),
            tempη,
            tempλ,
        )
    end

    function SlowLightIntensityPixel(met::Kerr{T}, longitude::T, latitude::T, θo::T, ro::T) where {T}
        @assert longitude >= -π/2 && longitude <= π/2 "Longitude must be in [-π/2, π/2]"
        a = met.spin
        red_α = sin(longitude)
        red_β = sin(latitude)
        p_local_u = [1, √abs(1-(red_α^2+red_β^2)), red_α, -red_β]
        p_bl_u = jac_bl_u_zamo_d(met, ro, θo) * p_local_u
        E, _, _, L = metric_dd(met, ro, θo) * p_bl_u
        tempλ = L/E
        tempη = (Σ(met,ro, θo)/E*p_bl_u[3])^2 - (a*cos(θo))^2 + (tempλ*cot(θo))^2

        roots = Krang.get_radial_roots(met, tempη, tempλ)
        numreals = sum(_isreal2.(roots))
        if (numreals == 2) && (abs(imag(roots[4]) / real(roots[4])) < eps(T))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        I1, I2, Ip, Im = radial_o_m_I0_terms_integrals(met, ro, roots, true)

        I0_o = Krang.Ir_s(met, ro, roots, true)
        I0_inf = Krang.Ir_inf(met, roots)
        new{T}(
            met,
            (longitude, latitude),
            roots,
            I0_o,
            I0_inf,
            total_mino_time(met, roots),
            Krang.Iϕ_o_m_I0_terms(met, ro, roots, tempλ, true),
            Krang.It_o_m_I0_terms(met, ro, roots, tempλ, true),
            I1,
            I2,
            Ip,
            Im,
            Krang._absGθo_Gθhat(met, θo, tempη, tempλ),
            Krang._absGϕo_Gϕhat(met, θo, tempη, tempλ),
            Krang._absGto_Gthat(met, θo, tempη, tempλ),
            θo,
            ro,
            tempη,
            tempλ,
        )
    end

    function SlowLightIntensityPixel(met::Kerr{T}, p_local_u::Vector{T}, θo::T, ro::T) where {T}
        a = met.spin
        p_bl_u = jac_bl_u_zamo_d(met, ro, θo) * p_local_u
        E, _, _, L = metric_dd(met, ro, θo) * p_bl_u
        tempλ = -L/E
        tempη = (Σ(met, ro, θo)/E*p_bl_u[3])^2 - (a*cos(θo))^2 + (tempλ*cot(θo))^2

        roots = Krang.get_radial_roots(met, tempη, tempλ)
        numreals = sum(_isreal2.(roots))
        if (numreals == 2) && (abs(imag(roots[4]) / real(roots[4])) < eps(T))
            roots = (roots[1], roots[4], roots[2], roots[3])
        end
        I1, I2, Ip, Im = radial_o_m_I0_terms_integrals(met, ro, roots, true)

        I0_o = Krang.Ir_s(met, ro, roots, true)
        I0_inf = Krang.Ir_inf(met, roots)
        new{T}(
            met,
            (p_local_u[3]/p_local_u[1], p_local_u[4]/p_local_u[1]),
            roots,
            I0_o,
            I0_inf,
            total_mino_time(met, roots),
            Krang.Iϕ_o_m_I0_terms(met, ro, roots, tempλ, true),
            Krang.It_o_m_I0_terms(met, ro, roots, tempλ, true),
            I1,
            I2,
            Ip,
            Im,
            Krang._absGθo_Gθhat(met, θo, tempη, tempλ),
            Krang._absGϕo_Gϕhat(met, θo, tempη, tempλ),
            Krang._absGto_Gthat(met, θo, tempη, tempλ),
            θo,
            ro,
            tempη,
            tempλ,
        )
    end

end

"""
    $TYPEDEF

Screen made of `SlowLightIntensityPixel`s.
"""
struct SlowLightIntensityScreen{T,A<:AbstractMatrix} <: AbstractScreen
    "Minimum and Maximum Bardeen α values"
    αrange::NTuple{2,T}

    "Minimum and Maximum Bardeen β values"
    βrange::NTuple{2,T}

    "Data type that stores screen pixel information"
    pixels::A

    @kernel function _generate_screen!(
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

    @kernel function _generate_screen!(
        screen,
        met::Kerr{T},
        longitude_min,
        longitude_max,
        latitude_min,
        latitude_max,
        θo,
        ro,
        res,
    ) where {T}
        I, J = @index(Global, NTuple)
        long = longitude_min + (longitude_max - longitude_min) * (T(I) - 1) / (res - 1)
        lat = latitude_min + (latitude_max - latitude_min) * (T(J) - 1) / (res - 1)
        screen[I, J] = SlowLightIntensityPixel(met, long, lat, θo, ro)
    end
    @doc """
        SlowLightIntensityScreen(met::Kerr{T}, αmin, αmax, βmin, βmax, θo, res; A=Matrix) where {T}

    Construct a `SlowLightIntensityScreen` object.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric object.
    - `αmin`: Minimum Bardeen α value.
    - `αmax`: Maximum Bardeen α value.
    - `βmin`: Minimum Bardeen β value.
    - `βmax`: Maximum Bardeen β value.
    - `θo`: Observer's inclination angle. θo ∈ (0, π).
    - `res`: Resolution of the screen (number of pixels along one dimension).
    - `A=Matrix`: Data type that stores screen pixel information (default is `Matrix`). A GPUMatrix can be used for GPU computations.

    # Returns
    A `SlowLightIntensityScreen` object.
    """
    function SlowLightIntensityScreen(
        met::Kerr{T},
        αmin,
        αmax,
        βmin,
        βmax,
        θo,
        res;
        A = Matrix,
    ) where {T}
        screen = A(Matrix{SlowLightIntensityPixel{T}}(undef, res, res))

        backend = get_backend(screen)

        _generate_screen!(backend)(
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
    function SlowLightIntensityScreen(
        met::Kerr{T},
        longitude_min,
        longitude_max,
        latitude_min,
        latitude_max,
        θo,
        ro,
        res;
        A = Matrix,
    ) where {T}
        screen = A(Matrix{SlowLightIntensityPixel{T}}(undef, res, res))

        backend = get_backend(screen)

        _generate_screen!(backend)(
            screen,
            met,
            longitude_min,
            longitude_max,
            latitude_min,
            latitude_max,
            θo,
            ro,
            res,
            ndrange = (res, res),
        )

        new{T,typeof(screen)}((longitude_min, longitude_max), (latitude_min, latitude_max), screen)
    end

end

"""
    $TYPEDEF

Camera that caches slow light raytracing information for an observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct SlowLightIntensityCamera{T,A} <: AbstractCamera
    metric::Kerr{T}
    "Data type that stores screen pixel information"
    screen::SlowLightIntensityScreen{T,A}
    "Observer screen_coordinate"
    screen_coordinate::NTuple{2,T}
    @doc """
        SlowLightIntensityCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res; A=Matrix) where {T}

    Constructs a `SlowLightIntensityCamera` object.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric.
    - `θo`: The Observer's inclination angle. θo ∈ (0, π).
    - `αmin`: Minimum α coordinate on the screen.
    - `αmax`: Maximum α coordinate on the screen.
    - `βmin`: Minimum β coordinate on the screen.
    - `βmax`: Maximum β coordinate on the screen.
    - `res`: Resolution of the screen (number of pixels along one dimension).
    - `A`: Data type that stores screen pixel information (default is `Matrix`). A GPUMatrix can be used for GPU computations.

    # Returns
    - A `SlowLightIntensityCamera` object.
    """
    function SlowLightIntensityCamera(
        met::Kerr{T},
        θo,
        αmin,
        αmax,
        βmin,
        βmax,
        res;
        A = Matrix,
    ) where {T}
        screen = SlowLightIntensityScreen(met, αmin, αmax, βmin, βmax, θo, res; A)
        new{T,typeof(screen.pixels)}(met, screen, (T(Inf), θo))
    end

    function SlowLightIntensityCamera(
        met::Kerr{T},
        θo,
        ro,
        longitude_min,
        longitude_max,
        latitude_min,
        latitude_max,
        res;
        A = Matrix,
    ) where {T}
        screen = SlowLightIntensityScreen(met, longitude_min, longitude_max, latitude_min, latitude_max, θo, ro, res; A)
        new{T,typeof(screen.pixels)}(met, screen, (ro, θo))
    end
end

function η(pix::SlowLightIntensityPixel)
    return pix.η
end
function λ(pix::SlowLightIntensityPixel)
    return pix.λ
end
function roots(pix::SlowLightIntensityPixel)
    return pix.roots
end
function screen_coordinate(pix::SlowLightIntensityPixel)
    return pix.screen_coordinate
end
function inclination(pix::SlowLightIntensityPixel)
    return pix.θo
end
function I0_inf(pix::SlowLightIntensityPixel)
    return pix.I0_inf
end
function total_mino_time(pix::SlowLightIntensityPixel)
    return pix.total_mino_time
end
function Ir_inf(pix::SlowLightIntensityPixel)
    return pix.I0_inf
end
function I1_o_m_I0_terms(pix::SlowLightIntensityPixel)
    return pix.I1_o_m_I0_terms
end
function I2_o_m_I0_terms(pix::SlowLightIntensityPixel)
    return pix.I2_o_m_I0_terms
end
function Ip_o_m_I0_terms(pix::SlowLightIntensityPixel)
    return pix.Ip_o_m_I0_terms
end
function Im_o_m_I0_terms(pix::SlowLightIntensityPixel)
    return pix.Im_o_m_I0_terms
end
function radial_o_integrals_m_I0_terms(pix::SlowLightIntensityPixel)
    return I1_o_m_I0_terms(pix),
    I2_o_m_I0_terms(pix),
    Ip_o_m_I0_terms(pix),
    Im_o_m_I0_terms(pix)
end
function Iϕ_inf(pix::SlowLightIntensityPixel)
    return pix.Iϕ_inf
end
function It_inf(pix::SlowLightIntensityPixel)
    return pix.It_inf
end
function absGθo_Gθhat(pix::SlowLightIntensityPixel)
    return pix.absGθo_Gθhat
end
function absGϕo_Gϕhat(pix::SlowLightIntensityPixel)
    return pix.absGϕo_Gϕhat
end
function absGto_Gthat(pix::SlowLightIntensityPixel)
    return pix.absGto_Gthat
end
