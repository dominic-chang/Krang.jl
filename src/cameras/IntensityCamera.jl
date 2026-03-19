export IntensityCamera
"""
    $TYPEDEF

Intensity Pixel Type. 
Each Pixel is associated with a single ray. 
Pixels cache some information about the ray.
"""
struct IntensityPixel{T, SCREEN_TYPE, ROOTS_TYPE, G0_TYPE} <: AbstractPixel{T}
	metric::Kerr{T}
	"Bardeen coordiantes if ro is \`Inf`, otherwise longitude and latitude"
	screen_coordinate::SCREEN_TYPE#AbstractArray{T}
	"Radial roots"
	roots::ROOTS_TYPE#SCREEN_TYPEAbstractArray{Complex{T}}
	"Radial anti-derivative at observer"
	I0_o::T
	"Radial anti-derivative at infinity"
	I0_inf::T
	"Total possible Mino time"
	total_mino_time::T
	"Angular antiderivative"
	absGθo_Gθhat::G0_TYPE
	"Inclination"
	θo::T
	"radius"
	ro::T
	η::T
	λ::T
	function _IntensityPixel(met::Kerr{T}, tempη, tempλ, θo, ro, x, y)  where T
		roots = Krang.get_radial_roots(met, tempη, tempλ)
		numreals = sum(_isreal2, roots)
		cond = (numreals == 2) & (abs(imag(roots[4]) / real(roots[4])) < eps(met.spin))
		roots = Base.ifelse(cond, (roots[1], roots[4], roots[2], roots[3]), roots)
		_I0_o = Krang.Ir_s(met, ro, roots, true)
		I0_inf = Krang.Ir_inf(met, roots)
		total_mino_time = _I0_o + I0_inf
        
        if numreals != 4 
            total_mino_time = _I0_o - Krang.Ir_s(met, Krang.horizon(met), roots, true)
        end

		screen_coords = [x, y]
		Go = collect(Krang._absGθo_Gθhat(met, θo, tempη, tempλ))
		new{T, typeof(screen_coords), typeof(roots), typeof(Go)}(
			met,
			screen_coords,
			roots,
			_I0_o,
			I0_inf,
			total_mino_time,
			Go,
			θo,
			ro,
			tempη,
			tempλ,
		)
	end

    @doc """
        IntensityPixel(met::Kerr{T}, α::T, β::T, θo::T) where {T}

    Construct an `IntensityPixel` object with the given Kerr metric, screen coordinates, and inclination.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric.
    - `α::T`: The Bardeen α value (screen coordinate).
    - `β::T`: The Bardeen β value (screen coordinate).
    - `θo::T`: The inclination angle.

    # Returns
    - An `IntensityPixel` object initialized with the given parameters.

	# Details
	This function calculates the η and λ values using the provided Kerr metric and screen coordinates. 
	It then computes the radial roots and adjusts them if necessary. 
	Finally, it initializes an `IntensityPixel` object with the calculated values and the provided parameters.
	"""
	function IntensityPixel(met::Kerr{A}, α, β, θo) where A
		tempη = Krang.η(met, α, β, θo)
		tempλ = Krang.λ(met, α, θo)
		_IntensityPixel(met, tempη, tempλ, θo, Inf, α, β)
	end




end

"""
    $TYPEDEF

Screen made of `IntensityPixel`s.
"""
struct IntensityScreen{T,A<:AbstractMatrix} <: AbstractScreen
    "Minimum and Maximum x coordinate values"
    xrange::NTuple{2,T}

    "Minimum and Maximum y coordinate values"
    yrange::NTuple{2,T}

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
        screen[I, J] = IntensityPixel(met, α, β, θo)
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
        long = latitude_min + (latitude_max - latitude_min) * (T(I) - 1) / (res - 1)
        lat = longitude_min + (longitude_max - longitude_min) * (T(J) - 1) / (res - 1)
        screen[I, J] = IntensityPixel(met, long, lat, θo, ro)
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
        res;
        A = Matrix,
    ) where {T}
        screen = A(Matrix{IntensityPixel{T}}(undef, res, res))

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
    function IntensityScreen(
        met::Kerr{T},
        longitude_min::T,
        longitude_max::T,
        latitude_min::T,
        latitude_max::T,
        θo::T,
        ro::T,
        res;
        A = Matrix,
    ) where {T}
        screen = A(Matrix{IntensityPixel{T}}(undef, res, res))

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

Camera that caches fast light raytracing information for an observer sitting at radial infinity.
The frame of this observer is alligned with the Boyer-Lindquist frame.
"""
struct IntensityCamera{T,A} <: AbstractCamera
    metric::Kerr{T}
    "Data type that stores screen pixel information"
    screen::IntensityScreen{T,A}
    "camera location"
    camera_location::NTuple{2,T}

    @doc """
        IntensityCamera(met::Kerr{T}, θo, αmin, αmax, βmin, βmax, res::Int; A=Matrix) where {T}

    Constructor for an Intensity Camera.

    # Arguments
    - `met::Kerr{T}`: The Kerr metric.
    - `θo`: Observer's inclination angle. θo ∈ (0, π).
    - `αmin`: Minimum value of α.
    - `αmax`: Maximum value of α.
    - `βmin`: Minimum value of β.
    - `βmax`: Maximum value of β.
    - `res::Int`: Resolution of the screen.
    - `A=Matrix`: Optional argument to specify the type of matrix to use. A GPUMatrix can be used for GPU computations.

    # Returns
    - `IntensityCamera{T, A}`: An intensity camera object.
    """
    function IntensityCamera(
        met::Kerr{T},
        θo,
        αmin,
        αmax,
        βmin,
        βmax,
        res::Int;
        A = Matrix,
    ) where {T}

        screen = IntensityScreen(met, αmin, αmax, βmin, βmax, θo, res; A = A)
        new{T,typeof(screen.pixels)}(
            met,
            screen,
            (T(Inf), θo),
        )
    end
    function IntensityCamera(
        met::Kerr{T},
        θo,
        ro,
        longitude_min,
        longitude_max,
        latitude_min,
        latitude_max,
        res::Int;
        A = Matrix,
    ) where {T}

        screen = IntensityScreen(met, longitude_min, longitude_max, latitude_min, latitude_max, θo, ro, res; A = A)
        new{T,typeof(screen.pixels)}(
            met,
            screen,
            (ro, θo),
        )
    end

end

function η(pix::IntensityPixel)
    return pix.η
end
function λ(pix::IntensityPixel)
    return pix.λ
end
function roots(pix::IntensityPixel)
    return pix.roots
end
function screen_coordinate(pix::IntensityPixel)
    return pix.screen_coordinate
end
function inclination(pix::IntensityPixel)
    return pix.θo
end
function I0_o(pix::IntensityPixel)
    return pix.I0_o
end
function total_mino_time(pix::IntensityPixel)
    return pix.total_mino_time
end

function absGθo_Gθhat(pix::IntensityPixel)
    return pix.absGθo_Gθhat
end
