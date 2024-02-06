struct CoordinatePoint<: AbstractMaterial end

function (prof::CoordinatePoint)(pix::AbstractPixel, geometry::ConeGeometry{T,A}) where {T, A}
    n, rmin, rmax = geometry.attriributes
    θs = geometry.opening_angle

    observation = zeros(T, 4)

    isindir = false 
    for _ in 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        ts, rs, θs, ϕs =  emission_coordinates(pix, geometry.opening_angle, isindir, n)
        if rs ≤ rmin || rs ≥ rmax
            continue
        end
        observation = isnan(rs) ? observation :  (ts, rs, θs, ϕs)
    end
    return observation
end
