struct CoordinateRadius<: AbstractMaterial end

function (prof::CoordinateRadius)(pix::AbstractPixel, geometry::ConeGeometry{T,A}) where {T, A}
    n, rmin, rmax = geometry.attriributes
    observation = zero(T)

    isindir = false
    for _ in 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        rs, _, _ =  emission_radius(pix, geometry.opening_angle, isindir, n)
        if rs ≤ rmin || rs ≥ rmax
            continue
        end
        observation = isnan(rs) ? observation : rs
    end
    return observation
end
