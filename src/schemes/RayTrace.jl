export generate_ray

function generate_ray!(ray::Matrix{T}, pixel::Krang.AbstractPixel, res::Int) where T
    numreals = Int(sum((Krang._isreal2.(Krang.roots(pixel)))))
    τf = zero(T)
    actual_res = res
    if numreals == 4 
        τf = 2Krang.I0_inf(pixel)
        actual_res += 1 # The last point on scattering rays is at infinity
    else
        rh = Krang.horizon(pixel.metric)
        radial_roots = roots(pixel)
        r1, r2, _, r4 = radial_roots
        r21, r31, r32, r41, r42, _ = _get_root_diffs(radial_roots...)
        r1, r2, r21 = real.((r1, r2, r21))

        if numreals == 2
            A = √abs(r32 * r42)
            B = √abs(r31 * r41)
            k = (((A + B)^2 - r21^2) / (4 * A * B))
            temprat = B * (rh - r2) / (A * (rh - r1))
            x3_s = clamp(((one(T) - temprat) / (one(T) + temprat)), -one(T), one(T))
            coef = one(T) * √inv(A * B)
            τf = I0_inf(pixel) - coef * JacobiElliptic.F((acos(x3_s)), k)
        else
            C = √abs(r31 * r42)
            D = √abs(r32 * r41)
            k4 = 4C * D / (C + D)^2
            a2 = abs(imag(r1))
            b1 = real(r4)
        
            k4 = T(4) * C * D / (C + D)^2
        
            go = √max((T(4)a2^2 - (C - D)^2) / ((C + D)^2 - T(4)a2^2), zero(T))
            x4_s = (rh + b1) / a2
            coef = 2 / (C + D)
            τf = I0_inf(pixel) - coef*JacobiElliptic.F(atan(x4_s) + atan(go), k4)
        end
    end

    Δτ = T(τf/(actual_res))
    for I in range(1, res)

        _, curr_rad, curr_inc, curr_az, _, _ = emission_coordinates(pixel, Δτ*I)
        ray[1,I] = curr_rad * sin(curr_inc)*cos(curr_az)
        ray[2,I] = curr_rad * sin(curr_inc)*sin(curr_az)
        ray[3,I] = curr_rad * cos(curr_inc)

    end

end

function generate_ray(pixel::Krang.AbstractPixel{T}, res::Int) where T
    ray = zeros(T, 3, res)
    generate_ray!(ray, pixel, res)
    return ray
end
