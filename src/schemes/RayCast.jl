function generate_trajectory(pixel::Krang.IntensityPixel{T}, res::Int, n::Int) where T
    curr_rad, curr_inc, curr_az = T(Inf), Krang.inclination(pixel), zero(T)
    prev_x = curr_rad * sin(curr_inc)*cos(curr_az)
    prev_y = curr_rad * sin(curr_inc)*sin(curr_az)
    prev_z = curr_rad * cos(curr_inc)
    θo = Krang.inclination(pixel)
    lines = []
    for θs in 0:T(π/res):T(π)
        if (θs == θo || (π-θs) == θo) continue end
        curr_rad, curr_inc, curr_az, _, _ = Krang.emission_coordinates_fast_light(pixel, θs, Krang.screen_coordinate(pixel)[2] ≥ 0, n)
        if iszero(curr_rad) continue end
        curr_x = (curr_rad * sin(curr_inc)*cos(curr_az))
        curr_y = (curr_rad * sin(curr_inc)*sin(curr_az))
        curr_z = (curr_rad * cos(curr_inc))
        curr_line = Line(Point(prev_x, prev_y, prev_z), Point(curr_x, curr_y, curr_z))
        
        append!(lines, curr_line)
        prev_x, prev_y, prev_z = curr_x, curr_y, curr_z
    end

    return lines
end
