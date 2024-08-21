function line_intersection(line::Line, cone::ConeGeometry)
    origin, lp2 = line.points
    return atan(hypot(origin[1], origin[2]), origin[3]) ≤ cone.opening_angle ≤ atan(hypot(lp2[1], lp2[2]), lp2[3])
end