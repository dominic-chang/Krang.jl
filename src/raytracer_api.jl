function tree_function(x, y, a)
    return max(
        min(
            max(
                min(
                    (x/a - 47/48)^2 + (y/a + 19/29)^2 - 1, 
                    -(x/a - 31/77)^2 - (y/a + 23/22)^2 + 25/22, 
                    x/a
                ), 
                min(
                    (x/a + 47/48)^2 + (y/a + 19/29)^2 - 1, 
                    -(x/a + 27/67)^2 - (y/a + 23/22)^2 + 25/22, 
                    -x/a
                )
            ), 
                -y/a - 1/2
        ), 
        min(
            max(
                min(
                    (x/a - 16/19)^2 + (y/a + 56/67)^2 - 49/25, 
                    -(x/a - 8/19)^2 - (y/a + 27/20)^2 + 31/16, 
                    x/a
                ), 
                min(
                    (x/a + 16/19)^2 + (y/a + 56/67)^2 - 49/25, 
                    -(x/a + 8/19)^2 - (y/a + 27/20)^2 + 31/16, 
                    -x/a
                )
            ), 
            -y/a - 1/2
        ), 
        min(
            -(x/a + 17/37)^2 - (y/a + 3/5)^2 + 13/28, 
            y^2/a^2 + (x/a + 1)^2 - 1, 
            -x/a
        ), 
        min(
            -(x/a - 17/37)^2 - (y/a + 3/5)^2 + 13/28, 
            y^2/a^2 + (x/a - 1)^2 - 1, 
            x/a, 
            -y/a
        ), 
        min(
            1/6 - x/a, 
            x/a + 1/6, 
            1/3*(-y/a - 5/2), 
            1/3*(y/a + 7/2)
        )
    )
end
function raytrace_w_inclination(mesh, pix, resolution, nrange)
    op = mesh.opening_angle
    θo = inclination(pix)
    θoc = π - θo
    issouth = θo > π / 2
    ending = issouth ? π : 0
    Δθ = π / resolution
    Δθ, curr = issouth ? (-Δθ, Δθ) : (Δθ, π - Δθ)
    i1 = ceil((curr - θoc) / Δθ)
    curr - i1 * Δθ
    θoc
    if (curr - (i1 - 1) * Δθ ≥ θoc) ⊻ (θo > π / 2)
        #i1 -= 1
        Δθ *= 1 - 1e-5
    end
    if (curr - (i1) * Δθ ≥ θoc) ⊻ (θo > π / 2)
        Δθ *= 1 + 1e-5
    end
    i2 = ceil((curr - i1 * Δθ - θo) / Δθ)
    i3 = ceil((curr - (i1 + i2) * Δθ - ending - Δθ) / Δθ)

    ans = 0
    for n in nrange[end:-1:begin]
        boundaries = n % 2 == 0 ? [
            [(curr - (i1 - 1) * Δθ, true, 1, i1), (curr, false, -1, i1)],
            [(curr - i1 * Δθ, false, -1, i2), (curr - (i1) * Δθ, true, -1, i2)],
            [(curr - (i1 + i2) * Δθ, true, -1, i3), (curr - (i1 + i2 + i3 - 1) * Δθ, false, 1, i3)]
        ] : [
            [(curr - (i1 + i2) * Δθ, true, -1, i3), (curr - (i1 + i2 + i3 - 1) * Δθ, false, 1, i3)],
            [(curr - (i1 + i2 - 1) * Δθ, false, 1, i2), (curr - (i1 + i2 - 1) * Δθ, true, 1, i2)],
            [(curr - (i1 - 1) * Δθ, true, 1, i1), (curr, false, -1, i1)],
        ]
        for boundary in boundaries
            for (temp, flag, sign, curri) in boundary
                count = 0
                while count < curri
                    x2 = Krang.emission_radius(pix, temp, flag, n)
                    #ans = (!isnan(x2[1]) && x2[1] < 10) ? x2[1] : ans
                    ans = (!isnan(x2[1]) && tree_function(x2[1]*sin(temp), x2[1]*cos(temp) - 3, 3) > 0) ? x2[1] : ans
                    count += 1
                    temp += sign * Δθ
                end
            end
        end
    end

    return ans
end

