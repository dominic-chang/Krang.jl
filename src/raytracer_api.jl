function raytrace_w_inclination(mesh, pix, resolution, nrange)
    θo = inclination(pix)
    θoc = π - θo
    issouth = θo > π / 2
    int1 = []
    int2 = []
    int3 = []

    i1temp = 0
    i2temp = 0
    i3temp = 0

    Δθ = π/resolution 
    Δθ, curr = issouth ? (-Δθ, 0 + Δθ) : (Δθ, π -Δθ)
    while 0 < curr < π
        if (curr > θoc) ⊻ (θo > π / 2)
            i1temp += 1
            append!(int1, curr)
        elseif (curr > θo) ⊻ (θo > π / 2)
            i2temp += 1
            append!(int2, curr)
        else
            i3temp += 1
            append!(int3, curr)
        end
        curr -= Δθ
    end

    i1 = 0
    i2 = 0
    i3 = 0
    i4 = 0
    i5 = 0
    θsvals = []

    ans = 0
    for n in nrange[end:-1:begin]
        if n % 2 == 0
            i1 = i1temp
            i2 = 2i1
            i3 = i2 + i2temp
            i4 = i3 + i2temp
            i5 = i4 + i3temp
            θsvals = (int1[end:-1:begin]..., int1..., int2..., int2..., int3..., int3[end:-1:begin]...)
        else
            i1 = i3temp
            i2 = 2i1
            i3 = i2 + i2temp
            i4 = i3 + i2temp
            i5 = i4 + i1temp
            θsvals = (int3..., int3[end:-1:begin]..., int2[end:-1:begin]..., int2[end:-1:begin]..., int1[end:-1:begin]..., int1...)
        end

        flag = true
        for (i, θs) in enumerate(θsvals)
            curr_rs, _, _, _ = emission_radius(pix, θs, flag, n)
            ans = isnan(curr_rs) || (curr_rs > 10) ? ans : curr_rs
            if i ∈ (i1, i2, i3, i4, i5)
                flag ⊻= true
            end
        end
    end
    return ans
end