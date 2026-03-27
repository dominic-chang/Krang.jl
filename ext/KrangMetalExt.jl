module KrangMetalExt
using Krang
using Metal

function _custom_pow(x::ComplexF32, y::Float32) 
    s,c = sincos(y*angle(x)*one(x))
    return (abs(x)^y)*(c-s*one(x)im)
end

function Krang.get_radial_roots(metric::Krang.Kerr{T}, η::T, λ::T) where {T<:Float32}
    a = metric.spin

    a2 = a^2
    A = a2 - η - λ^2
    A2 = A + A
    B = T(2) * (η + (λ - a)^2)
    C = -a2 * η

    P = -A^2 / T(12) - C
    Q = -A / T(3) * (A * A / T(36) + zero(T)im - C) - B^2 / T(8)

    negΔ3 = T(4) * P * P * P + T(27) * (Q^2)
    ωp =  _custom_pow(-Q / T(2) + sqrt(negΔ3 / T(108)) + zero(T)im,  T(1 / 3))

    #C = ((-1+0im)^(2/3), (-1+0im)^(4/3), 1) .* ωp
    C = (-T(1 / 2) + (√T(3) / 2)im, -T(1 / 2) - T(√T(3) / 2)im, one(T) + zero(T)im) .* ωp

    v = -P .* _custom_pow.(T(3) .* C, -one(T))

    ξ0 = argmax(real, (C .+ v)) - A / T(3)
    ξ02 = 2ξ0

    predet1 = A2 + ξ02
    predet2 = (sqrt(T(2)) * B) * _custom_pow(sqrt(ξ0), -one(T))
    det1 = sqrt(-(predet1 - predet2))
    det2 = sqrt(-(predet1 + predet2))

    sqrtξ02 = sqrt(ξ02)

    r1 = (-sqrtξ02 - det1) / 2
    r2 = (-sqrtξ02 + det1) / 2
    r3 = (sqrtξ02 - det2) / 2
    r4 = (sqrtξ02 + det2) / 2

    roots = (r1, r2, r3, r4)
    if (sum(Krang._isreal2, roots) == 2) && (abs(imag(roots[4])) < sqrt(eps(T)))
        roots = (roots[1], roots[4], roots[2], roots[3])
    end
    return roots
end

end