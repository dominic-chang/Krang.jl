module KrangMetalExt
using Krang
using Metal

@inline _complex_inv(z::ComplexF32) = conj(z) / abs2(z)
@inline _complex_cbrt(z::ComplexF32) = abs(z)^(1f0/3) * cis(angle(z) / 3f0)

@inline function _argmax_real3(x1, x2, x3)
    if real(x1) >= real(x2)
        return real(x1) >= real(x3) ? x1 : x3
    end
    return real(x2) >= real(x3) ? x2 : x3
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
    ωp = _complex_cbrt(-Q / T(2) + sqrt(negΔ3 / T(108)) + zero(T)im)

    #C = ((-1+0im)^(2/3), (-1+0im)^(4/3), 1) .* ωp
    C = (-T(1 / 2) + (√T(3) / 2)im, -T(1 / 2) - T(√T(3) / 2)im, one(T) + zero(T)im) .* ωp

    v = -P .* _complex_inv.(T(3) .* C)

    ξ0 = _argmax_real3(C[1] + v[1], C[2] + v[2], C[3] + v[3]) - A / T(3)
    ξ02 = 2ξ0

    predet1 = A2 + ξ02
    predet2 = (sqrt(T(2)) * B) * _complex_inv(sqrt(ξ0))
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
