export p_bl_d,
    penrose_walker,
    screen_polarization,
    PowerLaw,
    evpa,
    jac_zamo_d_bl_u,
    jac_bl_d_zamo_u,
    jac_zamo_u_bl_d,
    jac_bl_u_zamo_d,
    jac_fluid_u_zamo_d
##----------------------------------------------------------------------------------------------------------------------
#Polarization stuff
##----------------------------------------------------------------------------------------------------------------------
"""
    Returns the momentum form in the Boyer-Lindquist basis.
"""
function p_bl_d(metric::Kerr{T}, r, θ, η, λ, νr::Bool, νθ::Bool) where {T}
    return @SVector[
        -one(T),
        (νr ? one(T) : -one(T)) * √max(zero(T), r_potential(metric, η, λ, r)) /
        Δ(metric, r),
        (νθ ? one(T) : -one(T)) * √max(zero(T), θ_potential(metric, η, λ, θ)),
        λ,
    ]
end

"""
Returns the tetrad that transforms the Boyer-Lindquist covector to the ZAMO basis.

The tetrad is (Eq 3.1 1972ApJ...178..347B):

```
    √(A/(ΣΔ))   0        0        2ar/√(AΔΣ)
    0           √(Δ/Σ)   0        0
    0           0        0        √(Σ/A)csc(θ)
    0           0       -1/√Σ     0
```
"""
function jac_zamo_d_bl_u(metric::Kerr{T}, r, θ) where {T}
    a = metric.spin
    # coords = {t, r, θ, ϕ}
    Σt = Σ(metric, r, θ)
    Δt = Δ(metric, r)
    At = A(metric, r, θ)

    return @SMatrix [# Eq 3.1 1972ApJ...178..347B
        √(At / (Σt * Δt)) zero(T) zero(T) 2*a*r/√(At * Δt * Σt)
        zero(T) √(Δt / Σt) zero(T) zero(T)
        zero(T) zero(T) zero(T) √(Σt / At)*csc(θ)
        zero(T) zero(T) -inv(√Σt) zero(T)
    ]
end

"""
Transforms the Boyer-Lindquist covector to the ZAMO basis.

Applies the tetrad above without forming it.
"""
@inline function jac_zamo_d_bl_u(metric::Kerr{T}, r, θ, v_d) where {T}
    Δt = Δ(metric, r)
    Σt = Σ(metric, r, θ)
    Ξt = Ξ(metric, r, θ; Δ = Δt)
    ωt = T(2) * metric.spin * r / Ξt

    vt_d, vr_d, vθ_d, vϕ_d = v_d
    return @SVector [
        (vt_d + ωt * vϕ_d) * √(Ξt / (Σt * Δt)),
        √(Δt / Σt) * vr_d,
        csc(θ) * √(Σt / Ξt) * vϕ_d,
        -vθ_d / √Σt,
    ]
end

"""
Returns the tetrad that takes the ZAMO covector to the Boyer-Lindquist basis.

The tetrad is:

```
    √((ΣΔ)/A)   0        -2ar sin(θ)/√(AΣ)   0
    0           √(Σ/Δ)    0                  0
    0           0         0                 -√Σ
    0           0         √(A/Σ)sin(θ)       0
```
"""
function jac_bl_d_zamo_u(metric::Kerr{T}, r, θ) where {T}
    a = metric.spin
    Σt = Σ(metric, r, θ)
    Δt = Δ(metric, r)
    At = A(metric, r, θ)

    return @SMatrix [
        # coords = {t, r, ϕ, θ}
        √((Σt * Δt) / At) zero(T) -2*a*r*sin(θ)/√(At * Σt) zero(T)
        zero(T) √(Σt / Δt) zero(T) zero(T)
        zero(T) zero(T) zero(T) -√Σt
        zero(T) zero(T) √(At / Σt)*sin(θ) zero(T)
    ]
end

"""
Returns the ZAMO covector in the Boyer-Lindquist basis.

Applies the tetrad above without forming it.
"""
@inline function jac_bl_d_zamo_u(metric::Kerr{T}, r, θ, v_d) where {T}
    Δt = Δ(metric, r)
    Σt = Σ(metric, r, θ)
    Ξt = Ξ(metric, r, θ; Δ = Δt)
    ωt = T(2) * metric.spin * r / Ξt
    ϖ = sin(θ) * √(Ξt / Σt)

    vt_d, vr_d, vϕ_d, vθ_d = v_d
    return @SVector [
        √((Σt * Δt) / Ξt) * vt_d - ωt * ϖ * vϕ_d,
        √(Σt / Δt) * vr_d,
        -√Σt * vθ_d,
        ϖ * vϕ_d,
    ]
end

"""
Returns the tetrad that transforms the Boyer-Lindquist vector to the ZAMO basis.

The tetrad is (Eq 3.2 1972ApJ...178..347B):

```
    √((ΣΔ)/A)           0        0    0
    0                   √(Σ/Δ)   0    0
    -(2ar sin(θ))/√(AΣ) 0        0    √(A/Σ)sin(θ)
    0                   0       -√Σ   0
```
"""
function jac_zamo_u_bl_d(metric::Kerr{T}, r, θ) where {T}
    a = metric.spin
    Σt = Σ(metric, r, θ)
    Δt = Δ(metric, r)
    At = A(metric, r, θ)

    return @SMatrix [#  Eq 3.2 1972ApJ...178..347B
        # coords = {t, r, ϕ, θ}
        √((Σt * Δt) / At) zero(T) zero(T) zero(T)
        zero(T) √(Σt / Δt) zero(T) zero(T)
        -(2a * r * sin(θ))/√(At * Σt) zero(T) zero(T) √(At / Σt)*sin(θ)
        zero(T) zero(T) -√Σt zero(T)
    ]
end

"""
Transforms the Boyer-Lindquist vector to the ZAMO basis.

Applies the tetrad above without forming it.
"""
@inline function jac_zamo_u_bl_d(metric::Kerr{T}, r, θ, v_u) where {T}
    Δt = Δ(metric, r)
    Σt = Σ(metric, r, θ)
    Ξt = Ξ(metric, r, θ; Δ = Δt)
    ωt = T(2) * metric.spin * r / Ξt
    ϖ = sin(θ) * √(Ξt / Σt)

    vt_u, vr_u, vθ_u, vϕ_u = v_u
    return @SVector [
        √((Σt * Δt) / Ξt) * vt_u,
        √(Σt / Δt) * vr_u,
        ϖ * (vϕ_u - ωt * vt_u),
        -√Σt * vθ_u,
    ]
end

"""
Returns the tetrad that takes the ZAMO vector to the Boyer-Lindquist basis.

The tetrad is:

```
    √(A/(ΣΔ))    0        0             0
    0            √(Δ/Σ)   0             0
    0            0        0            -1/√Σ
    2ar/√(ΣΔA)   0        √(Σ/A)csc(θ)  0
```
"""
function jac_bl_u_zamo_d(metric::Kerr{T}, r, θ) where {T}
    a = metric.spin
    Σt = Σ(metric, r, θ)
    Δt = Δ(metric, r)
    At = A(metric, r, θ)

    return @SMatrix [
        # coords = {t, r, ϕ, θ}
        √(At / (Σt * Δt)) zero(T) zero(T) zero(T)
        zero(T) √(Δt / Σt) zero(T) zero(T)
        zero(T) zero(T) zero(T) -inv(√Σt)
        2a*r/√(Σt * Δt * At) zero(T) √(Σt / At)*csc(θ) zero(T)
    ]
end

"""
Returns the ZAMO vector in the Boyer-Lindquist basis.

Applies the tetrad above without forming it.
"""
@inline function jac_bl_u_zamo_d(metric::Kerr{T}, r, θ, v_u) where {T}
    Δt = Δ(metric, r)
    Σt = Σ(metric, r, θ)
    Ξt = Ξ(metric, r, θ; Δ = Δt)
    ωt = T(2) * metric.spin * r / Ξt

    vt_u, vr_u, vϕ_u, vθ_u = v_u
    vt_bl_u = √(Ξt / (Σt * Δt)) * vt_u
    return @SVector [
        vt_bl_u,
        √(Δt / Σt) * vr_u,
        -vθ_u / √Σt,
        ωt * vt_bl_u + csc(θ) * √(Σt / Ξt) * vϕ_u,
    ]
end

"""
Returns the Jacobian that takes the ZAMO vector into the frame of a fluid moving with
speed `β` along the direction `(θ, φ)`. Pass `-β` for the inverse boost.

The Jacobian is:

```
    γ            -βγcosφsinθ                   -βγsinφsinθ                -βγcosθ
    -βγcosφsinθ  cos²θcos²φ+γcos²φsin²θ+sin²φ  (γ-1)cosφsin²θsinφ         (γ-1)cosθcosφsinθ
    -βγsinθsinφ  (γ-1)cosφsin²θsinφ            cos²φ+(cos²θ+γsin²θ)sin²φ  (γ-1)cosθsinθsinφ
    -βγcosθ      (γ-1)cosθcosφsinθ             (γ-1)cosθsinθsinφ          γcos²θ+sin²θ
```
"""
function jac_fluid_u_zamo_d(::Kerr{T}, β, θ, φ) where {T}
    γ = inv(√(1 - β^2))
    sinφ = sin(φ)
    cosφ = cos(φ)
    sinθ = sin(θ)
    cosθ = cos(θ)

    return @SMatrix [
        γ -β*γ*cosφ*sinθ -β*γ*sinφ*sinθ -β*γ*cosθ
        -β*γ*cosφ*sinθ cosθ^2*cosφ^2+γ*cosφ^2*sinθ^2+sinφ^2 (γ-1)*cosφ*sinθ^2*sinφ (γ-1)*cosθ*cosφ*sinθ
        -β*γ*sinθ*sinφ (γ-1)*cosφ*sinθ^2*sinφ cosφ^2+(cosθ^2+γ*sinθ^2)*sinφ^2 (γ-1)*cosθ*sinθ*sinφ
        -β*γ*cosθ (γ-1)*cosθ*cosφ*sinθ (γ-1)*cosθ*sinθ*sinφ γ*cosθ^2+sinθ^2
    ]
end

"""
Returns the ZAMO vector in the frame of a fluid moving with speed `β` along the direction
`(θ, φ)`. Pass `-β` for the inverse boost.

Applies the Jacobian above without forming it. Its spatial block is `δᵢⱼ + (γ-1)n̂ᵢn̂ⱼ` for
`n̂ = (sin(θ)cos(φ), sin(θ)sin(φ), cos(θ))`, so the product is the boost
`γ(vᵗ - β n̂⋅v⃗)`, `v⃗ + ((γ-1)n̂⋅v⃗ - γβvᵗ)n̂`.
"""
@inline function jac_fluid_u_zamo_d(::Kerr{T}, β, θ, φ, v_u) where {T}
    γ = inv(√(one(T) - β^2))
    sinθ = sin(θ)
    nx = sinθ * cos(φ)
    ny = sinθ * sin(φ)
    nz = cos(θ)

    vt_u, vx_u, vy_u, vz_u = v_u
    ndotv = nx * vx_u + ny * vy_u + nz * vz_u
    scale = (γ - one(T)) * ndotv - γ * β * vt_u

    return @SVector [
        γ * (vt_u - β * ndotv),
        vx_u + scale * nx,
        vy_u + scale * ny,
        vz_u + scale * nz,
    ]
end

##----------------------------------------------------------------------------------------------------------------------
# Composed transformations
##----------------------------------------------------------------------------------------------------------------------
"""
Returns the photon momentum in the ZAMO basis.

Composes `p_bl_d`, `metric_uu` and `jac_zamo_u_bl_d`. The frame-dragging term the metric
puts into `p^ϕ` is the one the tetrad subtracts off again, so the composition is:

```
    p^t̂ = (1 - ωλ)/α
    p^r̂ = ±√R/√(ΣΔ)
    p^ϕ̂ = λ/ϖ
    p^θ̂ = ∓√Θ/√Σ
```

with lapse `α = √(ΣΔ/A)`, frame-dragging rate `ω = 2ar/A` and cylindrical radius
`ϖ = sin(θ)√(A/Σ)`.
"""
@inline function p_zamo_u(metric::Kerr{T}, r, θ, η, λ, νr::Bool, νθ::Bool) where {T}
    a = metric.spin
    Δt = Δ(metric, r)
    Σt = Σ(metric, r, θ)
    Ξt = Ξ(metric, r, θ; Δ = Δt)
    ωt = T(2) * a * r / Ξt

    # `p_bl_d` carries a 1/Δ on the radial component that the tetrad's √(Δ/Σ) undoes,
    # so take the square roots of the potentials undivided.
    pr = (νr ? one(T) : -one(T)) * √max(zero(T), r_potential(metric, η, λ, r))
    pθ = (νθ ? one(T) : -one(T)) * √max(zero(T), θ_potential(metric, η, λ, θ))

    return @SVector [
        (one(T) - ωt * λ) * √(Ξt / (Σt * Δt)),
        pr / √(Σt * Δt),
        λ * csc(θ) * √(Σt / Ξt),
        -pθ / √Σt,
    ]
end

"""
Returns the photon momentum and synchrotron emission vector in the fluid frame, shared by
[`Krang.synchrotronIntensity`](@ref) and [`Krang.synchrotronPolarization`](@ref).

Returns `(p_zamo_u, p_fluid_u, vec, norm)`, where `vec = p̂_fluid × B` and
`norm = |vec| + eps(T)`. 
"""
@inline function synchrotron_emission_frame(
    metric::Kerr{T},
    α,
    β,
    ri,
    θs,
    θo,
    magnetic_field::SVector{3,T},
    βfluid::SVector{3,T},
    νr::Bool,
    θsign::Bool,
) where {T}
    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)

    # qualified so the call resolves to the function rather than the local being assigned
    p_zamo_u = Krang.p_zamo_u(metric, ri, θs, ηtemp, λtemp, νr, θsign)
    p_fluid_u = jac_fluid_u_zamo_d(metric, βfluid[1], βfluid[2], βfluid[3], p_zamo_u)

    magnetic_fieldx, magnetic_fieldy, magnetic_fieldz = magnetic_field
    # the emission vector only cares about the photon's direction in the fluid frame
    _, p_fluid_ux, p_fluid_uy, p_fluid_uz = p_fluid_u ./ p_fluid_u[1]
    vec = @SVector [
        p_fluid_uy * magnetic_fieldz - p_fluid_uz * magnetic_fieldy,
        p_fluid_uz * magnetic_fieldx - p_fluid_ux * magnetic_fieldz,
        p_fluid_ux * magnetic_fieldy - p_fluid_uy * magnetic_fieldx,
    ]

    return p_zamo_u, p_fluid_u, vec, √sum(vec .* vec) + eps(T)
end

"""
Returns the Penrose walker constant for a photon with momentum p_u emitted from a fluid particle with momentum f_u.
"""
function penrose_walker(
    metric::Kerr{T},
    r,
    θ,
    p_u::AbstractVector,
    f_u::AbstractVector,
) where {T}# Eq 6 arXiv:2001.08750v1
    a = metric.spin
    pt, pr, pϕ, pθ = p_u
    ft, fr, fϕ, fθ = f_u
    sinθ = sin(θ)
    cosθ = cos(θ)

    A = pt * fr - pr * ft + a * sinθ * sinθ * (pr * fϕ - pϕ * fr)
    B = ((r * r + a * a) * (pϕ * fθ - pθ * fϕ) - a * (pt * fθ - pθ * ft)) * sinθ
    return A * r - B * a * cosθ, -(A * a * cosθ + B * r)
end
