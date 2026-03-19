export Kerr, horizon, metric_dd, metric_uu

"""
    $TYPEDEF

Kerr Metric in Boyer Lindquist Coordinates

$FIELDS
"""
struct Kerr{T} <: AbstractMetric
    "M = mass"
    mass::T
    "a = J/M, where J is the angular momentum and M is the mass of the black hole."
    spin::T

    @doc """
        Kerr(spin::T) where {T}

    Constructs a `Kerr` object representing the Kerr metric.

    # Arguments
    - `spin::T`: The spin parameter `a = J/M`, where `J` is the angular momentum and `M` is the mass of the black hole. spin ∈ (-1, 0) ∪ (0, 1).

    # Returns
    - A `Kerr` object with the given spin and a default mass of 1.
    """
    function Kerr(spin::T) where {T}
        new{T}(one(T), spin)
    end
end

"""
Outer Horizon for the Kerr metric.
"""
function horizon(metric::Kerr{T}) where {T}
    spin = metric.spin
    return 1 + √(1 - spin^2)
end

function Δ(metric::Kerr{T}, r) where {T}
    return r^2 - T(2)r + metric.spin^2
end
function Σ(metric::Kerr{T}, r, θ) where {T}
    return r^2 + metric.spin^2 * cos(θ)^2
end
function A(metric::Kerr{T}, r, θ; Δ = Δ(metric, r)) where {T}
    a = metric.spin
    return (r^2 + a^2)^2 - a^2 * Δ * sin(θ)^2
end
function Ξ(metric::Kerr{T}, r, θ; Δ = Δ(metric, r)) where {T}
    a = metric.spin
    return (r^2 + a^2)^2 - Δ * a^2 * sin(θ)^2
end
function ω(metric::Kerr{T}, r, θ; Ξ = Ξ(metric, r, θ)) where {T}
    return T(2) * metric.spin * r / Ξ
end


"""
Inverse Kerr metric in Boyer Lindquist (BL) coordinates.
"""
@inline function metric_uu(metric::Kerr{T}, r, θ) where {T}
    Δt = Δ(metric, r)
    Ξt = Ξ(metric, r, θ; Δ = Δt)
    Σt = Σ(metric, r, θ)
    ωt = ω(metric, r, θ; Ξ = Ξt)
    z = zero(T)

    return @SMatrix [ #Eq 1 2105.09440
        -Ξt/(Σt*Δt) z z -Ξt*ωt/(Σt*Δt)
        z Δt/Σt z z
        z z inv(Σt) z
        -Ξt*ωt/(Σt*Δt) z z Σt*csc(θ)^2/Ξt-Ξt*ωt^2/(Σt*Δt)
    ]
end

"""
Inverse Kerr metric in Boyer Lindquist (BL) coordinates.

# Arguments

- `metric` : Kerr metric
- `coordinates` : Coordinates (t, r, θ, ϕ)
"""
@inline function metric_uu(metric::Kerr{T}, coordinates) where {T}
    _, r, θ, _ = coordinates
    metric_uu(metric, r, θ)
end

"""
Kerr metric in Boyer Lindquist (BL) coordinates.
"""
@inline function metric_dd(metric::Kerr{T}, r, θ) where {T}
    Δt = Δ(metric, r)
    Ξt = Ξ(metric, r, θ; Δ = Δt)
    Σt = Σ(metric, r, θ)
    ωt = ω(metric, r, θ; Ξ = Ξt)
    sin2 = sin(θ)^2
    z = zero(T)

    return @SMatrix [ #Eq 1 2105.09440
        -(Δt * Σt / Ξt)+(Ξt*ωt^2*sin2/Σt) z z -Ξt*ωt*sin2/Σt
        z Σt/Δt z z
        z z Σt z
        -Ξt*ωt*sin2/Σt z z Ξt*sin2/Σt
    ]
end


"""
Kerr metric in Boyer Lindquist (BL) coordinates.

# Arguments

- `metric` : Kerr metric
- `coordinates` : Coordinates (t, r, θ, ϕ)
"""
@inline function metric_dd(metric::Kerr{T}, coordinates) where {T}
    _, r, θ, _ = coordinates
    return metric_dd(metric, r, θ)
end
