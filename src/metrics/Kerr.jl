"""
Kerr Metric in Boyer Lindquist Coordinates
"""
struct Kerr{T} <: AbstractMetric
    mass::T
    spin::T
    function Kerr(spin::T) where T
      new{T}(one(T), spin)
    end
end

"""
    horizon(metric::Kerr{T})

Params

`metric` : Kerr{T} metric

Outer Horizon for the Kerr{T} metric.
"""
function horizon(metric::Kerr{T}) where T 
  m = metric.mass
  return m + √(m - metric.spin^2) 
end

function Δ(metric::Kerr{T}, r) where T return r^2 - T(2)r + metric.spin^2 end
function Σ(metric::Kerr{T}, r, θ) where T return r^2 + metric.spin^2 * cos(θ)^2 end
function A(metric::Kerr{T}, r, θ) where T 
  a = metric.spin
  return (r^2 + a^2)^2 - a^2 * Δ(metric, r) * sin(θ)^2 
end
function Ξ(metric::Kerr{T}, r, θ) where T 
  a = metric.spin
  return (r^2 + a^2)^2 - Δ(metric, r) * a^2 * sin(θ)^2 
end
function ω(metric::Kerr{T}, r, θ) where T return T(2) * metric.spin * r / Ξ(metric, r, θ) end


"""
  metric_uu(metric::Kerr{T}, r, θ)

# Params

`metric` : Kerr{T} metric

`r` : Radius

`θ` : Inclination 

Inverse Kerr{T} metric in Boyer Lindquist (BL) coordinates.
"""
function metric_uu(metric::Kerr{T}, r, θ) where T
  Ξt = Ξ(metric, r, θ)
  Σt = Σ(metric, r, θ)
  Δt = Δ(metric, r)
  ωt = ω(metric, r, θ)
  z = zero(T)

  return @SMatrix [ #Eq 1 2105.09440
    -Ξt/(Σt*Δt)     z       z                               -Ξt*ωt/(Σt*Δt)                  
    z               Δt/Σt   z                               z
    z               z       inv(Σt)                         z
    -Ξt*ωt/(Σt*Δt)  z       z                               Σt*csc(θ)^2/Ξt-Ξt*ωt^2/(Σt*Δt)
  ]
end

"""
  metric_uu(metric::Kerr{T}, coordinates)

# Params

`metric` : Kerr{T} metric

`coordinates` : Coordinates (t, r, θ, ϕ)

Inverse Kerr{T} metric in Boyer Lindquist (BL) coordinates.
"""
function metric_uu(metric::Kerr{T}, coordinates) where T
  _,r,θ,_ = coordinates
  metric_uu(metric, r, θ)
end

"""
    metric_dd(metric::Kerr{T}, r, θ)

# Params

`metric` : Kerr{T} metric

`r` : Radius

`θ` : Inclination 

Inverse Kerr{T} metric in Boyer Lindquist (BL) coordinates.
"""
function metric_dd(metric::Kerr{T}, r, θ) where T
  Ξt = Ξ(metric, r, θ)
  Σt = Σ(metric, r, θ)
  Δt = Δ(metric, r)
  ωt = ω(metric, r, θ)
  sin2 = sin(θ)^2
  z = zero(T)

  return @SMatrix [ #Eq 1 2105.09440
    -(Δt*Σt/Ξt)+(Ξt*ωt^2*sin2/Σt)     z         z                  -Ξt*ωt*sin2/Σt      
    z                                 Σt/Δt     z                  z
    z                                 z         Σt                 z
    -Ξt*ωt*sin2/Σt                    z         z                  Ξt*sin2/Σt
  ]
end


"""
  metric_dd(metric::Kerr{T}, coordinates)

# Params

`metric` : Kerr{T} metric

`coordinates` : Coordinates (t, r, θ, ϕ)

Kerr{T} metric in Boyer Lindquist (BL) coordinates.
"""
function metric_dd(metric::Kerr{T}, coordinates) where T
  _, r, θ, _ = coordinates
  return metric_dd(metric, r, θ)
end

