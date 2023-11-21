# Krang api

```@index
Pages = ["api.md"]
```

### Raytracing Functions
```@docs
Krang.emission_radius(metric::Kerr{T}, α, β, θs, θo, isindir, n) where T
Krang.emission_radius(metric::Kerr{T}, α, β, τ, θo) where T
Krang.emission_inclination(metric::Kerr{T}, α, β, θo, rs, νr) where T
Krang.emission_inclination(metric::Kerr{T}, α, β, τ, θo) where T
Krang.emission_coordinates_fast_light
Krang.emission_coordinates
Krang.raytrace
Krang.mino_time
```

### Metric Functions
```@docs
Krang.AbstractMetric
Krang.Kerr
Krang.metric_uu
Krang.metric_dd
Krang.horizon
```

```@docs
Krang.λ
Krang.η
```

### Radial Integrals
```@docs
Krang.r_potential
Krang.get_radial_roots
Krang.Ir(metric::Kerr, νr::Bool, θo, rs, α, β)
Krang.Ir(metric::Kerr{T}, νr::Bool, rs, η, λ) where T
Krang.radial_integrals_case2
Krang.radial_integrals_case3
Krang.radial_integrals_case4
```

### Angular Integrals
```@docs
Krang.θ_potential
Krang.Gθ
```

### Screen Coordinates
```@docs
Krang.α
Krang.β
Krang.αboundary
Krang.βboundary
```

### Polarization Functions
```@docs
Krang.p_bl_d
Krang.jac_bl_u_zamo_d
Krang.jac_zamo_u_bl_d
Krang.jac_bl_d_zamo_u
Krang.jac_zamo_d_bl_u
Krang.jac_fluid_u_zamo_d
Krang.screen_polarisation
Krang.penrose_walker
Krang.calcPol
```

### Misc
```@docs
Krang._isreal2
Krang.regularized_Pi
```