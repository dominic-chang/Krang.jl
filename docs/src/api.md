# Kang api

```@index
Pages = ["api.md"]
```

### Raytracing Functions
```@docs
Kang.emission_radius(metric::Kerr{T}, α, β, θs, θo, isindir, n) where T
Kang.emission_radius(metric::Kerr{T}, α, β, τ, θo) where T
Kang.emission_inclination(metric::Kerr{T}, α, β, θo, rs, νr) where T
Kang.emission_inclination(metric::Kerr{T}, α, β, τ, θo) where T
Kang.emission_coordinates_fast_light
Kang.emission_coordinates
Kang.raytrace
Kang.mino_time
```

### Metric Functions
```@docs
Kang.AbstractMetric
Kang.Kerr
Kang.metric_uu
Kang.metric_dd
Kang.horizon
```

```@docs
Kang.λ
Kang.η
```

### Radial Integrals
```@docs
Kang.r_potential
Kang.get_radial_roots
Kang.Ir(metric::Kerr, νr::Bool, θo, rs, α, β)
Kang.Ir(metric::Kerr{T}, νr::Bool, rs, η, λ) where T
Kang.radial_integrals_case2
Kang.radial_integrals_case3
Kang.radial_integrals_case4
```

### Angular Integrals
```@docs
Kang.θ_potential
Kang.Gθ
```

### Screen Coordinates
```@docs
Kang.α
Kang.β
Kang.αboundary
Kang.βboundary
```

### Polarization Functions
```@docs
Kang.p_bl_d
Kang.jac_bl_u_zamo_d
Kang.jac_zamo_u_bl_d
Kang.jac_bl_d_zamo_u
Kang.jac_zamo_d_bl_u
Kang.jac_fluid_u_zamo_d
Kang.screen_polarisation
Kang.penrose_walker
Kang.calcPol
```

### Misc
```@docs
Kang._isreal2
Kang.regularized_Pi
```