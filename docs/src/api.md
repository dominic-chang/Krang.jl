# Krang api

```@index
Pages = ["api.md"]
```

### Raytracing Functions
```@docs
Krang.emission_radius
Krang.emission_inclination
Krang.emission_coordinates_fast_light
Krang.emission_coordinates
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
Krang.Ir
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
Krang.screen_polarization
Krang.penrose_walker
Krang.synchrotronIntensity
Krang.synchrotronPolarization
```

### Raytracing API Related Functions
```@docs
Krang.Mesh
Krang.ElectronSynchrotronPowerLawPolarization
Krang.UnionGeometry
```

### Misc
```@docs
Krang._isreal2
Krang.regularized_Pi
Krang.radial_inf_integrals_case4 
Krang.It_w_I0_terms 
Krang.Ir_inf 
Krang.emission_azimuth 
Krang.ConeGeometry
Krang.Iϕ 
Krang.SlowLightIntensityCamera
Krang.Iϕ_w_I0_terms 
Krang.SlowLightIntensityPixel
Krang.Iϕ_inf 
Krang.radial_inf_integrals_case3 
Krang.AbstractPixel
Krang.radial_w_I0_terms_integrals_case2 
Krang.SlowLightIntensityScreen
Krang.radial_w_I0_terms_integrals_case4 
Krang.radial_integrals 
Krang.AbstractCamera
Krang.IntensityScreen
Krang.It 
Krang.IntensityPixel
Krang.AbstractMaterial
Krang.AbstractScreen
Krang.radial_w_I0_terms_integrals_case3 
Krang.IntensityCamera
Krang.Ir_s 
Krang.AbstractGeometry
Krang.radial_inf_integrals_case2 
Krang.It_inf 
```