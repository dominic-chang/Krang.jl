
# Krang api {#Krang-api}
- [`Krang.AbstractCamera`](#Krang.AbstractCamera)
- [`Krang.AbstractGeometry`](#Krang.AbstractGeometry)
- [`Krang.AbstractMaterial`](#Krang.AbstractMaterial)
- [`Krang.AbstractMetric`](#Krang.AbstractMetric)
- [`Krang.AbstractPixel`](#Krang.AbstractPixel)
- [`Krang.AbstractScreen`](#Krang.AbstractScreen)
- [`Krang.ConeGeometry`](#Krang.ConeGeometry)
- [`Krang.ElectronSynchrotronPowerLawPolarization`](#Krang.ElectronSynchrotronPowerLawPolarization)
- [`Krang.IntensityCamera`](#Krang.IntensityCamera)
- [`Krang.IntensityPixel`](#Krang.IntensityPixel)
- [`Krang.IntensityScreen`](#Krang.IntensityScreen)
- [`Krang.Kerr`](#Krang.Kerr)
- [`Krang.Mesh`](#Krang.Mesh)
- [`Krang.SlowLightIntensityCamera`](#Krang.SlowLightIntensityCamera)
- [`Krang.SlowLightIntensityPixel`](#Krang.SlowLightIntensityPixel)
- [`Krang.SlowLightIntensityScreen`](#Krang.SlowLightIntensityScreen)
- [`Krang.UnionGeometry`](#Krang.UnionGeometry)
- [`Krang.Gθ`](#Krang.Gθ)
- [`Krang.Ir`](#Krang.Ir)
- [`Krang.Ir_inf`](#Krang.Ir_inf)
- [`Krang.Ir_s`](#Krang.Ir_s)
- [`Krang.It`](#Krang.It)
- [`Krang.It_inf`](#Krang.It_inf)
- [`Krang.It_m_I0_terms`](#Krang.It_m_I0_terms)
- [`Krang.It_w_I0_terms`](#Krang.It_w_I0_terms)
- [`Krang.Iϕ`](#Krang.Iϕ)
- [`Krang.Iϕ_inf`](#Krang.Iϕ_inf)
- [`Krang.Iϕ_m_I0_terms`](#Krang.Iϕ_m_I0_terms)
- [`Krang.Iϕ_w_I0_terms`](#Krang.Iϕ_w_I0_terms)
- [`Krang._isreal2`](#Krang._isreal2)
- [`Krang.emission_azimuth`](#Krang.emission_azimuth)
- [`Krang.emission_coordinates`](#Krang.emission_coordinates)
- [`Krang.emission_coordinates_fast_light`](#Krang.emission_coordinates_fast_light)
- [`Krang.emission_inclination`](#Krang.emission_inclination)
- [`Krang.emission_radius`](#Krang.emission_radius)
- [`Krang.get_radial_roots`](#Krang.get_radial_roots)
- [`Krang.horizon`](#Krang.horizon)
- [`Krang.jac_bl_d_zamo_u`](#Krang.jac_bl_d_zamo_u)
- [`Krang.jac_bl_u_zamo_d`](#Krang.jac_bl_u_zamo_d)
- [`Krang.jac_fluid_u_zamo_d`](#Krang.jac_fluid_u_zamo_d)
- [`Krang.jac_zamo_d_bl_u`](#Krang.jac_zamo_d_bl_u)
- [`Krang.jac_zamo_u_bl_d`](#Krang.jac_zamo_u_bl_d)
- [`Krang.metric_dd`](#Krang.metric_dd)
- [`Krang.metric_uu`](#Krang.metric_uu)
- [`Krang.mino_time`](#Krang.mino_time)
- [`Krang.p_bl_d`](#Krang.p_bl_d)
- [`Krang.penrose_walker`](#Krang.penrose_walker)
- [`Krang.r_potential`](#Krang.r_potential)
- [`Krang.radial_inf_integrals_case2`](#Krang.radial_inf_integrals_case2)
- [`Krang.radial_inf_integrals_case3`](#Krang.radial_inf_integrals_case3)
- [`Krang.radial_inf_integrals_case4`](#Krang.radial_inf_integrals_case4)
- [`Krang.radial_integrals`](#Krang.radial_integrals)
- [`Krang.radial_m_I0_terms_integrals_case2`](#Krang.radial_m_I0_terms_integrals_case2)
- [`Krang.radial_m_I0_terms_integrals_case3`](#Krang.radial_m_I0_terms_integrals_case3)
- [`Krang.radial_m_I0_terms_integrals_case4`](#Krang.radial_m_I0_terms_integrals_case4)
- [`Krang.radial_w_I0_terms_integrals_case2`](#Krang.radial_w_I0_terms_integrals_case2)
- [`Krang.radial_w_I0_terms_integrals_case3`](#Krang.radial_w_I0_terms_integrals_case3)
- [`Krang.radial_w_I0_terms_integrals_case4`](#Krang.radial_w_I0_terms_integrals_case4)
- [`Krang.raytrace`](#Krang.raytrace)
- [`Krang.regularized_Pi`](#Krang.regularized_Pi)
- [`Krang.screen_polarisation`](#Krang.screen_polarisation)
- [`Krang.synchrotronIntensity`](#Krang.synchrotronIntensity)
- [`Krang.synchrotronPolarization`](#Krang.synchrotronPolarization)
- [`Krang.α`](#Krang.α)
- [`Krang.αboundary`](#Krang.αboundary)
- [`Krang.β`](#Krang.β)
- [`Krang.βboundary`](#Krang.βboundary)
- [`Krang.η`](#Krang.η)
- [`Krang.θ_potential`](#Krang.θ_potential)
- [`Krang.λ`](#Krang.λ)


### Raytracing Functions {#Raytracing-Functions}
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.emission_radius' href='#Krang.emission_radius'>#</a>&nbsp;<b><u>Krang.emission_radius</u></b> &mdash; <i>Function</i>.




```julia
emission_radius(
    pix::Krang.AbstractPixel,
    θs,
    isindir,
    n
) -> Tuple

```


Emission radius for point originating at inclination θs whose nth order image appears at the screen coordinate (`α`, `β`).  Returns NaN if the emission coordinates do not exist for that screen coordinate.

**Arguments**
- `pix` : Pixel information 
  
- `θs` : Emission inclination
  
- `isindir` : Is emission to observer direct or indirect
  
- `n` : Image index
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/raytracer.jl#L5)



```julia
emission_radius(
    pix::Krang.AbstractPixel,
    τ
) -> Tuple{Any, Any, Any}

```


Emission radius for point originating at at Mino time τ whose image appears at the screen coordinate (`α`, `β`).  Returns NaN if the emission coordinates do not exist for that screen coordinate.

**Arguments**

-`pix` : Pixel information
- `τ` : Mino time
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/raytracer.jl#L41)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.emission_inclination' href='#Krang.emission_inclination'>#</a>&nbsp;<b><u>Krang.emission_inclination</u></b> &mdash; <i>Function</i>.




```julia
emission_inclination(pix::Krang.AbstractPixel, rs, νr)

```


Emission inclination for point originating at inclination rs whose nth order image appears at screen coordinate (`α`, `β`).

**Arguments**
- `pix` : Pixel information
  
- `rs` : Emission radius
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/raytracer.jl#L70)



```julia
emission_inclination(
    pix::Krang.AbstractPixel,
    τ
) -> NTuple{6, Any}

```


Emission inclination for point at Mino time τ whose image appears at screen coordinate (`α`, `β`).

**Arguments**
- `pix` : Pixel information
  
- `τ` : Mino Time
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/raytracer.jl#L83)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.emission_coordinates_fast_light' href='#Krang.emission_coordinates_fast_light'>#</a>&nbsp;<b><u>Krang.emission_coordinates_fast_light</u></b> &mdash; <i>Function</i>.




```julia
emission_coordinates_fast_light(
    pix::Krang.AbstractPixel,
    θs,
    isindir,
    n
) -> NTuple{5, Any}

```


Emission radius and azimuthal angle for point originating at inclination θs whose nth order image appears at the screen coordinate (`α`, `β`).  for an observer located at inclination θo.

**Arguments**
- `pix` : Pixel information
  
- `θs` : Emission Inclination
  
- `isindir` : Whether emission to observer is direct or indirect
  
- `n` : Image index
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/raytracer.jl#L125)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.emission_coordinates' href='#Krang.emission_coordinates'>#</a>&nbsp;<b><u>Krang.emission_coordinates</u></b> &mdash; <i>Function</i>.




```julia
emission_coordinates(
    pix::Krang.AbstractPixel,
    θs,
    isindir,
    n
) -> NTuple{6, Any}

```


Emission radius and azimuthal angle for point originating at inclination θs whose nth order image appears at the screen coordinate (`α`, `β`). for an observer located at inclination θo.

**Arguments**
- `pix` : Pixel information
  
- `θs` : Emission Inclination
  
- `isindir` : Whether emission to observer is direct or indirect
  
- `n` : Image index
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/raytracer.jl#L160)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.raytrace' href='#Krang.raytrace'>#</a>&nbsp;<b><u>Krang.raytrace</u></b> &mdash; <i>Function</i>.




```julia
raytrace(pix::Krang.AbstractPixel, τ) -> NTuple{6, Any}

```


Raytrace a point that appears at the screen coordinate (`α`, `β`) for an observer located at inclination θo

**Arguments**
- `pix` : Pixel information
  
- `τ` : Mino Time
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/raytracer.jl#L221)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.mino_time' href='#Krang.mino_time'>#</a>&nbsp;<b><u>Krang.mino_time</u></b> &mdash; <i>Function</i>.




```julia
mino_time(pix, rs, isindir) -> Any

```


Mino time of trajectory between an observer at infinity and point at radius rs

**Arguments**
- `pix` : Pixel information
  
- `rs` : Emission radius
  
- `isindir` : Is the path direct or indirect?
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L283)



```julia
mino_time(pix::Krang.AbstractPixel, θs, isindir, n) -> Any

```


Mino time of trajectory between two inclinations for a given screen coordinate

**Arguments**
- `pix` : Pixel information
  
- `θs` : Emission inclination
  
- `isindir` : Is the path direct or indirect?
  
- `n` : nth image in orde of amount of minotime traversed
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1906)

</div>
<br>

### Metric Functions {#Metric-Functions}
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.AbstractMetric' href='#Krang.AbstractMetric'>#</a>&nbsp;<b><u>Krang.AbstractMetric</u></b> &mdash; <i>Type</i>.




```julia
abstract type AbstractMetric
```


Abstract Metric Type


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/AbstractMetric.jl#L1)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Kerr' href='#Krang.Kerr'>#</a>&nbsp;<b><u>Krang.Kerr</u></b> &mdash; <i>Type</i>.




```julia
struct Kerr{T} <: Krang.AbstractMetric
```


Kerr Metric in Boyer Lindquist Coordinates
- `mass`: M = mass
  
- `spin`: a = J/M, where J is the angular momentum and M is the mass of the blackhole.
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/Kerr.jl#L3)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.metric_uu' href='#Krang.metric_uu'>#</a>&nbsp;<b><u>Krang.metric_uu</u></b> &mdash; <i>Function</i>.




```julia
metric_uu(metric::Krang.AbstractMetric, args...) -> Any

```


Returns the inverse metric in some representation (usually as an nxn matrix).


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/AbstractMetric.jl#L15)



```julia
metric_uu(metric::Kerr{T}, r, θ) -> Any

```


Inverse Kerr metric in Boyer Lindquist (BL) coordinates.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/Kerr.jl#L47)



```julia
metric_uu(metric::Kerr{T}, coordinates) -> Any

```


Inverse Kerr metric in Boyer Lindquist (BL) coordinates.

**Arguments**
- `metric` : Kerr metric
  
- `coordinates` : Coordinates (t, r, θ, ϕ)
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/Kerr.jl#L65)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.metric_dd' href='#Krang.metric_dd'>#</a>&nbsp;<b><u>Krang.metric_dd</u></b> &mdash; <i>Function</i>.




```julia
metric_dd(metric::Krang.AbstractMetric, args...) -> Any

```


Returns the metric in some representation (usually as an nxn matrix).


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/AbstractMetric.jl#L8)



```julia
metric_dd(metric::Kerr{T}, r, θ) -> Any

```


Inverse Kerr metric in Boyer Lindquist (BL) coordinates.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/Kerr.jl#L78)



```julia
metric_dd(metric::Kerr{T}, coordinates) -> Any

```


Kerr metric in Boyer Lindquist (BL) coordinates.

**Arguments**
- `metric` : Kerr metric
  
- `coordinates` : Coordinates (t, r, θ, ϕ)
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/Kerr.jl#L98)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.horizon' href='#Krang.horizon'>#</a>&nbsp;<b><u>Krang.horizon</u></b> &mdash; <i>Function</i>.




```julia
horizon(metric::Kerr{T}) -> Any

```


Outer Horizon for the Kerr metric.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/Kerr.jl#L20)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.λ' href='#Krang.λ'>#</a>&nbsp;<b><u>Krang.λ</u></b> &mdash; <i>Function</i>.




```julia
λ(_::Kerr, α, θo) -> Any

```


Energy reduced azimuthal angular momentum

**Arguments**
- `metric`: Kerr
  
- `α`: Horizontal Bardeen screen coordinate
  
- `θo`: Observer inclination
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L112)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.η' href='#Krang.η'>#</a>&nbsp;<b><u>Krang.η</u></b> &mdash; <i>Function</i>.




```julia
η(metric::Kerr, α, β, θo) -> Any

```


Energy reduced Carter integral

**Arguments**
- `metric`: Kerr
  
- `α`: Horizontal Bardeen screen coordinate
  
- `β`: Bardeen vertical coordinate
  
- `θo`: Observer inclination
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L125)

</div>
<br>

### Radial Integrals {#Radial-Integrals}
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.r_potential' href='#Krang.r_potential'>#</a>&nbsp;<b><u>Krang.r_potential</u></b> &mdash; <i>Function</i>.




```julia
r_potential(metric::Kerr{T}, η, λ, r) -> Any

```


Radial potential of spacetime

**Arguments**
- `metric`: Kerr{T} metric
  
- `η`  : Reduced Carter constant
  
- `λ`  : Reduced azimuthal angular momentum
  
- `r`  : Boyer Lindquist radius
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L194)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.get_radial_roots' href='#Krang.get_radial_roots'>#</a>&nbsp;<b><u>Krang.get_radial_roots</u></b> &mdash; <i>Function</i>.




```julia
get_radial_roots(metric::Kerr{T}, η, λ) -> NTuple{4, Any}

```


Returns roots of $r^4 + (a^2-η-λ^2)r^2 + 2(η+(a-λ)^2)r - a^2η$

**Arguments**
- `metric`: Kerr{T} metric
  
- `η`  : Reduced Carter constant
  
- `λ`  : Reduced azimuthal angular momentum
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L228)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Ir' href='#Krang.Ir'>#</a>&nbsp;<b><u>Krang.Ir</u></b> &mdash; <i>Function</i>.




```julia
Ir(pix::Krang.AbstractPixel, νr::Bool, rs) -> Any

```


Returns the antiderivative $I_r=\int\frac{dr}{\sqrt{\mathcal{R(r)}}}$. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `pix`  : Pixel information
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  
- `rs` : Emission radius
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1754)

</div>
<br>

### Angular Integrals {#Angular-Integrals}
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.θ_potential' href='#Krang.θ_potential'>#</a>&nbsp;<b><u>Krang.θ_potential</u></b> &mdash; <i>Function</i>.




```julia
θ_potential(metric::Kerr{T}, η, λ, θ) -> Any

```


Theta potential of a Kerr blackhole

**Arguments**
- `metric`: Kerr{T} metric
  
- `η`  : Reduced Carter constant
  
- `λ`  : Reduced azimuthal angular momentum
  
- `θ`  : Boyer Lindquist inclination
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L210)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Gθ' href='#Krang.Gθ'>#</a>&nbsp;<b><u>Krang.Gθ</u></b> &mdash; <i>Function</i>.




```julia
Gθ(
    pix::Krang.AbstractPixel,
    θs,
    isindir,
    n
) -> NTuple{5, Any}

```


Returns the antiderivative $G_\theta=\int\frac{d\theta}{\sqrt{\Theta(\theta)}}$. See [`θ_potential(x)`](/api#Krang.θ_potential) for an implementation of $\Theta(	heta)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `pix` : Pixel information
  
- `θs` : Emission inclination
  
- `θo` : Observer inclination
  
- `isindir` : Is the path direct or indirect?
  
- `n` : nth image ordered by minotime
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1921)

</div>
<br>

### Screen Coordinates {#Screen-Coordinates}
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.α' href='#Krang.α'>#</a>&nbsp;<b><u>Krang.α</u></b> &mdash; <i>Function</i>.




```julia
α(_::Kerr, λ, θo) -> Any

```


Horizontal Bardeen Screen Coordinate

**Arguments**
- `metric`: Kerr
  
- `α`: Horizontal Bardeen screen coordinate
  
- `θo`: Observer inclination
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L139)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.β' href='#Krang.β'>#</a>&nbsp;<b><u>Krang.β</u></b> &mdash; <i>Function</i>.




```julia
β(metric::Kerr, λ, η, θo) -> Any

```


Horizontal Bardeen Screen Coordinate

**Arguments**
- `metric`: Kerr
  
- `λ`: Energy reduced Azimuthal angular momentul
  
- `η`: Energy reduced Carter integral 
  
- `θo`: Observer inclination
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L152)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.αboundary' href='#Krang.αboundary'>#</a>&nbsp;<b><u>Krang.αboundary</u></b> &mdash; <i>Function</i>.




```julia
αboundary(metric::Kerr, θs) -> Any

```


Defines a horizontal boundary on the assmyptotic observers screen where emission that originates from θs must fall within.

**Arguments**
- `metric`: Kerr metric
  
- `θs`  : Emission Inclination
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L166)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.βboundary' href='#Krang.βboundary'>#</a>&nbsp;<b><u>Krang.βboundary</u></b> &mdash; <i>Function</i>.




```julia
βboundary(metric::Kerr{T}, α, θo, θs) -> Any

```


Defines a vertical boundary on the Assyptotic observers screen where emission that originates from θs must fall within.

**Arguments**
- `metric`: Kerr{T} metric
  
- `α`   : Horizontal Bardeen screen coordinate
  
- `θo`  : Observer inclination
  
- `θs`  : Emission Inclination
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L178)

</div>
<br>

### Polarization Functions {#Polarization-Functions}
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.p_bl_d' href='#Krang.p_bl_d'>#</a>&nbsp;<b><u>Krang.p_bl_d</u></b> &mdash; <i>Function</i>.




```julia
p_bl_d(
    metric::Kerr{T},
    r,
    θ,
    η,
    λ,
    νr::Bool,
    νθ::Bool
) -> Any

```


```
Returns the momentum form in the Boyer-Lindquist basis.
```



[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/physicsUtils.jl#L5)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.jac_bl_u_zamo_d' href='#Krang.jac_bl_u_zamo_d'>#</a>&nbsp;<b><u>Krang.jac_bl_u_zamo_d</u></b> &mdash; <i>Function</i>.




```julia
jac_bl_u_zamo_d(metric::Kerr{T}, r, θ) -> Any

```


Jacobian which converts ZAMO vector to a Boyer-Lindquist basis


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/physicsUtils.jl#L71)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.jac_zamo_u_bl_d' href='#Krang.jac_zamo_u_bl_d'>#</a>&nbsp;<b><u>Krang.jac_zamo_u_bl_d</u></b> &mdash; <i>Function</i>.




```julia
jac_zamo_u_bl_d(metric::Kerr{T}, r, θ) -> Any

```


Jacobian which converts Boyer-Lindquist vector to a ZAMO basis


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/physicsUtils.jl#L53)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.jac_bl_d_zamo_u' href='#Krang.jac_bl_d_zamo_u'>#</a>&nbsp;<b><u>Krang.jac_bl_d_zamo_u</u></b> &mdash; <i>Function</i>.




```julia
jac_bl_d_zamo_u(metric::Kerr{T}, r, θ) -> Any

```


Jacobian which converts ZAMO covector to a Boyer-Lindquist basis


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/physicsUtils.jl#L35)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.jac_zamo_d_bl_u' href='#Krang.jac_zamo_d_bl_u'>#</a>&nbsp;<b><u>Krang.jac_zamo_d_bl_u</u></b> &mdash; <i>Function</i>.




```julia
jac_zamo_d_bl_u(metric::Kerr{T}, r, θ) -> Any

```


Returns the Jacobian which converts a Boyer-Lindquist covector to ZAMO basis.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/physicsUtils.jl#L17)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.jac_fluid_u_zamo_d' href='#Krang.jac_fluid_u_zamo_d'>#</a>&nbsp;<b><u>Krang.jac_fluid_u_zamo_d</u></b> &mdash; <i>Function</i>.




```julia
jac_fluid_u_zamo_d(_::Kerr{T}, β, θ, φ) -> Any

```


Jacobian which expreases ZAMO vector in the fluid frame


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/physicsUtils.jl#L89)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.screen_polarisation' href='#Krang.screen_polarisation'>#</a>&nbsp;<b><u>Krang.screen_polarisation</u></b> &mdash; <i>Function</i>.




```julia
screen_polarisation(
    metric::Kerr{T},
    κ::Complex,
    θ,
    α,
    β
) -> Tuple{Any, Any}

```


Returns the screen polarization associated with a killing spinor κ as seen seen by an assymptotic observer.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/ElectronSynchrotronPowerLawPolarization.jl#L1)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.penrose_walker' href='#Krang.penrose_walker'>#</a>&nbsp;<b><u>Krang.penrose_walker</u></b> &mdash; <i>Function</i>.




```julia
penrose_walker(
    metric::Kerr{T},
    r,
    θ,
    p_u::AbstractVector,
    f_u::AbstractVector
) -> Tuple{Any, Any}

```


Returns the Penrose walker constant for a photon with momentum p_u emitted from a fluid particle with momentum f_u.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/physicsUtils.jl#L107)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.synchrotronIntensity' href='#Krang.synchrotronIntensity'>#</a>&nbsp;<b><u>Krang.synchrotronIntensity</u></b> &mdash; <i>Function</i>.




```julia
synchrotronIntensity(
    metric::Kerr{T},
    α,
    β,
    ri,
    θs,
    θo,
    magfield::StaticArraysCore.SArray{Tuple{3}, T, 1, 3},
    βfluid::StaticArraysCore.SArray{Tuple{3}, T, 1, 3},
    νr::Bool,
    θsign::Bool
) -> Tuple{Any, Any, Any}

```


Calculates the intensity of a photon emitted from a fluid particle with momentum f_u and observed by an asymptotic observer.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/ElectronSynchrotronPowerLawIntensity.jl#L1)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.synchrotronPolarization' href='#Krang.synchrotronPolarization'>#</a>&nbsp;<b><u>Krang.synchrotronPolarization</u></b> &mdash; <i>Function</i>.




```julia
synchrotronPolarization(
    metric::Kerr{T},
    α,
    β,
    ri,
    θs,
    θo,
    magfield::StaticArraysCore.SArray{Tuple{3}, T, 1, 3},
    βfluid::StaticArraysCore.SArray{Tuple{3}, T, 1, 3},
    νr::Bool,
    θsign::Bool
) -> NTuple{4, Any}

```


Calculates the polarization of a photon emitted from a fluid particle with momentum f_u and observed by an asymptotic observer.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/ElectronSynchrotronPowerLawPolarization.jl#L19)

</div>
<br>

### Raytracing API Related Functions {#Raytracing-API-Related-Functions}
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Mesh' href='#Krang.Mesh'>#</a>&nbsp;<b><u>Krang.Mesh</u></b> &mdash; <i>Type</i>.




```julia
struct Mesh{G<:Krang.AbstractGeometry, M<:Krang.AbstractMaterial}
```



[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/geometries/geometry_types.jl#L46)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.ElectronSynchrotronPowerLawPolarization' href='#Krang.ElectronSynchrotronPowerLawPolarization'>#</a>&nbsp;<b><u>Krang.ElectronSynchrotronPowerLawPolarization</u></b> &mdash; <i>Type</i>.




```julia
struct ElectronSynchrotronPowerLawPolarization <: Krang.AbstractMaterial
```


Linear polarization material from https://doi.org/10.3847/1538-4357/abf117


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/ElectronSynchrotronPowerLawPolarization.jl#L68)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.UnionGeometry' href='#Krang.UnionGeometry'>#</a>&nbsp;<b><u>Krang.UnionGeometry</u></b> &mdash; <i>Type</i>.




```julia
struct UnionGeometry{G1, G2} <: Krang.AbstractGeometry
```


Geometry that is comprised of the union of two geometries.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/geometries/geometry_types.jl#L32)

</div>
<br>

### Misc {#Misc}
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang._isreal2' href='#Krang._isreal2'>#</a>&nbsp;<b><u>Krang._isreal2</u></b> &mdash; <i>Function</i>.




```julia
_isreal2(num::Complex{T}) -> Any

```


Checks if a complex number is real to some tolerance


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L20)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.regularized_Pi' href='#Krang.regularized_Pi'>#</a>&nbsp;<b><u>Krang.regularized_Pi</u></b> &mdash; <i>Function</i>.




```julia
regularized_Pi(n, ϕ, k) -> Any

```


Regularized elliptic integral of the third kind

**Arguments**
- `n`: Parameter
  
- `ϕ`: Arguments
  
- `k`: Parameter
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L30)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_m_I0_terms_integrals_case3' href='#Krang.radial_m_I0_terms_integrals_case3'>#</a>&nbsp;<b><u>Krang.radial_m_I0_terms_integrals_case3</u></b> &mdash; <i>Function</i>.




```julia
radial_m_I0_terms_integrals_case3(
    metric::Kerr{T},
    rs,
    roots::NTuple{4, T} where T
) -> Tuple

```


Returns the radial integrals for the case where there are two real roots in the radial potential


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1662)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_inf_integrals_case4' href='#Krang.radial_inf_integrals_case4'>#</a>&nbsp;<b><u>Krang.radial_inf_integrals_case4</u></b> &mdash; <i>Function</i>.




```julia
radial_inf_integrals_case4(
    metric::Kerr{T},
    roots::NTuple{4, T} where T
) -> NTuple{4, Any}

```


Returns the radial integrals for the case where there are no real roots in the radial potential


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1397)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.It_m_I0_terms' href='#Krang.It_m_I0_terms'>#</a>&nbsp;<b><u>Krang.It_m_I0_terms</u></b> &mdash; <i>Function</i>.




```julia
It_m_I0_terms(
    metric::Kerr{T},
    rs,
    roots::NTuple{4, T} where T,
    λ,
    νr
) -> Any

```


Returns the antiderivative $I_t=\int\frac{r^2\Delta+2Mr(r^2+a^2-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$  without I0 terms.

See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `roots` : Radial roots
  
- `λ`  : Reduced azimuthal angular momentum
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1148)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.It_w_I0_terms' href='#Krang.It_w_I0_terms'>#</a>&nbsp;<b><u>Krang.It_w_I0_terms</u></b> &mdash; <i>Function</i>.




```julia
It_w_I0_terms(
    metric::Kerr{T},
    rs,
    τ,
    roots::NTuple{4, T} where T,
    λ,
    νr
) -> Any

```


Returns the antiderivative $I_t=\int\frac{r^2\Delta+2Mr(r^2+a^2-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$ with  I0 terms.

See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `roots` : Radial roots
  
- `λ`  : Reduced azimuthal angular momentum
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L982)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Ir_inf' href='#Krang.Ir_inf'>#</a>&nbsp;<b><u>Krang.Ir_inf</u></b> &mdash; <i>Function</i>.




```julia
Ir_inf(metric::Kerr{T}, roots) -> Any

```


Returns the antiderivative $I_r=\int\frac{dr}{\sqrt{\mathcal{R(r)}}}$ evaluated at infinity. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `roots`  : Roots of the radial potential
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L296)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.emission_azimuth' href='#Krang.emission_azimuth'>#</a>&nbsp;<b><u>Krang.emission_azimuth</u></b> &mdash; <i>Function</i>.




```julia
emission_azimuth(
    pix::Krang.AbstractPixel,
    θs,
    rs,
    τ,
    νr,
    isindir,
    n
) -> Any

```


Emission azimuth for point at Mino time τ whose image appears at screen coordinate (`α`, `β`).

**Arguments**
- `pix` : Pixel information
  
- `θs` : Emission Inclination
  
- `rs` :Emission radius
  
- `τ` : Mino Time
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/raytracer.jl#L99)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.ConeGeometry' href='#Krang.ConeGeometry'>#</a>&nbsp;<b><u>Krang.ConeGeometry</u></b> &mdash; <i>Type</i>.




```julia
struct ConeGeometry{T, A} <: Krang.AbstractGeometry
```


Cone Geometry with half opening angle `opening_angle`.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/geometries/geometry_types.jl#L18)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Iϕ' href='#Krang.Iϕ'>#</a>&nbsp;<b><u>Krang.Iϕ</u></b> &mdash; <i>Function</i>.




```julia
Iϕ(pix::Krang.AbstractPixel, rs, τ, νr) -> Any

```


Returns the antiderivative $I_ϕ=\int\frac{a(2Mr-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `pix`: SlowLightIntensityPixel
  
- `rs` : Emission radius
  
- `τ` : Mino time
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1768)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.SlowLightIntensityCamera' href='#Krang.SlowLightIntensityCamera'>#</a>&nbsp;<b><u>Krang.SlowLightIntensityCamera</u></b> &mdash; <i>Type</i>.




```julia
struct SlowLightIntensityCamera{T} <: Krang.AbstractCamera
```


Observer sitting at radial infinity. The frame of this observer is alligned with the Boyer-Lindquist frame.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/cameras/SlowLightIntensityCamera.jl#L87)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Iϕ_w_I0_terms' href='#Krang.Iϕ_w_I0_terms'>#</a>&nbsp;<b><u>Krang.Iϕ_w_I0_terms</u></b> &mdash; <i>Function</i>.




```julia
Iϕ_w_I0_terms(metric::Kerr{T}, rs, τ, roots, νr, λ) -> Any

```


Returns the antiderivative $I_ϕ=\int\frac{a(2Mr-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$ with full I0 terms.

See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `τ`: Minotime 
  
- `roots` : Radial roots
  
- `λ`  : Reduced azimuthal angular momentum
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L699)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.SlowLightIntensityPixel' href='#Krang.SlowLightIntensityPixel'>#</a>&nbsp;<b><u>Krang.SlowLightIntensityPixel</u></b> &mdash; <i>Type</i>.




```julia
struct SlowLightIntensityPixel{T} <: Krang.AbstractPixel
```


Intensity Pixel Type.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/cameras/SlowLightIntensityCamera.jl#L3)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_m_I0_terms_integrals_case4' href='#Krang.radial_m_I0_terms_integrals_case4'>#</a>&nbsp;<b><u>Krang.radial_m_I0_terms_integrals_case4</u></b> &mdash; <i>Function</i>.




```julia
radial_m_I0_terms_integrals_case4(
    metric::Kerr{T},
    rs,
    roots::NTuple{4, T} where T
) -> Tuple

```


Returns the radial integrals for the case where there are no real roots in the radial potential


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1707)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Iϕ_inf' href='#Krang.Iϕ_inf'>#</a>&nbsp;<b><u>Krang.Iϕ_inf</u></b> &mdash; <i>Function</i>.




```julia
Iϕ_inf(metric::Kerr{T}, roots, λ) -> Any

```


Returns the antiderivative $I_ϕ=\int\frac{a(2Mr-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$ evaluated at infinity without I0 terms.

See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `roots` : Radial roots
  
- `λ`  : Reduced azimuthal angular momentum
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L444)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_m_I0_terms_integrals_case2' href='#Krang.radial_m_I0_terms_integrals_case2'>#</a>&nbsp;<b><u>Krang.radial_m_I0_terms_integrals_case2</u></b> &mdash; <i>Function</i>.




```julia
radial_m_I0_terms_integrals_case2(
    metric::Kerr{T},
    rs,
    roots::NTuple{4, T} where T,
    νr
) -> Tuple

```


Returns the radial integrals for the case where there are four real roots in the radial potential, with roots outside the horizon.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1612)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_inf_integrals_case3' href='#Krang.radial_inf_integrals_case3'>#</a>&nbsp;<b><u>Krang.radial_inf_integrals_case3</u></b> &mdash; <i>Function</i>.




```julia
radial_inf_integrals_case3(
    metric::Kerr{T},
    roots::NTuple{4, T} where T
) -> NTuple{4, Any}

```


Returns the radial integrals for the case where there are two real roots in the radial potential


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1359)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.AbstractPixel' href='#Krang.AbstractPixel'>#</a>&nbsp;<b><u>Krang.AbstractPixel</u></b> &mdash; <i>Type</i>.




```julia
abstract type AbstractPixel
```


Abstract Pixel Type


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/cameras/camera_types.jl#L15)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_w_I0_terms_integrals_case2' href='#Krang.radial_w_I0_terms_integrals_case2'>#</a>&nbsp;<b><u>Krang.radial_w_I0_terms_integrals_case2</u></b> &mdash; <i>Function</i>.




```julia
radial_w_I0_terms_integrals_case2(
    metric::Kerr{T},
    rs,
    roots::NTuple{4, T} where T,
    τ,
    νr
) -> Tuple

```


Returns the radial integrals for the case where there are four real roots in the radial potential, with roots outside the horizon.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1457)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.SlowLightIntensityScreen' href='#Krang.SlowLightIntensityScreen'>#</a>&nbsp;<b><u>Krang.SlowLightIntensityScreen</u></b> &mdash; <i>Type</i>.




```julia
struct SlowLightIntensityScreen{T} <: Krang.AbstractScreen
```


Screen made of Intensity Pixels.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/cameras/SlowLightIntensityCamera.jl#L59)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_w_I0_terms_integrals_case4' href='#Krang.radial_w_I0_terms_integrals_case4'>#</a>&nbsp;<b><u>Krang.radial_w_I0_terms_integrals_case4</u></b> &mdash; <i>Function</i>.




```julia
radial_w_I0_terms_integrals_case4(
    metric::Kerr{T},
    rs,
    roots::NTuple{4, T} where T,
    τ
) -> Tuple

```


Returns the radial integrals for the case where there are no real roots in the radial potential


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1556)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_integrals' href='#Krang.radial_integrals'>#</a>&nbsp;<b><u>Krang.radial_integrals</u></b> &mdash; <i>Function</i>.




```julia
radial_integrals(
    pix::Krang.AbstractPixel,
    rs,
    τ,
    νr
) -> NTuple{5, Any}

```


Return the radial integrals
- `pix`: SlowLightIntensityPixel
  
- `rs` : Emission radius
  
- `τ` : Mino time
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1808)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.AbstractCamera' href='#Krang.AbstractCamera'>#</a>&nbsp;<b><u>Krang.AbstractCamera</u></b> &mdash; <i>Type</i>.




```julia
abstract type AbstractCamera
```


Abstract Observer Type


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/cameras/camera_types.jl#L1)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.IntensityScreen' href='#Krang.IntensityScreen'>#</a>&nbsp;<b><u>Krang.IntensityScreen</u></b> &mdash; <i>Type</i>.




```julia
struct IntensityScreen{T} <: Krang.AbstractScreen
```


Screen made of Intensity Pixels.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/cameras/IntensityCamera.jl#L40)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.It' href='#Krang.It'>#</a>&nbsp;<b><u>Krang.It</u></b> &mdash; <i>Function</i>.




```julia
It(pix::Krang.AbstractPixel, rs, τ, νr) -> Any

```


Returns the antiderivative $I_t=\int\frac{a(2Mr-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `pix`: SlowLightIntensityPixel
  
- `rs` : Emission radius
  
- `τ` : Mino time
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1788)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.IntensityPixel' href='#Krang.IntensityPixel'>#</a>&nbsp;<b><u>Krang.IntensityPixel</u></b> &mdash; <i>Type</i>.




```julia
struct IntensityPixel{T} <: Krang.AbstractPixel
```


Intensity Pixel Type.  Each Pixel is associated with a single ray, and caches some information about the ray.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/cameras/IntensityCamera.jl#L2)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Iϕ_m_I0_terms' href='#Krang.Iϕ_m_I0_terms'>#</a>&nbsp;<b><u>Krang.Iϕ_m_I0_terms</u></b> &mdash; <i>Function</i>.




```julia
Iϕ_m_I0_terms(metric::Kerr{T}, rs, roots, νr, λ) -> Any

```


Returns the antiderivative $I_ϕ=\int\frac{a(2Mr-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$ with I0 terms. 

See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `τ`: Minotime 
  
- `roots` : Radial roots
  
- `λ`  : Reduced azimuthal angular momentum
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L570)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.AbstractMaterial' href='#Krang.AbstractMaterial'>#</a>&nbsp;<b><u>Krang.AbstractMaterial</u></b> &mdash; <i>Type</i>.




```julia
abstract type AbstractMaterial
```


Abstract Material


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/materials/material_types.jl#L1)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.AbstractScreen' href='#Krang.AbstractScreen'>#</a>&nbsp;<b><u>Krang.AbstractScreen</u></b> &mdash; <i>Type</i>.




```julia
abstract type AbstractScreen
```


Abstract Screen Type


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/cameras/camera_types.jl#L8)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_w_I0_terms_integrals_case3' href='#Krang.radial_w_I0_terms_integrals_case3'>#</a>&nbsp;<b><u>Krang.radial_w_I0_terms_integrals_case3</u></b> &mdash; <i>Function</i>.




```julia
radial_w_I0_terms_integrals_case3(
    metric::Kerr{T},
    rs,
    roots::NTuple{4, T} where T,
    τ
) -> Tuple

```


Returns the radial integrals for the case where there are two real roots in the radial potential


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1509)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.IntensityCamera' href='#Krang.IntensityCamera'>#</a>&nbsp;<b><u>Krang.IntensityCamera</u></b> &mdash; <i>Type</i>.




```julia
struct IntensityCamera{T} <: Krang.AbstractCamera
```


Observer sitting at radial infinity. The frame of this observer is alligned with the Boyer-Lindquist frame.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/cameras/IntensityCamera.jl#L68)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.Ir_s' href='#Krang.Ir_s'>#</a>&nbsp;<b><u>Krang.Ir_s</u></b> &mdash; <i>Function</i>.




```julia
Ir_s(metric::Kerr{T}, rs, roots, νr) -> Any

```


Returns the antiderivative $I_r=\int\frac{dr}{\sqrt{\mathcal{R(r)}}}$. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `rs` : Emission radius
  
- `roots`  : Roots of the radial potential
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L368)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.AbstractGeometry' href='#Krang.AbstractGeometry'>#</a>&nbsp;<b><u>Krang.AbstractGeometry</u></b> &mdash; <i>Type</i>.




```julia
abstract type AbstractGeometry
```


Abstract Geometry


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/geometries/geometry_types.jl#L2)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.radial_inf_integrals_case2' href='#Krang.radial_inf_integrals_case2'>#</a>&nbsp;<b><u>Krang.radial_inf_integrals_case2</u></b> &mdash; <i>Function</i>.




```julia
radial_inf_integrals_case2(
    metric::Kerr{T},
    roots::NTuple{4, T} where T
) -> NTuple{4, Any}

```


Returns the radial integrals for the case where there are four real roots in the radial potential, with roots outside the horizon.


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L1321)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='Krang.It_inf' href='#Krang.It_inf'>#</a>&nbsp;<b><u>Krang.It_inf</u></b> &mdash; <i>Function</i>.




```julia
It_inf(
    metric::Kerr{T},
    roots::NTuple{4, T} where T,
    λ
) -> Any

```


Returns the antiderivative $I_t=\int\frac{r^2\Delta+2Mr(r^2+a^2-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$  evaluated at infinity without I0 terms.

See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `roots` : Radial roots
  
- `λ`  : Reduced azimuthal angular momentum
  


[source](https://github.com/dchang10/Krang.jl/blob/1915c3885ae3114730d1f3b649882f5c7b2f6934/src/metrics/Kerr/misc.jl#L828)

</div>
<br>
