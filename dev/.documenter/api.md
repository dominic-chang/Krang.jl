
# Krang api {#Krang-api}
- [`Krang.AbstractCamera`](#Krang.AbstractCamera)
- [`Krang.AbstractGeometry`](#Krang.AbstractGeometry)
- [`Krang.AbstractLevelSetGeometry`](#Krang.AbstractLevelSetGeometry)
- [`Krang.AbstractMaterial`](#Krang.AbstractMaterial)
- [`Krang.AbstractMetric`](#Krang.AbstractMetric)
- [`Krang.AbstractPixel`](#Krang.AbstractPixel)
- [`Krang.AbstractScreen`](#Krang.AbstractScreen)
- [`Krang.ConeGeometry`](#Krang.ConeGeometry)
- [`Krang.ElectronSynchrotronPowerLawIntensity`](#Krang.ElectronSynchrotronPowerLawIntensity)
- [`Krang.ElectronSynchrotronPowerLawPolarization`](#Krang.ElectronSynchrotronPowerLawPolarization)
- [`Krang.IntensityCamera`](#Krang.IntensityCamera)
- [`Krang.IntensityPixel`](#Krang.IntensityPixel)
- [`Krang.IntensityScreen`](#Krang.IntensityScreen)
- [`Krang.Kerr`](#Krang.Kerr)
- [`Krang.Mesh`](#Krang.Mesh)
- [`Krang.MeshGeometry`](#Krang.MeshGeometry)
- [`Krang.SlowLightIntensityCamera`](#Krang.SlowLightIntensityCamera)
- [`Krang.SlowLightIntensityPixel`](#Krang.SlowLightIntensityPixel)
- [`Krang.SlowLightIntensityScreen`](#Krang.SlowLightIntensityScreen)
- [`Krang.Gθ`](#Krang.Gθ)
- [`Krang.I0_inf`](#Krang.I0_inf)
- [`Krang.I1_inf_m_I0_terms`](#Krang.I1_inf_m_I0_terms)
- [`Krang.I2_inf_m_I0_terms`](#Krang.I2_inf_m_I0_terms)
- [`Krang.Im_inf_m_I0_terms`](#Krang.Im_inf_m_I0_terms)
- [`Krang.Ip_inf_m_I0_terms`](#Krang.Ip_inf_m_I0_terms)
- [`Krang.Ir`](#Krang.Ir)
- [`Krang.Ir_inf`](#Krang.Ir_inf)
- [`Krang.Ir_s`](#Krang.Ir_s)
- [`Krang.It`](#Krang.It)
- [`Krang.It_inf`](#Krang.It_inf)
- [`Krang.It_w_I0_terms`](#Krang.It_w_I0_terms)
- [`Krang.Iϕ`](#Krang.Iϕ)
- [`Krang.Iϕ_inf`](#Krang.Iϕ_inf)
- [`Krang.Iϕ_w_I0_terms`](#Krang.Iϕ_w_I0_terms)
- [`Krang._isreal2`](#Krang._isreal2)
- [`Krang.absGto_Gthat`](#Krang.absGto_Gthat)
- [`Krang.absGθo_Gθhat`](#Krang.absGθo_Gθhat)
- [`Krang.absGϕo_Gϕhat`](#Krang.absGϕo_Gϕhat)
- [`Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild`](#Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild)
- [`Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light`](#Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light)
- [`Krang.emission_azimuth`](#Krang.emission_azimuth)
- [`Krang.emission_coordinates`](#Krang.emission_coordinates)
- [`Krang.emission_coordinates_fast_light`](#Krang.emission_coordinates_fast_light)
- [`Krang.emission_inclination`](#Krang.emission_inclination)
- [`Krang.emission_radius`](#Krang.emission_radius)
- [`Krang.get_radial_roots`](#Krang.get_radial_roots)
- [`Krang.horizon`](#Krang.horizon)
- [`Krang.inclination`](#Krang.inclination)
- [`Krang.jac_bl_d_zamo_u`](#Krang.jac_bl_d_zamo_u)
- [`Krang.jac_bl_u_zamo_d`](#Krang.jac_bl_u_zamo_d)
- [`Krang.jac_fluid_u_zamo_d`](#Krang.jac_fluid_u_zamo_d)
- [`Krang.jac_zamo_d_bl_u`](#Krang.jac_zamo_d_bl_u)
- [`Krang.jac_zamo_u_bl_d`](#Krang.jac_zamo_u_bl_d)
- [`Krang.metric`](#Krang.metric)
- [`Krang.metric_dd`](#Krang.metric_dd)
- [`Krang.metric_uu`](#Krang.metric_uu)
- [`Krang.mino_time`](#Krang.mino_time)
- [`Krang.p_bl_d`](#Krang.p_bl_d)
- [`Krang.penrose_walker`](#Krang.penrose_walker)
- [`Krang.r_potential`](#Krang.r_potential)
- [`Krang.radial_inf_integrals_case2`](#Krang.radial_inf_integrals_case2)
- [`Krang.radial_inf_integrals_case3`](#Krang.radial_inf_integrals_case3)
- [`Krang.radial_inf_integrals_case4`](#Krang.radial_inf_integrals_case4)
- [`Krang.radial_inf_integrals_m_I0_terms`](#Krang.radial_inf_integrals_m_I0_terms)
- [`Krang.radial_integrals`](#Krang.radial_integrals)
- [`Krang.radial_w_I0_terms_integrals_case2`](#Krang.radial_w_I0_terms_integrals_case2)
- [`Krang.radial_w_I0_terms_integrals_case3`](#Krang.radial_w_I0_terms_integrals_case3)
- [`Krang.radial_w_I0_terms_integrals_case4`](#Krang.radial_w_I0_terms_integrals_case4)
- [`Krang.regularized_Pi`](#Krang.regularized_Pi)
- [`Krang.roots`](#Krang.roots)
- [`Krang.screen_coordinate`](#Krang.screen_coordinate)
- [`Krang.screen_polarization`](#Krang.screen_polarization)
- [`Krang.synchrotronIntensity`](#Krang.synchrotronIntensity)
- [`Krang.synchrotronPolarization`](#Krang.synchrotronPolarization)
- [`Krang.total_mino_time`](#Krang.total_mino_time)
- [`Krang.α`](#Krang.α)
- [`Krang.αboundary`](#Krang.αboundary)
- [`Krang.β`](#Krang.β)
- [`Krang.βboundary`](#Krang.βboundary)
- [`Krang.η`](#Krang.η)
- [`Krang.θ_potential`](#Krang.θ_potential)
- [`Krang.λ`](#Krang.λ)


### Raytracing Functions {#Raytracing-Functions}
<details class='jldocstring custom-block' open>
<summary><a id='Krang.emission_radius' href='#Krang.emission_radius'><span class="jlbinding">Krang.emission_radius</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
emission_radius(
    pix::Krang.AbstractPixel,
    θs,
    isindir,
    n
) -> Tuple{Any, Any, Any, Any, Bool}

```


Emission radius for point originating at inclination θs whose nth order image appears at the screen coordinate (`α`, `β`).  Returns 0 if the emission coordinates do not exist for that screen coordinate.

**Arguments**
- `pix` : Pixel information 
  
- `θs` : Emission inclination
  
- `isindir` : Is emission to observer direct or indirect
  
- `n` : Image index
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/emission_coordinates.jl#L21" target="_blank" rel="noreferrer">source</a></Badge>



```julia
emission_radius(
    pix::Krang.AbstractPixel,
    τ
) -> Tuple{Any, Any, Any, Bool}

```


Emission radius for point originating at at Mino time τ whose image appears at the screen coordinate (`α`, `β`).  Returns 0 if the emission coordinates do not exist for that screen coordinate.

**Arguments**

-`pix` : Pixel information
- `τ` : Mino time
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/emission_coordinates.jl#L61" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.emission_inclination' href='#Krang.emission_inclination'><span class="jlbinding">Krang.emission_inclination</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
emission_inclination(
    pix::Krang.AbstractPixel,
    rs,
    νr
) -> NTuple{6, Any}

```


Emission inclination for point originating at inclination rs whose nth order image appears at screen coordinate (`α`, `β`).

**Arguments**
- `pix` : Pixel information
  
- `rs` : Emission radius
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/emission_coordinates.jl#L90" target="_blank" rel="noreferrer">source</a></Badge>



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
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/emission_coordinates.jl#L103" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.emission_coordinates_fast_light' href='#Krang.emission_coordinates_fast_light'><span class="jlbinding">Krang.emission_coordinates_fast_light</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
emission_coordinates_fast_light(
    pix::Krang.AbstractPixel,
    θs,
    isindir,
    n
) -> Tuple{Any, Any, Any, Any, Bool}

```


Emission radius and azimuthal angle for point originating at inclination θs whose nth order image appears at the screen  coordinate (`α`, `β`) for an observer located at inclination θo.

**Arguments**
- `pix` : Pixel information
  
- `θs` : Emission Inclination
  
- `isindir` : Whether emission to observer is direct or indirect
  
- `n` : Image index
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/emission_coordinates.jl#L145" target="_blank" rel="noreferrer">source</a></Badge>



```julia
emission_coordinates_fast_light(
    pix::Krang.AbstractPixel,
    τ
) -> Tuple{Any, Any, Any, Any, Any, Bool}

```


Ray trace a point that appears at the screen coordinate (`α`, `β`) for an observer located at inclination θo

**Arguments**
- `pix` : Pixel information
  
- `τ` : Mino Time
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/emission_coordinates.jl#L185" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.emission_coordinates' href='#Krang.emission_coordinates'><span class="jlbinding">Krang.emission_coordinates</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
emission_coordinates(
    pix::Krang.AbstractPixel,
    θs,
    isindir,
    n
) -> Tuple{Any, Any, Any, Any, Any, Bool}

```


Emission radius and azimuthal angle for point originating at inclination θs whose nth order image appears at the screen coordinate (`α`, `β`) for an observer located at inclination θo.

**Arguments**
- `pix` : Pixel information
  
- `θs` : Emission Inclination
  
- `isindir` : Whether emission to observer is direct or indirect
  
- `n` : Image index
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/emission_coordinates.jl#L231" target="_blank" rel="noreferrer">source</a></Badge>



```julia
emission_coordinates(
    pix::Krang.AbstractPixel,
    τ
) -> Tuple{Any, Any, Any, Any, Any, Any, Bool}

```


Ray trace a point that appears at the screen coordinate (`α`, `β`) for an observer located at inclination θo

**Arguments**
- `pix` : Pixel information
  
- `τ` : Mino Time
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/emission_coordinates.jl#L300" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.mino_time' href='#Krang.mino_time'><span class="jlbinding">Krang.mino_time</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
mino_time(pix, rs, isindir) -> Any

```


Mino time of trajectory between an observer at infinity and point at radius rs

**Arguments**
- `pix` : Pixel information
  
- `rs` : Emission radius
  
- `isindir` : Is the path direct or indirect?
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L300" target="_blank" rel="noreferrer">source</a></Badge>



```julia
mino_time(pix::Krang.AbstractPixel, θs, isindir, n) -> Any

```


Mino time of trajectory between two inclinations for a given screen coordinate

**Arguments**
- `pix` : Pixel information
  
- `θs` : Emission inclination
  
- `isindir` : Is the path direct or indirect?
  
- `n` : nth image in orde of amount of minotime traversed
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1680" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Metric Functions {#Metric-Functions}
<details class='jldocstring custom-block' open>
<summary><a id='Krang.AbstractMetric' href='#Krang.AbstractMetric'><span class="jlbinding">Krang.AbstractMetric</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type AbstractMetric
```


Abstract Metric Type


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/AbstractMetric.jl#L1" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Kerr' href='#Krang.Kerr'><span class="jlbinding">Krang.Kerr</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct Kerr{T} <: Krang.AbstractMetric
```


Kerr Metric in Boyer Lindquist Coordinates
- `mass`: M = mass
  
- `spin`: a = J/M, where J is the angular momentum and M is the mass of the black hole.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/Kerr.jl#L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.metric_uu' href='#Krang.metric_uu'><span class="jlbinding">Krang.metric_uu</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
metric_uu(metric::Krang.AbstractMetric, args...) -> Any

```


Returns the inverse metric in some representation (usually as an nxn matrix).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/AbstractMetric.jl#L15" target="_blank" rel="noreferrer">source</a></Badge>



```julia
metric_uu(metric::Kerr{T}, r, θ) -> Any

```


Inverse Kerr metric in Boyer Lindquist (BL) coordinates.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/Kerr.jl#L59" target="_blank" rel="noreferrer">source</a></Badge>



```julia
metric_uu(metric::Kerr{T}, coordinates) -> Any

```


Inverse Kerr metric in Boyer Lindquist (BL) coordinates.

**Arguments**
- `metric` : Kerr metric
  
- `coordinates` : Coordinates (t, r, θ, ϕ)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/Kerr.jl#L77" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.metric_dd' href='#Krang.metric_dd'><span class="jlbinding">Krang.metric_dd</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
metric_dd(metric::Krang.AbstractMetric, args...) -> Any

```


Returns the metric in some representation (usually as an nxn matrix).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/AbstractMetric.jl#L8" target="_blank" rel="noreferrer">source</a></Badge>



```julia
metric_dd(metric::Kerr{T}, r, θ) -> Any

```


Kerr metric in Boyer Lindquist (BL) coordinates.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/Kerr.jl#L90" target="_blank" rel="noreferrer">source</a></Badge>



```julia
metric_dd(metric::Kerr{T}, coordinates) -> Any

```


Kerr metric in Boyer Lindquist (BL) coordinates.

**Arguments**
- `metric` : Kerr metric
  
- `coordinates` : Coordinates (t, r, θ, ϕ)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/Kerr.jl#L110" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.horizon' href='#Krang.horizon'><span class="jlbinding">Krang.horizon</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
horizon(metric::Kerr{T}) -> Any

```


Outer Horizon for the Kerr metric.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/Kerr.jl#L32" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.λ' href='#Krang.λ'><span class="jlbinding">Krang.λ</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
λ(pix::Krang.AbstractPixel) -> Any

```


```
λ(pix::AbstractPixel)
```


Calculate the λ value for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The λ value of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L85" target="_blank" rel="noreferrer">source</a></Badge>



```julia
λ(_::Kerr, α, θo) -> Any

```


Energy reduced azimuthal angular momentum

**Arguments**
- `metric`: Kerr
  
- `α`: Horizontal Bardeen screen coordinate
  
- `θo`: Observer inclination
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L125" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.η' href='#Krang.η'><span class="jlbinding">Krang.η</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
η(pix::Krang.AbstractPixel) -> Any

```


```
η(pix::AbstractPixel)
```


Calculate the η value for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The η value of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L67" target="_blank" rel="noreferrer">source</a></Badge>



```julia
η(metric::Kerr, α, β, θo) -> Any

```


Energy reduced Carter integral

**Arguments**
- `metric`: Kerr
  
- `α`: Horizontal Bardeen screen coordinate
  
- `β`: Bardeen vertical coordinate
  
- `θo`: Observer inclination
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L138" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Radial Integrals {#Radial-Integrals}
<details class='jldocstring custom-block' open>
<summary><a id='Krang.r_potential' href='#Krang.r_potential'><span class="jlbinding">Krang.r_potential</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
r_potential(metric::Kerr{T}, η, λ, r) -> Any

```


Radial potential of spacetime

**Arguments**
- `metric`: Kerr{T} metric
  
- `η`  : Reduced Carter constant
  
- `λ`  : Reduced azimuthal angular momentum
  
- `r`  : Boyer Lindquist radius
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L210" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.get_radial_roots' href='#Krang.get_radial_roots'><span class="jlbinding">Krang.get_radial_roots</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_radial_roots(metric::Kerr{T}, η, λ) -> NTuple{4, Any}

```


Returns roots of $r^4 + (a^2-η-λ^2)r^2 + 2(η+(a-λ)^2)r - a^2η$

**Arguments**
- `metric`: Kerr{T} metric
  
- `η`  : Reduced Carter constant
  
- `λ`  : Reduced azimuthal angular momentum
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L245" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Ir' href='#Krang.Ir'><span class="jlbinding">Krang.Ir</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Ir(pix::Krang.AbstractPixel, νr::Bool, rs) -> Any

```


Returns the antiderivative $I_r=\int\frac{dr}{\sqrt{\mathcal{R(r)}}}$. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `pix`  : Pixel information
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  
- `rs` : Emission radius
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1526" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Angular Integrals {#Angular-Integrals}
<details class='jldocstring custom-block' open>
<summary><a id='Krang.θ_potential' href='#Krang.θ_potential'><span class="jlbinding">Krang.θ_potential</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
θ_potential(metric::Kerr{T}, η, λ, θ) -> Any

```


Theta potential of a Kerr black hole

**Arguments**
- `metric`: Kerr{T} metric
  
- `η`  : Reduced Carter constant
  
- `λ`  : Reduced azimuthal angular momentum
  
- `θ`  : Boyer Lindquist inclination
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L227" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Gθ' href='#Krang.Gθ'><span class="jlbinding">Krang.Gθ</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Gθ(
    pix::Krang.AbstractPixel,
    θs,
    isindir,
    n
) -> Tuple{Any, Any, Any, Any, Any, Bool}

```


Returns the antiderivative $G_\theta=\int\frac{d\theta}{\sqrt{\Theta(\theta)}}$. See [`θ_potential(x)`](/api#Krang.θ_potential) for an implementation of $\Theta(	heta)$.

**Arguments**
- `pix` : Pixel information
  
- `θs` : Emission inclination
  
- `isindir` : Is the path direct or indirect?
  
- `n` : nth image ordered by minotime
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1694" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Screen Coordinates {#Screen-Coordinates}
<details class='jldocstring custom-block' open>
<summary><a id='Krang.α' href='#Krang.α'><span class="jlbinding">Krang.α</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
α(_::Kerr, λ, θo) -> Any

```


Horizontal Bardeen Screen Coordinate

**Arguments**
- `metric`: Kerr
  
- `α`: Horizontal Bardeen screen coordinate
  
- `θo`: Observer inclination
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L152" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.β' href='#Krang.β'><span class="jlbinding">Krang.β</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
β(metric::Kerr, λ, η, θo) -> Any

```


Vertical Bardeen Screen Coordinate

**Arguments**
- `metric`: Kerr
  
- `λ`: Energy reduced Azimuthal angular momentul
  
- `η`: Energy reduced Carter integral 
  
- `θo`: Observer inclination
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L165" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.αboundary' href='#Krang.αboundary'><span class="jlbinding">Krang.αboundary</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
αboundary(metric::Kerr, θs) -> Any

```


Defines a horizontal boundary on the assymptotic observer&#39;s screen that emission that from θs must fall within.

**Arguments**
- `metric`: Kerr metric
  
- `θs`  : Emission Inclination
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L179" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.βboundary' href='#Krang.βboundary'><span class="jlbinding">Krang.βboundary</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
βboundary(metric::Kerr{T}, α, θo, θs) -> Any

```


Defines a vertical boundary on the assymptotic observer&#39;s screen that emission that from θs must fall within.

**Arguments**
- `metric`: Kerr{T} metric
  
- `α`   : Horizontal Bardeen screen coordinate
  
- `θo`  : Observer inclination
  
- `θs`  : Emission Inclination
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L191" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Local Frame transformations {#Local-Frame-transformations}
<details class='jldocstring custom-block' open>
<summary><a id='Krang.p_bl_d' href='#Krang.p_bl_d'><span class="jlbinding">Krang.p_bl_d</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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



<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/physicsUtils.jl#L14" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.jac_bl_u_zamo_d' href='#Krang.jac_bl_u_zamo_d'><span class="jlbinding">Krang.jac_bl_u_zamo_d</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
jac_bl_u_zamo_d(metric::Kerr{T}, r, θ) -> Any

```


Jacobian which converts ZAMO vector to a Boyer-Lindquist basis


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/physicsUtils.jl#L81" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.jac_zamo_u_bl_d' href='#Krang.jac_zamo_u_bl_d'><span class="jlbinding">Krang.jac_zamo_u_bl_d</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
jac_zamo_u_bl_d(metric::Kerr{T}, r, θ) -> Any

```


Jacobian which converts Boyer-Lindquist vector to a ZAMO basis


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/physicsUtils.jl#L63" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.jac_bl_d_zamo_u' href='#Krang.jac_bl_d_zamo_u'><span class="jlbinding">Krang.jac_bl_d_zamo_u</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
jac_bl_d_zamo_u(metric::Kerr{T}, r, θ) -> Any

```


Jacobian which converts ZAMO covector to a Boyer-Lindquist basis


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/physicsUtils.jl#L45" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.jac_zamo_d_bl_u' href='#Krang.jac_zamo_d_bl_u'><span class="jlbinding">Krang.jac_zamo_d_bl_u</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
jac_zamo_d_bl_u(metric::Kerr{T}, r, θ) -> Any

```


Returns the Jacobian which converts a Boyer-Lindquist covector to ZAMO basis.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/physicsUtils.jl#L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.jac_fluid_u_zamo_d' href='#Krang.jac_fluid_u_zamo_d'><span class="jlbinding">Krang.jac_fluid_u_zamo_d</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
jac_fluid_u_zamo_d(_::Kerr{T}, β, θ, φ) -> Any

```


Jacobian which expreases ZAMO vector in the fluid frame


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/physicsUtils.jl#L99" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Polarization {#Polarization}
<details class='jldocstring custom-block' open>
<summary><a id='Krang.screen_polarization' href='#Krang.screen_polarization'><span class="jlbinding">Krang.screen_polarization</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
screen_polarization(
    metric::Kerr{T},
    κ::Complex,
    θ,
    α,
    β
) -> Tuple{Any, Any}

```


Returns the screen polarization associated with a killing spinor κ as seen seen by an asymptotic observer.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/ElectronSynchrotronPowerLawPolarization.jl#L1" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.penrose_walker' href='#Krang.penrose_walker'><span class="jlbinding">Krang.penrose_walker</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/physicsUtils.jl#L117" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.synchrotronIntensity' href='#Krang.synchrotronIntensity'><span class="jlbinding">Krang.synchrotronIntensity</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
synchrotronIntensity(
    metric::Kerr{T},
    α,
    β,
    ri,
    θs,
    θo,
    magnetic_field::StaticArraysCore.SArray{Tuple{3}, T, 1, 3},
    βfluid::StaticArraysCore.SArray{Tuple{3}, T, 1, 3},
    νr::Bool,
    θsign::Bool
) -> Tuple{Any, Any, Any}

```


Calculates the intensity of a photon emitted from a fluid particle with momentum f_u and observed by an asymptotic observer.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/ElectronSynchrotronPowerLawIntensity.jl#L1" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.synchrotronPolarization' href='#Krang.synchrotronPolarization'><span class="jlbinding">Krang.synchrotronPolarization</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
synchrotronPolarization(
    metric::Kerr{T},
    α,
    β,
    ri,
    θs,
    θo,
    magnetic_field::StaticArraysCore.SArray{Tuple{3}, T, 1, 3},
    βfluid::StaticArraysCore.SArray{Tuple{3}, T, 1, 3},
    νr::Bool,
    θsign::Bool
) -> NTuple{4, Any}

```


Calculates the polarization of a photon emitted from a fluid particle with momentum f_u and observed by an asymptotic observer.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/ElectronSynchrotronPowerLawPolarization.jl#L19" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Meshes, Geometries and Materials {#Meshes,-Geometries-and-Materials}
<details class='jldocstring custom-block' open>
<summary><a id='Krang.Mesh' href='#Krang.Mesh'><span class="jlbinding">Krang.Mesh</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct Mesh{G<:Krang.AbstractGeometry, M<:Krang.AbstractMaterial}
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/geometries/geometry_types.jl#L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.MeshGeometry' href='#Krang.MeshGeometry'><span class="jlbinding">Krang.MeshGeometry</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Load a mesh from a file.

# Arguments
- `filename::String`: The path to the file containing the mesh.

# Returns
- A `Mesh` object representing the mesh.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/geometries/mesh_geometry.jl#L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.AbstractLevelSetGeometry' href='#Krang.AbstractLevelSetGeometry'><span class="jlbinding">Krang.AbstractLevelSetGeometry</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Level Sets should be a functor
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/geometries/level_set_geometry.jl#L1" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.ElectronSynchrotronPowerLawIntensity' href='#Krang.ElectronSynchrotronPowerLawIntensity'><span class="jlbinding">Krang.ElectronSynchrotronPowerLawIntensity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct ElectronSynchrotronPowerLawIntensity{N, T} <: Krang.AbstractMaterial
```


A struct representing the linear polarization intensity of synchrotron radiation from electrons following a power-law energy distribution. This model is based on the material described in https://doi.org/10.3847/1538-4357/abf117.

**Fields**

```
- `magnetic_field::SVector{3, T}`: The magnetic field vector components (x, y, z).
- `fluid_velocity::SVector{3, T}`: The fluid velocity vector components (speed, inclination angle, azimuthal angle).
- `spectral_index::T`: The spectral index of the electron energy distribution.
- `R::T`: The characteristic radius of the emissivity profile.
- `p1::T`: The first power-law index.
- `p2::T`: The second power-law index.
- `subimgs::NTuple{N, Int}`: The sub-images to ray trace.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/ElectronSynchrotronPowerLawIntensity.jl#L42" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.ElectronSynchrotronPowerLawPolarization' href='#Krang.ElectronSynchrotronPowerLawPolarization'><span class="jlbinding">Krang.ElectronSynchrotronPowerLawPolarization</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct ElectronSynchrotronPowerLawPolarization{N, T} <: Krang.AbstractMaterial
```


A struct representing the linear polarization of synchrotron radiation from electrons following a power-law energy distribution. This model is based on the material described in https://doi.org/10.3847/1538-4357/abf117.

**Fields**

```
- `magnetic_field::SVector{3, T}`: The magnetic field vector components (x, y, z).
- `fluid_velocity::SVector{3, T}`: The fluid velocity vector components (speed, inclination angle, azimuthal angle).
- `spectral_index::T`: The spectral index of the electron energy distribution.
- `R::T`: The characteristic radius of the emissivity profile.
- `p1::T`: The first power-law index.
- `p2::T`: The second power-law index.
- `subimgs::NTuple{N, Int}`: The sub-images to ray trace.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/ElectronSynchrotronPowerLawPolarization.jl#L83" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Pixel Related Functions {#Pixel-Related-Functions}
<details class='jldocstring custom-block' open>
<summary><a id='Krang.absGθo_Gθhat' href='#Krang.absGθo_Gθhat'><span class="jlbinding">Krang.absGθo_Gθhat</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
absGθo_Gθhat(
    pix::Krang.AbstractPixel
) -> Tuple{T, T} where T

```


```
absGθo_Gθhat(pix::AbstractPixel)
```


Calculate the absolute value of Gθo divided by Gθhat for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The absolute value of Gθo divided by Gθhat of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L266" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.I1_inf_m_I0_terms' href='#Krang.I1_inf_m_I0_terms'><span class="jlbinding">Krang.I1_inf_m_I0_terms</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
I1_inf_m_I0_terms(pix::Krang.AbstractPixel) -> Any

```


```
I1_inf_m_I0_terms(pix::AbstractPixel)
```


Calculate the I1 infinity minus I0 terms for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The I1 infinity minus I0 terms of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L162" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.radial_inf_integrals_m_I0_terms' href='#Krang.radial_inf_integrals_m_I0_terms'><span class="jlbinding">Krang.radial_inf_integrals_m_I0_terms</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
radial_inf_integrals_m_I0_terms(
    pix::Krang.AbstractPixel
) -> NTuple{4, Any}

```


```
radial_inf_integrals_m_I0_terms(pix::AbstractPixel)
```


Calculate the radial infinity minus I0 terms for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The radial infinity minus I0 terms of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L221" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.inclination' href='#Krang.inclination'><span class="jlbinding">Krang.inclination</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
inclination(pix::Krang.AbstractPixel) -> Any

```


```
inclination(pix::AbstractPixel)
```


Get the inclination angle (θo) of a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The inclination angle (θo) of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L52" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.metric' href='#Krang.metric'><span class="jlbinding">Krang.metric</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
metric(pix::Krang.AbstractPixel) -> Any

```


```
metric(pix::AbstractPixel)
```


Get the space time metric.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The space time metric
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.absGto_Gthat' href='#Krang.absGto_Gthat'><span class="jlbinding">Krang.absGto_Gthat</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
absGto_Gthat(
    pix::Krang.AbstractPixel
) -> Tuple{T, T} where T

```


```
absGto_Gthat(pix::AbstractPixel)
```


Calculate the absolute value of Gto divided by Gthat for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The absolute value of Gto divided by Gthat of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L296" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.screen_coordinate' href='#Krang.screen_coordinate'><span class="jlbinding">Krang.screen_coordinate</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
screen_coordinate(
    pix::Krang.AbstractPixel
) -> Tuple{T, T} where T

```


```
screen_coordinate(pix::AbstractPixel)
```


Get the screen coordinates of a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel for which to get the screen coordinates.
  

**Returns**
- The screen coordinates of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L22" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Ip_inf_m_I0_terms' href='#Krang.Ip_inf_m_I0_terms'><span class="jlbinding">Krang.Ip_inf_m_I0_terms</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Ip_inf_m_I0_terms(pix::Krang.AbstractPixel) -> Any

```


```
Ip_inf_m_I0_terms(pix::AbstractPixel)
```


Calculate the Ip infinity minus I0 terms for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The Ip infinity minus I0 terms of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L191" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.I2_inf_m_I0_terms' href='#Krang.I2_inf_m_I0_terms'><span class="jlbinding">Krang.I2_inf_m_I0_terms</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
I2_inf_m_I0_terms(pix::Krang.AbstractPixel) -> Any

```


```
I2_inf_m_I0_terms(pix::AbstractPixel)
```


Calculate the I2 infinity minus I0 terms for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The I2 infinity minus I0 terms of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L177" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.roots' href='#Krang.roots'><span class="jlbinding">Krang.roots</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
roots(
    pix::Krang.AbstractPixel
) -> NTuple{4, Complex{T}} where T

```


```
roots(pix::AbstractPixel)
```


Calculate the radial roots for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The radial roots of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L103" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Im_inf_m_I0_terms' href='#Krang.Im_inf_m_I0_terms'><span class="jlbinding">Krang.Im_inf_m_I0_terms</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Im_inf_m_I0_terms(pix::Krang.AbstractPixel) -> Any

```


```
Im_inf_m_I0_terms(pix::AbstractPixel)
```


Calculate the Im infinity minus I0 terms for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The Im infinity minus I0 terms of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L206" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.absGϕo_Gϕhat' href='#Krang.absGϕo_Gϕhat'><span class="jlbinding">Krang.absGϕo_Gϕhat</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
absGϕo_Gϕhat(
    pix::Krang.AbstractPixel
) -> Tuple{T, T} where T

```


```
absGϕo_Gϕhat(pix::AbstractPixel)
```


Calculate the absolute value of Gϕo divided by Gϕhat for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The absolute value of Gϕo divided by Gϕhat of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L281" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.I0_inf' href='#Krang.I0_inf'><span class="jlbinding">Krang.I0_inf</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
I0_inf(pix::Krang.AbstractPixel) -> Any

```


```
I0_inf(pix::AbstractPixel)
```


Calculate the I0 infinity value for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The I0 infinity value of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L118" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.total_mino_time' href='#Krang.total_mino_time'><span class="jlbinding">Krang.total_mino_time</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
total_mino_time(pix::Krang.AbstractPixel) -> Any

```


```
total_mino_time(pix::AbstractPixel)
```


Return the total possible Mino time for a ray associated with a pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The total possible Mino time for a ray associated with the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L133" target="_blank" rel="noreferrer">source</a></Badge>



```julia
total_mino_time(
    metric::Kerr{T},
    roots::NTuple{4, T} where T
) -> Any

```


Returns the maximum mino time that can be accrued along a ray.

**Arguments**
- `metric` : Kerr metric
  
- `roots` : Roots of the radial potential
  
- `I0_inf` : Mino time at infinity
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1500" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Misc {#Misc}
<details class='jldocstring custom-block' open>
<summary><a id='Krang._isreal2' href='#Krang._isreal2'><span class="jlbinding">Krang._isreal2</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
_isreal2(num::Complex{T}) -> Any

```


Checks if a complex number is real to some tolerance


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.regularized_Pi' href='#Krang.regularized_Pi'><span class="jlbinding">Krang.regularized_Pi</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
regularized_Pi(n, ϕ, k) -> Any

```


Regularized elliptic integral of the third kind

**Arguments**
- `n`: Parameter
  
- `ϕ`: Arguments
  
- `k`: Parameter
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.radial_inf_integrals_case4' href='#Krang.radial_inf_integrals_case4'><span class="jlbinding">Krang.radial_inf_integrals_case4</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
radial_inf_integrals_case4(
    metric::Kerr{T},
    roots::NTuple{4, T} where T
) -> NTuple{4, Any}

```


Returns the radial integrals for the case where there are no real roots in the radial potential


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1253" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.It_w_I0_terms' href='#Krang.It_w_I0_terms'><span class="jlbinding">Krang.It_w_I0_terms</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
  
- `νr` : Radial emission direction (Only necessary for case 1&amp;2 geodesics)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L944" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Ir_inf' href='#Krang.Ir_inf'><span class="jlbinding">Krang.Ir_inf</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Ir_inf(pix::Krang.AbstractPixel) -> Any

```


```
Ir_inf(pix::AbstractPixel)
```


Calculate the Ir infinity value for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The Ir infinity value of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L147" target="_blank" rel="noreferrer">source</a></Badge>



```julia
Ir_inf(metric::Kerr{T}, roots) -> Any

```


Returns the antiderivative $I_r=\int\frac{dr}{\sqrt{\mathcal{R(r)}}}$ evaluated at infinity. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `roots`  : Roots of the radial potential
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L313" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.emission_azimuth' href='#Krang.emission_azimuth'><span class="jlbinding">Krang.emission_azimuth</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/emission_coordinates.jl#L119" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.ConeGeometry' href='#Krang.ConeGeometry'><span class="jlbinding">Krang.ConeGeometry</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct ConeGeometry{T, A} <: Krang.AbstractGeometry
```


Cone Geometry with half opening angle `opening_angle`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/geometries/geometry_types.jl#L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Iϕ' href='#Krang.Iϕ'><span class="jlbinding">Krang.Iϕ</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Iϕ(pix::Krang.AbstractPixel, rs, τ, νr) -> Any

```


Returns the antiderivative $I_ϕ=\int\frac{a(2Mr-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `pix`: SlowLightIntensityPixel
  
- `rs` : Emission radius
  
- `τ` : Mino time
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1540" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.SlowLightIntensityCamera' href='#Krang.SlowLightIntensityCamera'><span class="jlbinding">Krang.SlowLightIntensityCamera</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct SlowLightIntensityCamera{T, A} <: Krang.AbstractCamera
```


Camera that caches slow light raytracing information for an observer sitting at radial infinity. The frame of this observer is alligned with the Boyer-Lindquist frame.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/SlowLightIntensityCamera.jl#L167" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Iϕ_w_I0_terms' href='#Krang.Iϕ_w_I0_terms'><span class="jlbinding">Krang.Iϕ_w_I0_terms</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Iϕ_w_I0_terms(metric::Kerr{T}, rs, τ, roots, νr, λ) -> Any

```


Returns the antiderivative $I_ϕ=\int\frac{a(2Mr-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$ with full I0 terms.

See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `τ`: Minotime 
  
- `roots` : Radial roots
  
- `νr` : Radial emission direction (Only necessary for case 1&amp;2 geodesics)
  
- `λ`  : Reduced azimuthal angular momentum
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L603" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.SlowLightIntensityPixel' href='#Krang.SlowLightIntensityPixel'><span class="jlbinding">Krang.SlowLightIntensityPixel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct SlowLightIntensityPixel{T} <: Krang.AbstractPixel{T}
```


Intensity Pixel Type.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/SlowLightIntensityCamera.jl#L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Iϕ_inf' href='#Krang.Iϕ_inf'><span class="jlbinding">Krang.Iϕ_inf</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Iϕ_inf(pix::Krang.AbstractPixel) -> Any

```


```
Iϕ_inf(pix::AbstractPixel)
```


Calculate the Iϕ infinity terms for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The Iϕ infinity terms of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L236" target="_blank" rel="noreferrer">source</a></Badge>



```julia
Iϕ_inf(metric::Kerr{T}, roots, λ) -> Any

```


Returns the antiderivative $I_ϕ=\int\frac{a(2Mr-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$ evaluated at infinity without I0 terms.

See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `roots` : Radial roots
  
- `λ`  : Reduced azimuthal angular momentum
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L462" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.radial_inf_integrals_case3' href='#Krang.radial_inf_integrals_case3'><span class="jlbinding">Krang.radial_inf_integrals_case3</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
radial_inf_integrals_case3(
    metric::Kerr{T},
    roots::NTuple{4, T} where T
) -> NTuple{4, Any}

```


Returns the radial integrals for the case where there are two real roots in the radial potential


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1210" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.AbstractPixel' href='#Krang.AbstractPixel'><span class="jlbinding">Krang.AbstractPixel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type AbstractPixel{T}
```


Abstract Pixel Type


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.radial_w_I0_terms_integrals_case2' href='#Krang.radial_w_I0_terms_integrals_case2'><span class="jlbinding">Krang.radial_w_I0_terms_integrals_case2</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1322" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.SlowLightIntensityScreen' href='#Krang.SlowLightIntensityScreen'><span class="jlbinding">Krang.SlowLightIntensityScreen</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct SlowLightIntensityScreen{T, A<:(AbstractMatrix)} <: Krang.AbstractScreen
```


Screen made of `SlowLightIntensityPixel`s.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/SlowLightIntensityCamera.jl#L89" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.radial_w_I0_terms_integrals_case4' href='#Krang.radial_w_I0_terms_integrals_case4'><span class="jlbinding">Krang.radial_w_I0_terms_integrals_case4</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
radial_w_I0_terms_integrals_case4(
    metric::Kerr{T},
    rs,
    roots::NTuple{4, T} where T,
    τ
) -> NTuple{4, Any}

```


Returns the radial integrals for the case where there are no real roots in the radial potential


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1448" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.radial_integrals' href='#Krang.radial_integrals'><span class="jlbinding">Krang.radial_integrals</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1579" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.AbstractCamera' href='#Krang.AbstractCamera'><span class="jlbinding">Krang.AbstractCamera</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type AbstractCamera
```


Abstract Observer Type


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L1" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.IntensityScreen' href='#Krang.IntensityScreen'><span class="jlbinding">Krang.IntensityScreen</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct IntensityScreen{T, A<:(AbstractMatrix)} <: Krang.AbstractScreen
```


Screen made of `IntensityPixel`s.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/IntensityCamera.jl#L65" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.It' href='#Krang.It'><span class="jlbinding">Krang.It</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
It(pix::Krang.AbstractPixel, rs, τ, νr) -> Any

```


Returns the antiderivative $I_t=\int\frac{a(2Mr-a\lambda)}{\sqrt{\Delta\mathcal{R(r)}}}dr$. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `pix`: SlowLightIntensityPixel
  
- `rs` : Emission radius
  
- `τ` : Mino time
  
- `νr` : Sign of radial velocity direction at emission. This is always positive for case 3 and case 4 geodesics.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1560" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.IntensityPixel' href='#Krang.IntensityPixel'><span class="jlbinding">Krang.IntensityPixel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct IntensityPixel{T} <: Krang.AbstractPixel{T}
```


Intensity Pixel Type.  Each Pixel is associated with a single ray, and caches some information about the ray.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/IntensityCamera.jl#L2" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.AbstractMaterial' href='#Krang.AbstractMaterial'><span class="jlbinding">Krang.AbstractMaterial</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type AbstractMaterial
```


Abstract Material


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/materials/AbstractMaterialTypes.jl#L1" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.AbstractScreen' href='#Krang.AbstractScreen'><span class="jlbinding">Krang.AbstractScreen</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type AbstractScreen
```


Abstract Screen Type


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.radial_w_I0_terms_integrals_case3' href='#Krang.radial_w_I0_terms_integrals_case3'><span class="jlbinding">Krang.radial_w_I0_terms_integrals_case3</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
radial_w_I0_terms_integrals_case3(
    metric::Kerr{T},
    rs,
    roots::NTuple{4, T} where T,
    τ
) -> NTuple{4, Any}

```


Returns the radial integrals for the case where there are two real roots in the radial potential


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1390" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.IntensityCamera' href='#Krang.IntensityCamera'><span class="jlbinding">Krang.IntensityCamera</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct IntensityCamera{T, A} <: Krang.AbstractCamera
```


Camera that caches fast light raytracing information for an observer sitting at radial infinity. The frame of this observer is alligned with the Boyer-Lindquist frame.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/IntensityCamera.jl#L145" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.Ir_s' href='#Krang.Ir_s'><span class="jlbinding">Krang.Ir_s</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Ir_s(metric::Kerr{T}, rs, roots, νr) -> Any

```


Returns the antiderivative $I_r=\int\frac{dr}{\sqrt{\mathcal{R(r)}}}$. See [`r_potential(x)`](/api#Krang.r_potential) for an implementation of $\mathcal{R}(r)$.

**Arguments**
- `metric`: Kerr{T} metric
  
- `rs` : Emission radius
  
- `roots`  : Roots of the radial potential
  
- `νr` : Radial emission direction (Only necessary for case 1&amp;2 geodesics)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L385" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.AbstractGeometry' href='#Krang.AbstractGeometry'><span class="jlbinding">Krang.AbstractGeometry</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type AbstractGeometry
```


Abstract Geometry


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/geometries/geometry_types.jl#L2" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.radial_inf_integrals_case2' href='#Krang.radial_inf_integrals_case2'><span class="jlbinding">Krang.radial_inf_integrals_case2</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
radial_inf_integrals_case2(
    metric::Kerr{T},
    roots::NTuple{4, T} where T
) -> NTuple{4, Any}

```


Returns the radial integrals for the case where there are four real roots in the radial potential, with roots outside the horizon.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L1161" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.It_inf' href='#Krang.It_inf'><span class="jlbinding">Krang.It_inf</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
It_inf(pix::Krang.AbstractPixel) -> Any

```


```
It_inf(pix::AbstractPixel)
```


Calculate the It infinity terms for a given pixel.

**Arguments**
- `pix::AbstractPixel`: The pixel of a screen.
  

**Returns**
- The It infinity terms of the pixel.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/cameras/camera_types.jl#L251" target="_blank" rel="noreferrer">source</a></Badge>



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
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/metrics/Kerr/misc.jl#L751" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild' href='#Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild'><span class="jlbinding">Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
boyer_lindquist_to_quasi_cartesian_kerr_schild(
    metric::Kerr,
    tBL,
    rBL,
    θBL,
    ϕBL
) -> NTuple{4, Any}

```


```
Transforms from Boyer-Lindquist to Kerr-Schild coordinates
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/schemes/RayTrace.jl#L206" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light' href='#Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light'><span class="jlbinding">Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(
    metric::Kerr{T},
    rBL,
    θBL,
    ϕBL
) -> Any

```


```
Transforms from Boyer-Lindquist to Kerr-Schild coordinates ignoring time
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/dominic-chang/Krang.jl/blob/787e8de7bb6cb24071bd4ad077130c4d3f159686/src/schemes/RayTrace.jl#L220" target="_blank" rel="noreferrer">source</a></Badge>

</details>

