## Raytracing

Light rays in this module can be parameterized in terms of either the emission inclination ($\theta_s$), or the Mino time ($\Delta\tau$).
Parameterization in terms of emission inclination allows for images to be divided into sub images which are ray traced individually.

### Raytracing conical surfaces
Surfaces of constant $\theta_s$ define spin axis centered cones whose apex lie at the origin of the Boyer-Lindquist coordinate system.

![raytracing conical surfaces](examples/coordinate.gif)
```@raw html
<p style="text-align:center">n=0 and n=1 images of emission coordinates originating from conical surfaces.</p>
```

### Raytracing with rays parameterized by Mino time
Mino time, $\tau$, is a parameter monotonic in affine parameter, $\tau'$, defined by

```math
d\tau = \Sigma(r,\theta)d\tau',
```

where

```math
\Sigma(r,\theta) = r^2 +a^2\cos^2\theta.
```

![raytracing with Mino time](examples/raytrace.gif)
```@raw html
<p style="text-align:center">Coordinate evolution with Mino time.</p>
```

### Cameras
Cameras cache pre-computed information that is constant for a given camera location. 
There are currently two types of cameras which can be used for either 'slow light' or 'fast light' raytracing.

* `IntensityCamera` : Pre-computes geodesic information necessary to solve the 'fast light' raytracing problem.

* `SlowLightIntensityCamera` : Pre-computes geodesic information necessary to solve the 'slow light' raytracing problem.

The GPU arrays can be passed to the cameras on construction to raytrace enforce raytracing on the GPU.
A sketch of how to do this with a CUDA array is:

```julia
using CUDA
 
store = CUDA.fill(0.0, sze, sze)
camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze, A=CuArray)
Krang.render!(store, camera, scene)
```
