## Raytracing

Light rays in this module can be parameterized in terms of either the emission inclination ($\theta_s$), or the Mino time ($\Delta\tau$).
Parameterization in terms of emission inclination allows for images to be divided into sub images which are ray traced individually.

### Raytracing conical surfaces
Surfaces of constant $\theta_s$ define spin axis centered cones whose apex lie at the origin of the Boyer-Lindquist coordinate system.

See the [conical-surface example](examples/coordinate-example.md) for the $n=0$ and $n=1$ images of emission coordinates.

### Raytracing with rays parameterized by Mino time
Mino time, $\tau$, is a parameter monotonic in affine parameter, $\tau'$, defined by

```math
d\tau = \frac{d\tau'}{\Sigma(r,\theta)},
```

where

```math
\Sigma(r,\theta) = r^2 +a^2\cos^2\theta.
```

See the [Mino-time example](examples/mino-time-example.md) for coordinate evolution along a ray.

### Cameras
Cameras cache pre-computed information that is constant for a given camera location. 
There are currently two types of cameras which can be used for either 'slow light' or 'fast light' raytracing.

* `IntensityCamera` : Pre-computes geodesic information necessary to solve the 'fast light' raytracing problem.

* `SlowLightIntensityCamera` : Pre-computes geodesic information necessary to solve the 'slow light' raytracing problem.

With `KernelAbstractions` and a GPU backend loaded, construct a screen by passing the GPU array type as the final positional argument. For example, with CUDA:

```julia
using CUDA, KernelAbstractions
 
screen = Krang.SlowLightIntensityScreen(
    metric, -ρmax, ρmax, -ρmax, ρmax, θo, sze, CuArray
)
store = Krang.render.(screen.pixels, Ref(scene))
```
