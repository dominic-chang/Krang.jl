# What is Krang.jl?

Kerr Raytracer for Analytic Null Geodesics (Krang) is a raytracing code for geometries that are embedded in the Kerr spacetime.
Krang solves the Kerr geodesic problem using analytic solutions derived from the Hamilton-Jacobi formalism. 
This choice makes the code efficient and accurate, and is ideal for decomposing images into the characteristic sub-images which manifest from strong gravitational lensing. 

## Philosophy

Krang operates by ray tracing 'meshes' which can be thought of as a geometry that presents some sort of physics material. 
The meshes are viewed by cameras that are sensitive to particular observables (*e.g.intensity and polarization*).

## Examples
Some pedagogical example usage of the raytracing code in the [Examples](examples/coordinate-example.md) section. These examples primarily use `CairoMakie` and `GLMakie` for plotting.
You can install these by enering the `julia` package mode and running:

```julia
pkg> add CairoMakie GLMakie
```

