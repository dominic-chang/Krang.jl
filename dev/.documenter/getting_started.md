
# What is Krang.jl? {#What-is-Krang.jl?}

Kerr Raytracer for Analytic Null Geodesics (Krang) is a raytracing code for geometries that are embedded in the Kerr spacetime. Krang solves the Kerr geodesic problem using the analytic solutions that have been derived from the Hamilton-Jacobi Formalism. This choice makes the code efficient and accurate, and is an ideal formalism to isolated sub images that manifest from strong lensing.

## Philosophy

The current philosophy of Krang is to raytrace meshes. Meshes are made by &#39;painting&#39; materials onto geometries. These are then viewed by cameras that are sensitive to particular observables (_e.g. intensity and polarization_).

## Raytracing

Light rays in this module can be parameterized in terms of either the cones ($\theta_s$), or the minotime ($\Delta\tau$). Parameterization interms of cones allows for images to separated interms of sub images and ray traced individually.

### Raytracing conical surfaces {#Raytracing-conical-surfaces}

Surfaces of constant $\theta$ define spin axis centered cones whose apex lie at the origin of the coordinate system.


![](examples/coordinate.gif)

<p style="text-align:center">n=0 and n=1 images of emission coordinates originating from conical surfaces.</p>


### Raytracing with rays parameterized by Mino time {#Raytracing-with-rays-parameterized-by-Mino-time}


![](examples/raytrace.gif)

<p style="text-align:center">Coordinate evolution with Mino time.</p>

