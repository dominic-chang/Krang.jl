# # Creating a Custom Dual Cone Model

# We will define a custom model for low luminosity active galactic nuclei (LLAGN).
# A detailed description of the model can be found in this [reference](https://doi.org/10.48550/arXiv.2405.04749).
# We will show the emission of the n=0 (direct) and n=1 (indirect) photons as they are emitted from the 
# source, at a fixed inclination angle from the black hole's spin axis.
#
# First, we will import Krang and specify the spacetime. We will use a Kerr black hole of 0.94 spin.
using Krang
metric = Krang.Kerr(0.94);

# Let's create a camera with a resolution of 400x400 pixels viewed by an asymptotic observer at an inclination angle of $θo=17^\circ$. 
# The camera will have a field of view of $10 MG/c^2$.
θo = 17 * π / 180;
ρmax = 10.0;
camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, 400);

# We will need to create `Mesh` objects to render the scene.
# First, we will create the material for the mesh.
# Our material will be the `ElectronSynchrotronPowerLawPolarization` material with the following parameters.

χ = -1.7;
ι = 0.58;
βv = 0.87;
σ = 0.73;
η1 = 2.64;
η2 = π - η1;

# These will be used to define the magnetic field and fluid velocity.

magfield1 = Krang.SVector(sin(ι) * cos(η1), sin(ι) * sin(η1), cos(ι));
magfield2 = Krang.SVector(sin(ι) * cos(η2), sin(ι) * sin(η2), cos(ι));
vel = Krang.SVector(βv, (π / 2), χ);
R = 4.0;
p1 = 1.0;
p2 = 4.0;

# Finally, we define the materials for each of the cones
material1 = Krang.ElectronSynchrotronPowerLawPolarization(
    magfield1...,
    vel...,
    σ,
    R,
    p1,
    p2,
    (0, 1),
);
material2 = Krang.ElectronSynchrotronPowerLawPolarization(
    magfield2...,
    vel...,
    σ,
    R,
    p1,
    p2,
    (0, 1),
);

# Next we will define the geometries of each mesh. We will use a `ConeGeometry` with an opening angle of $75^\circ$.

θs = (75 * π / 180);
geometry1 = Krang.ConeGeometry(θs)
geometry2 = Krang.ConeGeometry(π - θs)


# We will create two meshes, one for each geometry anc create a scene with both meshes.
mesh1 = Krang.Mesh(geometry1, material1)
mesh2 = Krang.Mesh(geometry2, material2)

# Finally, we will render the scene with the camera

scene = Krang.Scene((mesh1, mesh2))
stokesvals = render(camera, scene)

# We will import CairoMakie for plotting the results.
import CairoMakie as CMk

curr_theme = CMk.Theme(
    Axis = (
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
    ),
)
CMk.set_theme!(merge!(curr_theme, CMk.theme_latexfonts()))

fig = CMk.Figure(resolution = (700, 700));
ax1 = CMk.Axis(fig[1, 1], aspect = 1, title = "I")
ax2 = CMk.Axis(fig[1, 3], aspect = 1, title = "Q")
ax3 = CMk.Axis(fig[2, 1], aspect = 1, title = "U")
ax4 = CMk.Axis(fig[2, 3], aspect = 1, title = "V")
colormaps = [:afmhot, :redsblues, :redsblues, :redsblues]

hms =
    zip(
        [ax1, ax2, ax3, ax4],
        [
            getproperty.(stokesvals, :I),
            getproperty.(stokesvals, :Q),
            getproperty.(stokesvals, :U),
            getproperty.(stokesvals, :V),
        ],
        colormaps,
    ) .|> x -> CMk.heatmap!(x[1], x[2], colormap = x[3], rasterize = true)
CMk.Colorbar(fig[1, 2], hms[1], labelsize = 20)
CMk.Colorbar(fig[1, 4], hms[2], labelsize = 20)
CMk.Colorbar(fig[2, 2], hms[3], labelsize = 20)
CMk.Colorbar(fig[2, 4], hms[4], labelsize = 20)


fig

CMk.save("polarization_example.png", fig)

# ![polarization of emission model](polarization_example.png)

# > [!IMPORTANT]
# > You can also save the rendered image as a fits file for further analysis.
# 
# ```julia
# using Comrade
# grid = imagepixels(μas2rad(120), μas2rad(120), 400, 400)
# Comrade.save_fits(joinpath((@__DIR__), "polarized_models.fits"), IntensityMap(stokesvals,grid))
# ```
