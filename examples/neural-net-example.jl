# # Neural Network Emission Model Example
# This is a pedagogical example that serves as a proof of concept.
# We will build a simple General Relativistic Neural Radiance Field (NeRF) model to be ray traced and optimize with the ADAM algorithm.

# ##Setup 
# We will first import the necessary packages and set the random seed.
# Our emission model will be a neural network built with Lux.jl
using Lux
using Krang
using Random
using Accessors
Random.seed!(123)
rng = Random.GLOBAL_RNG

# Our model will take in spacetime coordinates and return observed intensity value for a given pixel:
# We will use 0.99 spin Kerr metric with an observer sitting at 20 degrees inclination with respect to the spin axis in this example.
# These parameters are fixed for this example, but could be made to vary in the optimization process.

# Lets define an `KerrNeRF` which will be comprised of an emission layer that we will raytrace.
# We will do this by first creating a struct to represent our image model that will store our emission model as a layer.

struct KerrNeRF{T<:Chain}
    emission_layer::T
end

# The models in Lux are functors that take in features, parameters and model state, and return the output and model state.
# Lets define the function associated with our `KerrNeRF` type.
# We will assume that the emission is coming from emission that originates in the equatorial plane.
function (m::KerrNeRF)(x, ps, st)
    metric = Krang.Kerr(ps.spin)
    θo = ps.θo * π
    pixels = Krang.SlowLightIntensityPixel.(Ref(metric), x[1, :], x[2, :], θo)

    sze = unsafe_trunc(Int, sqrt(size(x)[2]))
    coords = zeros(Float64, 2, sze * sze)
    redshifts = zeros(Float64, 1, sze * sze)
    emission_vals = zeros(Float64, 1, sze * sze)
    for n = 0:1
        for i = 1:sze
            for j = 1:sze
                pix = pixels[i+(j-1)*sze]
                α, β = Krang.screen_coordinate(pix)
                ηtemp = η(metric, α, β, θo)
                λtemp = λ(metric, α, θo)
                rs, ϕs, νr, νθ, _ =
                    Krang.emission_coordinates_fast_light(pix, Float64(π / 2), β > 0, n)

                T = typeof(α)
                if rs ≤ Krang.horizon(metric)
                    coords[1, i+(j-1)*sze] = zero(T)
                    coords[2, i+(j-1)*sze] = zero(T)
                    redshifts[(i-1)*sze+j] = zero(T)
                else
                    ϕks = Krang.ϕ_kerr_schild(metric, rs, ϕs)
                    xs = rs * cos(ϕks)
                    ys = rs * sin(ϕks)

                    curr_p_bl_d = p_bl_d(metric, rs, π / 2, ηtemp, λtemp, νr, νθ)
                    curr_p_bl_u = metric_uu(metric, rs, π / 2) * curr_p_bl_d
                    p_zamo_u = jac_zamo_u_bl_d(metric, rs, π / 2) * curr_p_bl_u
                    coords[1, (i-1)*sze+j] = xs
                    coords[2, (i-1)*sze+j] = ys
                    redshifts[(i-1)*sze+j] = inv(p_zamo_u[1])
                end

            end
        end
        emission_vals .+= m.emission_layer(coords, ps, st)[1] .* redshifts
    end
    emission_vals, st
end

# Lets define an emission layer for our model as a simple fully connected neural network with 2 hidden layers.
# The emission layer will take in 2D coordinates on an equatorial disk in the bulk spacetime and return a scalar intensity value.

emission_model = Chain(
    Dense(2 => 20, Lux.tanh_fast),
    Dense(20 => 20, Lux.tanh_fast),
    Dense(20 => 1, Lux.sigmoid),
)

ps, st = Lux.setup(rng, emission_model); # Get the emission model parameters and state
ps = @insert ps.spin = 0.94 # Set the spin of the black hole
ps = @insert ps.θo = 20.0 / 180 # Set the inclination angle of the observer

# We can now create an image model with our emission layer.
image_model = KerrNeRF(emission_model)

## Plotting the model

# Lets create an 20x20 pixel image of the `image_model` with a field of view of $10 MG/c^2$.

sze = 30
ρmax = 10e0
pixels = zeros(Float64, 2, sze * sze)
for (iiter, i) in enumerate(range(-ρmax, ρmax, sze))
    for (jiter, j) in enumerate(range(-ρmax, ρmax, sze))
        pixels[2, iiter+(jiter-1)*sze] = Float64(i)
        pixels[1, iiter+(jiter-1)*sze] = Float64(j)
    end
end

# We can see the effects of raytracing on emission in the bulk spacetime by plotting an image of the emission model and the image model.
using CairoMakie
curr_theme = Theme(
    Axis = (
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
    ),
    Heatmap = (colormap = :afmhot,),
)
set_theme!(merge(curr_theme, theme_latexfonts()))

emitted_intensity = reshape(emission_model(pixels, ps, st)[1], sze, sze)
received_intensity = reshape(image_model(pixels, ps, st)[1], sze, sze)

fig = Figure();
heatmap!(Axis(fig[1, 1], aspect = 1, title = "Emission Model"), emitted_intensity)
heatmap!(
    Axis(fig[1, 2], aspect = 1, title = "Image Model (Lensed Emission Model)"),
    received_intensity,
)
CairoMakie.save("emission_model_and_target_model.png", fig)

# ![image](emission_model_and_target_model.png)

# ## Fitting the NeRF model
# This will be a toy example showing the mechanics of fitting our KerrNeRF to a target image using the normalized cross correlation as a kernel for our loss function.
# This will be the image we will try to fit our model to.
target_img = reshape(received_intensity, 1, sze * sze);

# Lets fit our model using the normalized cross correlation as a kernel for our loss function.

using Enzyme
using Optimization
using OptimizationOptimisers
using StatsBase
using ComponentArrays

function mse(img1::Matrix{T}, img2::Matrix{T}) where {T}
    mean(((img1 ./ sum(img1)) .- (img2 ./ sum(img2))) .^ 2)
end

function loss_function(pixels, y, ps, st)
    y_pred, st = image_model(pixels, ps, st)
    mse(y, y_pred), st
end

mse(target_img, target_img)

ps, st = Lux.setup(rng, emission_model);
ps = @insert ps.spin = 0.5
ps = @insert ps.θo = 60.0 / 180
image_model = KerrNeRF(emission_model);

emitted_intensity = reshape(emission_model(pixels, ps, st)[1], sze, sze)
received_intensity = reshape(image_model(pixels, ps, st)[1], sze, sze)
loss_function(pixels, target_img, ps, st)

fig = Figure();
heatmap!(
    Axis(fig[1, 1], aspect = 1, title = "Emission Model"),
    emitted_intensity,
    colormap = :afmhot,
)
heatmap!(
    Axis(fig[1, 2], aspect = 1, title = "Imgage Model (Lensed Emission Model)"),
    received_intensity,
    colormap = :afmhot,
)
CairoMakie.save("emission_model_and_image_model.png", fig)

# ![image](emission_model_and_image_model.png)

# Lets define callback function to print the loss as we optimize our model.
mutable struct Callback
    counter::Int
    stride::Int
    const f::Function
end
Callback(stride, f) = Callback(0, stride, f)
function (c::Callback)(state, loss, others...)
    c.counter += 1
    if c.counter % c.stride == 0
        @info "On step $(c.counter) loss = $(loss)"
        return false
    else
        return false
    end
end

# We can now optimize our model using the ADAM optimizer.
ps_trained, st_trained = let st = Ref(st), x = pixels, y = reshape(target_img, 1, sze * sze)

    optprob = Optimization.OptimizationProblem(
        Optimization.OptimizationFunction(
            function (ps, constants)
                loss, st[] = loss_function(x, y, ps, st[])
                loss
            end,
            Optimization.AutoEnzyme(),
        ),
        ComponentArrays.ComponentVector{Float64}(ps),
    )

    solution = Optimization.solve(
        optprob,
        OptimizationOptimisers.Adam(),
        maxiters = 25_000,
        callback = Callback(100, () -> nothing),
    )

    solution.u, st[]
end

# Let's plot the results of our optimization. and compare it to the target image.
received_intensity, st =
    ((x) -> (reshape(x[1], sze, sze), x[2]))(image_model(pixels, ps_trained, st_trained))
acc_intensity, st = ((x) -> (reshape(x[1], sze, sze), x[2]))(image_model(pixels, ps, st))
loss_function(pixels, target_img, ps, st)
loss_function(pixels, target_img, ps_trained, st_trained)

using Printf

fig = Figure(size = (700, 300))
heatmap!(Axis(fig[1, 1], aspect = 1, title = "Target Image"), reshape(target_img, sze, sze))
heatmap!(
    Axis(
        fig[1, 2],
        aspect = 1,
        title = "Starting State (loss=$(@sprintf("%0.2e", loss_function(pixels, target_img, ps, st)[1])))",
    ),
    acc_intensity,
)
heatmap!(
    Axis(
        fig[1, 3],
        aspect = 1,
        title = "Fitted State (loss=$(@sprintf("%0.2e", loss_function(pixels, target_img, ps_trained, st_trained)[1])))",
    ),
    received_intensity,
)
fig

save("neural_net_results.png", fig)

# ![image](neural_net_results.png)
