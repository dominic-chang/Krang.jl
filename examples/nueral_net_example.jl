# # Nueral Network Emission Model Example
# In this example we will build a simple Neural Network emission model which we will ray trace and optimize with LBFGS
using Pkg; Pkg.activate((@__DIR__))    
using Lux
using Krang
using Random
Random.seed!(123)
rng = Random.GLOBAL_RNG

# Our model with take in spacetime coordinates and return an intensity value:
# This is a simple fully connected Neural Network with 2 hidden layers

emission_model = Chain(
    Dense(2 => 20, Lux.sigmoid),  
    Dense(20 => 20, Lux.sigmoid),  
    Dense(20 => 1, Lux.sigmoid)
    ) 

ps, st = Lux.setup(rng, emission_model)

# Our image model will effectively act as a raytracing layer after our emission model.
# In our example, we will use 0.99 spin Kerr metric with an observer sitting at 20 degrees inclination with respect to
# the spin axis.
# We will create a struct for our image model and store our emission model as an argument.

struct ImageModel{T <: Chain}
    emission_layer::T
end

function (m::ImageModel)(x, ps, st)
    metric = Krang.Kerr(0.99e0)
    θo = Float64(20/180*π)
    pixels = Krang.IntensityPixel.(Ref(metric), x[1,:], x[2,:], θo)

    sze = unsafe_trunc(Int, sqrt(size(x)[2]))
    coords = zeros(Float64, 2,sze*sze)
    emission_vals = zeros(Float64, 1, sze*sze)
    for n in 0:1
        for i in 1:sze
            for j in 1:sze
                pix = pixels[i+(j-1)*sze]
                α, β = Krang.screen_coordinate(pix)
                T = typeof(α)
                rs, _, ϕs = Krang.emission_coordinates_fast_light(pix, Float64(π/2), β > 0, n)[1:3]
                xs = rs * cos(ϕs)
                ys = rs * sin(ϕs)
                if hypot(xs, ys) ≤ Krang.horizon(metric)
                    coords[1,i+(j-1)*sze] = zero(T)
                    coords[2,i+(j-1)*sze] = zero(T)
                else
                    coords[1,i+(j-1)*sze] = xs
                    coords[2,i+(j-1)*sze] = ys
                end

            end
        end
        emission_vals .+= m.emission_layer(coords, ps, st)[1]
    end
    emission_vals,st 
end

# Lets create an 20x20 pixel image with a field of view of $10 MG/c^2$.

image_model = ImageModel(emission_model)
sze = 20
ρmax = 10e0
pixels = zeros(Float64, 2, sze*sze)
for (iiter, i) in enumerate(range(-ρmax, ρmax, sze))
    for (jiter, j) in enumerate(range(-ρmax, ρmax, sze))
        pixels[1,iiter+(jiter-1)*sze] = Float64(i)
        pixels[2,iiter+(jiter-1)*sze] = Float64(j)
    end
end

# Lets see what our emission model looks like before and after raytracing.
using CairoMakie
curr_theme = Theme(
    Axis = (xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,),
    Heatmap =(colormap=:afmhot, ),
)
set_theme!(merge(curr_theme, theme_latexfonts()))

emitted_intensity = reshape(emission_model(pixels, ps, st)[1], sze, sze)
received_intensity = reshape(image_model(pixels, ps, st)[1], sze, sze)

fig = Figure();
heatmap!(Axis(fig[1,1], aspect=1, title="Emission Model"), emitted_intensity)
heatmap!(Axis(fig[1,2], aspect=1, title="Image Model (Lensed Emission Model)"), received_intensity)
save("emission_model_and_target_model.png", fig)

# ![image](emission_model_and_target_model.png)

# This will be the image we will try to fit our model to.
target_img = reshape(received_intensity, sze, sze)

# ## Fitting the model
# Lets fit our model using the normalized cross correlation as a kernel for our loss function.

using Enzyme
using Optimization 
using OptimizationOptimisers
using StatsBase
using ComponentArrays
Enzyme.API.runtimeActivity!(true)

Enzyme.Compiler.RunAttributor[] = false

function bhattacharyya(img1::Matrix{T}, img2::Matrix{T}) where T
    img1 = reshape(img1, sze, sze) ./ sum(img1)
    img2 = reshape(img2, sze, sze) ./ sum(img2)
    sum(sqrt.(img1 .* img2))
end

function loss_function(pixels, y, ps, st)
    y_pred, st = image_model(pixels, ps, st)
    -log(bhattacharyya(y, y_pred)), st
end

bhattacharyya(target_img, target_img)

ps, st = Lux.setup(rng, emission_model)
image_model = ImageModel(emission_model)

emitted_intensity = reshape(emission_model(pixels, ps, st)[1], sze, sze)
received_intensity = reshape(image_model(pixels, ps, st)[1], sze, sze)
loss_function(pixels, target_img, ps, st)

fig = Figure();
heatmap!(Axis(fig[1,1], aspect=1, title="Emission Model"), emitted_intensity, colormap=:afmhot)
heatmap!(Axis(fig[1,2], aspect=1, title="Imgage Model (Lensed Emission Model)"), received_intensity, colormap=:afmhot)
save("emission_model_and_image_model.png", fig)

# ![image](emission_model_and_image_model.png)


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

ps_trained, st_trained = let st=Ref(st), x=pixels, y=target_img
    
    optprob = Optimization.OptimizationProblem(
        Optimization.OptimizationFunction(
            function(ps, constants)
                loss, st[] = loss_function(x, y, ps, st[])
                loss
            end,
            Optimization.AutoEnzyme()
        ),
        ComponentArrays.ComponentVector{Float64}(ps)
    )
    
    solution = Optimization.solve(
        optprob,
        OptimizationOptimisers.ADAM(0.05),
        maxiters = 1_000, 
    callback=Callback(10,()->nothing)
    )
    
    solution.u, st[]
end

received_intensity, st = ((x) -> (reshape(x[1], sze, sze), x[2]))(image_model(pixels, ps_trained, st_trained))
acc_intensity, st = ((x) -> (reshape(x[1], sze, sze), x[2]))(image_model(pixels, ps, st))
loss_function(pixels, target_img, ps, st)
loss_function(pixels, target_img, ps_trained, st_trained)
using Printf
begin 
    fig = Figure(size=(700, 300));
    heatmap!(Axis(fig[1,1], aspect=1, title="Target Image"), target_img)
    heatmap!(Axis(fig[1,2], aspect=1, title="Starting State (loss=$(@sprintf("%0.2e", loss_function(pixels, target_img, ps, st)[1])))"), acc_intensity)  
    heatmap!(Axis(fig[1,3], aspect=1, title="Fitted State (loss=$(@sprintf("%0.2e", loss_function(pixels, target_img, ps_trained, st_trained)[1])))"), received_intensity)
    save("neural_net_results.png", fig)
end

# ![image](neural_net_results.png)

