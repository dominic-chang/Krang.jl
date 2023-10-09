using CairoMakie
using Kang
# # Equatorial
# Equatorial is an analytic model for axi-symmetric synchrotron emission around a Kerr Blackhole. The Emission region is
# defined by 2 cones whose emission intensity is described by a local fluid velocity, and magnetic field orienation.
# The tempered by a double power law enveloping function.
#
# We will begin by importing `Kang`, where our model is defined, and `CairoMakie`, for plotting.

spin = -0.5
model = Kang.Kerr(spin)
sze = 100
θo = π/4
rmax = 20
ρmax = 20
rmin = Kang.horizon(model)

function nan2num(x::T) where T 
    isnan(x) ? zero(T) : x
end

fig = Figure();

#record(fig, "emission_coordinates.mp4", range(0,π,length=50),framerate=12) do θs
for θs in range(0,π,length=100)
    empty!(fig)
    #for θs in range(0,π,length=100)

    # We will display the emission stokes parameters of the model as a set of 200 x 200 pixel images.

    αvals = range(-ρmax, ρmax, length=sze)
    βvals = range(-ρmax, ρmax, length=sze)
    tvals = zeros(sze, sze)
    rvals = zeros(sze, sze)
    θvals = zeros(sze, sze)
    ϕvals = zeros(sze, sze)

    #n = 1
    
    #Threads.@threads for i in 1:sze
    #    for j in 1:sze
    #        curr_ts, curr_rs, curr_θs, curr_ϕs = nan2num.(Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, true, n))# .+ nan2num.(Kang.emission_coordinates(model, αvals[i], βvals[j], false, n))
    #        tvals[i, j] = (rmax > curr_rs > rmin) ? curr_ts : tvals[i, j]
    #        rvals[i, j] = (rmax > curr_rs > rmin) ? curr_rs : rvals[i, j]
    #        θvals[i, j] = (rmax > curr_rs > rmin) ? curr_θs : θvals[i, j]
    #        ϕvals[i, j] = (rmax > curr_rs > rmin) ? curr_ϕs : ϕvals[i, j]

    #        curr_ts, curr_rs, curr_θs, curr_ϕs = nan2num.(Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, false, n))# .+ nan2num.(Kang.emission_coordinates(model, αvals[i], βvals[j], false, n))
    #        tvals[i, j] = (rmax > curr_rs > rmin && curr_ts != 0) ? curr_ts : tvals[i, j]
    #        rvals[i, j] = (rmax > curr_rs > rmin && curr_rs > 0) ? curr_rs : rvals[i, j]
    #        θvals[i, j] = (rmax > curr_rs > rmin && curr_θs > 0) ? curr_θs : θvals[i, j]
    #        ϕvals[i, j] = (rmax > curr_rs > rmin && curr_ϕs > 0) ? curr_ϕs : ϕvals[i, j]
    #    end
    #end

    n = 0
    Threads.@threads for i in 1:sze
    #for i in 1:sze
        for j in 1:sze
            curr_ts, curr_rs, curr_θs, curr_ϕs = nan2num.(Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, true, n))# .+ nan2num.(Kang.emission_coordinates(model, αvals[i], βvals[j], false, n))
            tvals[i, j] = (rmax > curr_rs > rmin) ? curr_ts : tvals[i, j]
            rvals[i, j] = (rmax > curr_rs > rmin) ? curr_rs : rvals[i, j]
            θvals[i, j] = (rmax > curr_rs > rmin) ? curr_θs : θvals[i, j]
            ϕvals[i, j] = (rmax > curr_rs > rmin) ? curr_ϕs : ϕvals[i, j]
            curr_ts, curr_rs, curr_θs, curr_ϕs = nan2num.(Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, false, n))# .+ nan2num.(Kang.emission_coordinates(model, αvals[i], βvals[j], false, n))
            tvals[i, j] = (rmax > curr_rs > rmin && curr_ts != 0) ? curr_ts : tvals[i, j]
            rvals[i, j] = (rmax > curr_rs > rmin && curr_rs > 0) ? curr_rs : rvals[i, j]
            θvals[i, j] = (rmax > curr_rs > rmin && curr_θs > 0) ? curr_θs : θvals[i, j]
            ϕvals[i, j] = (rmax > curr_rs > rmin && curr_ϕs > 0) ? curr_ϕs : ϕvals[i, j]

        end
    end

    # Lets plot the resulting image

    ax0 = Axis(fig[1, 1], aspect=1, title="time");
    ax1 = Axis(fig[1, 3], aspect=1, title="radius");
    ax2 = Axis(fig[2, 1], aspect=1, title="opening angle");
    ax3 = Axis(fig[2, 3], aspect=1, title="azimuth");

    hm0 = CairoMakie.heatmap!(ax0, tvals, colormap=:afmhot, colorrange=(-20,100));
    cb0 = Colorbar(fig[1, 2], hm0);
    hm1 = CairoMakie.heatmap!(ax1, rvals, colormap=:afmhot);
    cb1 = Colorbar(fig[1, 4], hm1);
    hm2 = CairoMakie.heatmap!(ax2, θvals, colormap=:afmhot)
    cb2 = Colorbar(fig[2, 2], hm2);
    hm3 = CairoMakie.heatmap!(ax3, ϕvals, colormap=:afmhot);
    cb3 = Colorbar(fig[2, 4], hm3);
    
    display(fig)
end

import DisplayAs
showable("image/png", fig)
img = DisplayAs.PNG(fig)

# Notice that there is no Stokes V since synchrotron emission produces linear polarization.
