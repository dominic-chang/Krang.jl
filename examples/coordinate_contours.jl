# # Raytracing with Kang.jl

# In this example, we will raytrace the region around a Kerr blackhole as seen by an observer stationed at infinity.
# We will show the emission coordinates of the n=0 (direct) and n=1 (indirect) photons as they are emitted from the 
# source, at a fixed inclination angle from the blackhole's spin axis.
#
# First, let's import Kang and CairoMakie for plotting.
using CairoMakie
using Kang
#
# We will use a 0.99 spin Kerr blackhole viewed by an assymptotic observer at an inclination angle of θo=π/4. 
# A region spanned by radii between the horizon and 20M at varying inclinations will be raytraced onto the 20Mx20M 
# screen of the observer.
model = Kang.Kerr(0.99);
θo = π/4;
sze = 100;
rmin = Kang.horizon(model)
rmax = 20;
ρmax = 20;

# Let us now create a figure to plot the emission coordinates on,
fig = Figure(resolution=(600, 600), fontfamily="Computer Modern", fontface="bold");
# and use this figure make an animation by looping over the inclination angle θs.
# This loop will plot the emission coordinates for each θs.
recording = record(fig, "emission_coordinates.gif", range(0,π,length=180),framerate=15) do θs
    empty!(fig);

    αvals = range(-ρmax, ρmax, length=sze);
    βvals = range(-ρmax, ρmax, length=sze);
    tvals = -Inf.*ones(sze, sze);
    rvals = -Inf.*ones(sze, sze);
    θvals = -Inf.*ones(sze, sze);
    ϕvals = -Inf.*ones(sze, sze);

    n = 1;
    Threads.@threads for i in 1:sze
        for j in 1:sze
            curr_ts, curr_rs, curr_θs, curr_ϕs, νr, νθ = Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, true, n);
            tvals[i, j] = (rmax > curr_rs > rmin) ? curr_ts : tvals[i, j];
            rvals[i, j] = (rmax > curr_rs > rmin) ? curr_rs : rvals[i, j];
            θvals[i, j] = (rmax > curr_rs > rmin) ? curr_θs : θvals[i, j];
            ϕvals[i, j] = (rmax > curr_rs > rmin) ? curr_ϕs : ϕvals[i, j];

            curr_ts, curr_rs, curr_θs, curr_ϕs, νr, νθ = Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, false, n);
            tvals[i, j] = (rmax > curr_rs > rmin && curr_ts != 0) ? curr_ts : tvals[i, j];
            rvals[i, j] = (rmax > curr_rs > rmin && curr_rs > 0) ? curr_rs : rvals[i, j];
            θvals[i, j] = (rmax > curr_rs > rmin && curr_θs > 0) ? curr_θs : θvals[i, j];
            ϕvals[i, j] = (rmax > curr_rs > rmin && curr_ϕs > 0) ? curr_ϕs : ϕvals[i, j];
        end
    end

    n = 0;
    Threads.@threads for i in 1:sze
        for j in 1:sze
            curr_ts, curr_rs, curr_θs, curr_ϕs, νr, νθ = Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, true, n);
            tvals[i, j] = (rmax > curr_rs > rmin) ? curr_ts : tvals[i, j];
            rvals[i, j] = (rmax > curr_rs > rmin) ? curr_rs : rvals[i, j];
            θvals[i, j] = (rmax > curr_rs > rmin) ? curr_θs : θvals[i, j];
            ϕvals[i, j] = (rmax > curr_rs > rmin) ? curr_ϕs : ϕvals[i, j];
            curr_ts, curr_rs, curr_θs, curr_ϕs, νr, νθ = Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, false, n);
            tvals[i, j] = (rmax > curr_rs > rmin && curr_ts != 0) ? curr_ts : tvals[i, j];
            rvals[i, j] = (rmax > curr_rs > rmin && curr_rs > 0) ? curr_rs : rvals[i, j];
            θvals[i, j] = (rmax > curr_rs > rmin && curr_θs > 0) ? curr_θs : θvals[i, j];
            ϕvals[i, j] = (rmax > curr_rs > rmin && curr_ϕs > 0) ? curr_ϕs : ϕvals[i, j];
        end
    end

    ax0 = Axis(fig[1, 1], aspect=1, title=L"\text{Radius increasing on emission}");
    ax1 = Axis(fig[1, 3], aspect=1, title=L"\text{inclination increasing on emission}");
    ax2 = Axis(fig[2, 1], aspect=1, title=L"\text{Opening angle } (\theta_s)");
    ax3 = Axis(fig[2, 3], aspect=1, title=L"\text{Azimuth } (\phi_s)");

    hm0 = CairoMakie.heatmap!(ax0, αvals, βvals, tvals, colormap=:afmhot, colorrange=(-100,100));
    cb0 = Colorbar(fig[1, 2], hm0);
    hm1 = CairoMakie.heatmap!(ax1, αvals, βvals, rvals, colormap=:afmhot, colorrange=(0,20));
    cb1 = Colorbar(fig[1, 4], hm1);
    hm2 = CairoMakie.heatmap!(ax2, αvals, βvals, θvals, colormap=:afmhot, colorrange=(0,π));
    cb2 = Colorbar(fig[2, 2], hm2);
    hm3 = CairoMakie.heatmap!(ax3, αvals, βvals, ϕvals, colormap=:afmhot, colorrange=(0,2π));
    cb3 = Colorbar(fig[2, 4], hm3);

    ax4 = Axis(
        fig[3, 1:3], 
        leftspinevisible=false, 
        rightspinevisible=false, 
        topspinevisible=false,
        bottomspinevisible=false,
        height=60,
    );
    hidedecorations!(ax4);
    CairoMakie.text!(ax4, L"θ_s=%$(Int(floor(θs*180/π)))^\circ", fontsize=30);
end

# ![image](emission_coordinates.gif)
