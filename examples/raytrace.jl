# # Raytracing with Minotime
#
# ### Coordinate evolution as a function of τ.

using CairoMakie
using Krang

model = Krang.Kerr(0.99);
θo = π/4;
sze = 500;
rmin = Krang.horizon(model)
rmax = 20;
ρmax = 20;
# Let us now create a figure to plot the emission coordinates on,
fig = Figure(resolution=(600, 600), fontfamily="Computer Modern", fontface="bold");
# This loop will plot the emission coordinates for each τ.
recording = record(fig, "raytrace.gif", range(0,3,length=180),framerate=15) do τ
    empty!(fig);

    αvals = range(-ρmax, ρmax, length=sze);
    βvals = range(-ρmax, ρmax, length=sze);
    tvals = -Inf.*ones(sze, sze);
    rvals = -Inf.*ones(sze, sze);
    θvals = -Inf.*ones(sze, sze);
    ϕvals = -Inf.*ones(sze, sze);

    Threads.@threads for i in 1:sze
        for j in 1:sze
            curr_ts, curr_rs, curr_θs, curr_ϕs, νr, νθ = Krang.raytrace(model, αvals[i], βvals[j], θo, τ);
            tvals[i, j] = (rmax > curr_rs > rmin) ? curr_ts : tvals[i, j];
            rvals[i, j] = (rmax > curr_rs > rmin) ? curr_rs : rvals[i, j];
            θvals[i, j] = (rmax > curr_rs > rmin) ? curr_θs : θvals[i, j];
            ϕvals[i, j] = (rmax > curr_rs > rmin) ? curr_ϕs : ϕvals[i, j];

        end
    end

    ax0 = Axis(fig[1, 1], aspect=1, title=L"\text{Regularized Time} (t_s)");
    ax1 = Axis(fig[1, 3], aspect=1, title=L"\text{Radius} (r_s)");
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
    CairoMakie.text!(ax4, L"τ=%$τ", fontsize=30);
end

# ![image](raytrace.gif)
