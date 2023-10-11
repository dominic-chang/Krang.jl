using CairoMakie
using Kang

spin = 0.99
model = Kang.Kerr(spin)
sze = 120
θo = π/4
rmax = 20
ρmax = 20
rmin = Kang.horizon(model)

fig = Figure(resolution=(600, 600), fontfamily="Computer Modern", fontface="bold");

record(fig, "emission_coordinates.mp4", range(0,π,length=180),framerate=15) do θs
    empty!(fig)

    αvals = range(-ρmax, ρmax, length=sze)
    βvals = range(-ρmax, ρmax, length=sze)
    tvals = -Inf.*ones(sze, sze)
    rvals = -Inf.*ones(sze, sze)
    θvals = -Inf.*ones(sze, sze)
    ϕvals = -Inf.*ones(sze, sze)

    n = 1
    
    Threads.@threads for i in 1:sze
        for j in 1:sze
            curr_ts, curr_rs, curr_θs, curr_ϕs, νr, νθ = Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, true, n) 
            tvals[i, j] = (rmax > curr_rs > rmin) ? νr ? 1 : -1 : tvals[i, j]
            rvals[i, j] = (rmax > curr_rs > rmin) ? νθ ? 1 : -1 : rvals[i, j]
            θvals[i, j] = (rmax > curr_rs > rmin) ? curr_θs : θvals[i, j]
            ϕvals[i, j] = (rmax > curr_rs > rmin) ? curr_ϕs : ϕvals[i, j]

            curr_ts, curr_rs, curr_θs, curr_ϕs, νr, νθ = Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, false, n)
            tvals[i, j] = (rmax > curr_rs > rmin && curr_ts != 0) ? νr ? 1 : -1 : tvals[i, j]
            rvals[i, j] = (rmax > curr_rs > rmin && curr_rs > 0) ? νθ ? 1 : -1 : rvals[i, j]
            θvals[i, j] = (rmax > curr_rs > rmin && curr_θs > 0) ? curr_θs : θvals[i, j]
            ϕvals[i, j] = (rmax > curr_rs > rmin && curr_ϕs > 0) ? curr_ϕs : ϕvals[i, j]
        end
    end

    n = 0
    Threads.@threads for i in 1:sze
    #for i in 1:sze
        for j in 1:sze
            curr_ts, curr_rs, curr_θs, curr_ϕs, νr, νθ = Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, true, n)
            tvals[i, j] = (rmax > curr_rs > rmin) ? νr ? 1 : -1 : tvals[i, j]
            rvals[i, j] = (rmax > curr_rs > rmin) ? νθ ? 1 : -1 : rvals[i, j]
            θvals[i, j] = (rmax > curr_rs > rmin) ? curr_θs : θvals[i, j]
            ϕvals[i, j] = (rmax > curr_rs > rmin) ? curr_ϕs : ϕvals[i, j]
            curr_ts, curr_rs, curr_θs, curr_ϕs, νr, νθ = Kang.emission_coordinates(model, αvals[i], βvals[j], θs, θo, false, n)
            tvals[i, j] = (rmax > curr_rs > rmin && curr_ts != 0) ? νr ? 1 : -1 : tvals[i, j]
            rvals[i, j] = (rmax > curr_rs > rmin && curr_rs > 0) ? νθ ? 1 : -1 : rvals[i, j]
            θvals[i, j] = (rmax > curr_rs > rmin && curr_θs > 0) ? curr_θs : θvals[i, j]
            ϕvals[i, j] = (rmax > curr_rs > rmin && curr_ϕs > 0) ? curr_ϕs : ϕvals[i, j]

        end
    end

    ax0 = Axis(fig[1, 1], aspect=1, title=L"\text{Radius increasing on emission}");
    ax1 = Axis(fig[1, 3], aspect=1, title=L"\text{inclination increasing on emission}");
    ax2 = Axis(fig[2, 1], aspect=1, title=L"\text{Opening angle } (\theta_s)");
    ax3 = Axis(fig[2, 3], aspect=1, title=L"\text{Azimuth } (\phi_s)");

    hm0 = CairoMakie.heatmap!(ax0, αvals, βvals, tvals, colormap=:afmhot, colorrange=(-1,1));
    cb0 = Colorbar(fig[1, 2], hm0);
    hm1 = CairoMakie.heatmap!(ax1, αvals, βvals, rvals, colormap=:afmhot, colorrange=(-1,1));
    cb1 = Colorbar(fig[1, 4], hm1);
    hm2 = CairoMakie.heatmap!(ax2, αvals, βvals, θvals, colormap=:afmhot, colorrange=(0,π))
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
    hidedecorations!(ax4)
    CairoMakie.text!(ax4, L"θ_s=%$(Int(floor(θs*180/π)))^\circ", fontsize=30) 
    display(fig)
end

